# HEEPS + ScopeSim end-to-end simulation

`simulation_e2e.py` runs local OPD/amplitude inputs through HEEPS, then passes the
resulting PSF intensity cube directly to ScopeSim to generate METIS science and
optional calibrator detector cubes. No notebook or `psf_workflow` is used.

## Run

From the HEEPS repository root, using an environment containing both HEEPS/PROPER
and ScopeSim:

```sh
conda activate heeps
python sim_heeps_pipe/simulation_e2e.py --config sim_heeps_pipe/configs/simulation_l6.yaml --dry-run
python sim_heeps_pipe/simulation_e2e.py --config sim_heeps_pipe/configs/simulation_l6.yaml --preview
python sim_heeps_pipe/simulation_e2e.py --config sim_heeps_pipe/configs/simulation_l6.yaml --full
```

Edit **`configs/simulation_l6.yaml`** for the combined workflow. The individual
HEEPS and ScopeSim YAML files are only needed when running those stages separately.
The combined example preserves phase + amplitude, five wavelengths, the Lp filter,
magnitude 7, the configured optical/noise effects, and one science/calibrator pair.
It sets an explicit noise seed and automatically downloads missing instrument
packages to `sim_heeps_pipe/inst_pkgs/`, reusing packages already there.

Update the phase/amplitude/pupil/mask paths for your machine. `scopesim.package_dir`
is both the download destination and the directory used to load packages. The example
uses `../inst_pkgs`, relative to the YAML, so no post-download path change is needed.
To reuse packages installed elsewhere, set this value to their parent directory
(containing `METIS/`, `ELT/`, and `Armazones/`). All YAML paths are relative to that YAML, including paths within sections;
`~` expands to your home directory. `--output-dir` is relative to the shell cwd.
For phase-only inputs set `heeps.amplitude_files: []`. Input phase is optical path
difference with explicit `phase_unit`, not phase in radians.

`--dry-run` validates configuration, local optical inputs, requested frame coverage,
group capacity and local package-directory presence. It prints the execution plan
without creating outputs or downloading anything. Instrument construction, effect
names, and actual PSF pixels are checked during execution. Missing package directories
are listed when downloads are enabled; otherwise they cause an error.

## Frame counts and timing

The top-level `run_mode` controls both stages. Stage-specific `run_mode` and output
paths are managed by the driver and cannot be set inside the sections. ScopeSim's
`input_cube` and input `pixel_scale_arcsec` are also managed automatically: the
driver checks the generated FITS sampling and converts HEEPS `PSCALE` from mas.
Detector sampling continues to come from the METIS instrument package.

For a star, the minimum HEEPS frame count is:

```text
nplanes × average_group × (2 if make_calibrator else 1)
```

In the **combined configuration**, `heeps.nframes: null` means compute this minimum.
It does not mean propagate the entire available phase cube. For example, 100 science
and 100 calibrator images averaging 10 PSFs each require 2,000 HEEPS frames. Set an
integer `heeps.nframes` to keep extra PSFs in full mode. The count must cover the
requested groups, with incomplete trailing groups rejected unless explicitly dropped.
The HEEPS-only entry point retains its original null = all available frames behavior.

Preview uses one detector plane, the first magnitude, and just enough HEEPS frames
for the selected grouping and calibrator setting. Full uses all requested detector
planes/magnitudes. Neither mode infers how many images you want from the input size.
Empty-sky mode needs only one PSF frame for geometry when `nframes` is null; PSF
pixel values are unused by ScopeSim in that mode.

Calibrators retain the existing midpoint split of complete PSF groups. Science uses
the first groups; calibrators use groups starting at `number_of_groups // 2`. They
are from the same optical-input sequence; this does not create an independently
configured calibration observation. Plan frame ranges are zero-based, stop excluded,
and local to the generated PSF cube. The HEEPS plan separately records global input
frame indices.

`cadence_s` is the interval between input samples; `average_group × cadence_s` is
the duration represented by a PSF group under uniform sampling. Detector exposure
is explicitly `dit × ndit`. The driver reports a warning when these differ and
does not modify them. The existing single averaged PSF per readout method remains.

## Scientific interface

### Automatic HEEPS band profiles

Both example YAML files now contain a `band_specs` table. Keep `spectrum: auto`
and `pixel_scale_mas: auto` to select the row matching `heeps.band` (or `band` in
HEEPS-only YAML). The same defaults apply when these two controls are omitted.

| Band | Center wavelength (microns) | PSF sampling (mas/pixel) | Default spectrum |
| --- | --- | --- | --- |
| L | 3.81 | 5.47 | Five equally weighted samples, fractional bandwidth 0.182 |
| M | 4.79 | 5.47 | Single wavelength |
| N1 | 8.70 | 6.79 | Single wavelength |
| N2 | 11.33 | 6.79 | Single wavelength |

Centers and sampling follow the existing HEEPS `read_config` band definitions,
using the configured GeoSnap N-band sampling. Only the existing L broadband settings
are preserved; no bandwidth has been inferred for another filter. For broadband
operation edit the selected band's `frac_bw` and `nlam` in the YAML table. For example:

```yaml
heeps:
  band: N2
  mask_name: N12c
  mask_file: ../inputs/calibration/mask/ls_SPM_N12c.fits
  spectrum: auto
  pixel_scale_mas: auto
  band_specs:
    N2:
      lambda_center_m: 11.33e-6
      pixel_scale_mas: 6.79
      frac_bw: 0.0             # Set the intended bandwidth for broadband work
      nlam: 1
```

This is a partial example; keep the other optical inputs and settings in your YAML.
Per-band entries can override any subset of these four values. Unknown band/profile
keys and invalid values are rejected, including entries for unselected bands.
An explicit `spectrum` mapping or numeric `pixel_scale_mas` takes precedence over
the band table, preserving existing configurations and custom wavelength weights.

N1/N2 accept N6c/N12c; L/M accept LM6c/LM12c. Change the mask file and
appropriate wavefront inputs when switching families. The actual HEEPS wavelength
may differ slightly from the requested center due to pupil-grid rounding and is
recorded in the FITS header. HEEPS's existing `flux_star`/`flux_bckg` definitions
remain in `read_config`; the profile table does not redefine absolute photometry.

The joint driver automatically transfers the resolved PSF sampling to ScopeSim.
It does **not** change ScopeSim mode, filter, detector settings or throughput when
the HEEPS band changes; those settings must be configured for the chosen instrument
setup separately.

All four bands with both compatible masks passed one-frame monochromatic, phase-only
HEEPS propagation checks (8 cases); see [the recorded result](../docs/heeps_band_validation.json).
This does not validate N-band ScopeSim integration, non-L broadband weighting, or
the physical suitability of a particular input dataset for another wavelength.

### PSF transfer to ScopeSim

HEEPS produces `(time, y, x)` relative intensity, with a weighted mean across its
wavelength samples. ScopeSim averages temporal groups and normalizes them to
`flux_scale`, then sets source photometry using its Vega spectrum and magnitude.
Absolute HEEPS throughput is therefore not automatically carried into photometry.
ScopeSim's own PSF effect must remain disabled and is rejected if enabled here.

This is a broadband spatial PSF passed into ScopeSim. It does not supply a separate
HEEPS PSF at every ScopeSim wavelength. Match HEEPS spectral weights/band to your
scientific setup; the driver cannot establish spectral equivalence from the filter
name alone. Both stages otherwise retain their existing propagation and detector
models, including the cold-stop transmission and configured effects.

## Outputs and retry

```text
outputs/<name>_<preview-or-full>_<configuration_hash>/
  psf_cube.fits
  science_mag7.fits
  calibrator_mag7.fits           # If requested
  config_input.yaml
  config_effective.yaml
  manifest.json
  .work/                        # Hidden internal data; retain for resume
    heeps/                      # PSF chunks, progress, manifest and result pointer
    scopesim_run01/              # Detector configuration and provenance
```

The top-level manifest lists absolute PSF and detector-product paths and records
each attempt and its failure stage. `config_input.yaml` is the source YAML;
`config_effective.yaml` is a resolved record of the two stage configurations,
not a replacement for the input schema. ScopeSim also writes its own resolved,
rerunnable stage configuration and instrument/dependency provenance.

The PSF is written directly as `psf_cube.fits` and passed by path without copying it.
No `input_chunks/` is generated. HEEPS retains its chunks and resume metadata in
`.work/heeps/`. On failure rerun the same command:
valid HEEPS chunks/full cube are checked and reused, while ScopeSim starts a fresh
run. Detector products from the second attempt have suffix `_run02`, then `_run03`,
and so on; existing products are never overwritten. Each attempt's internal records
are in `.work/scopesim_runNN/`; the top-level manifest lists the current products
and the stage-manifest paths for successful attempts. ScopeSim has no per-frame
resume. A successful rerun also starts a new detector run; with `seed: null` it
produces fresh noise. Changes to configuration or inputs can select a new run
directory and recompute HEEPS. A lock prevents simultaneous runs of the same plan.
Standalone stage entry points retain their original output layouts. Historical
nested results are left in place; only new joint runs use this flatter layout.

## Dependencies and sharing

Share this folder's three Python files (`simulation_e2e.py`, `heeps_e2e.py`, and
`scopesim_e2e_Lband.py`), preserving their directory layout, plus
the combined YAML, optical inputs/calibrations, and compatible instrument packages.
The scripts work when the folder is copied outside this repository provided the
compatible HEEPS checkout is importable, for example through `PYTHONPATH`.

The unchanged official HEEPS checkout at commit `c7dc898` is supported. This
workflow reads pupil/mask FITS explicitly and does not need the local mask-generation
extensions or a custom `read_config(download_inputs=False)` option. Its initialization
adapter uses that switch if available; otherwise it invokes the upstream reader with
private global bindings that disable only its legacy example downloader. Installed
HEEPS files, shared module bindings, and other HEEPS callers remain unchanged.
The adapter preserves the optical defaults of the HEEPS revision being used and
rejects unknown readers without a recognized offline path.

To publish on top of upstream HEEPS, place `sim_heeps_pipe/` next to upstream's
`heeps/` directory. The application selects that adjacent checkout automatically.
Copying into a different repository can change the HEEPS source fingerprints and
therefore select new run directories.

Existing ScopeSim requirements files describe
the separately tested 0.11.2 and 0.12.0 environments; they are not complete combined
HEEPS dependency specifications. The joint numerical check used the existing local
environment with ScopeSim 0.12.0, NumPy 2.5.4 and Astropy 8.0.2. A fresh combined
environment has not been independently installed or validated.

The example uses `download_missing_packages: true` and `package_release: stable`.
Missing METIS, ELT and Armazones packages are downloaded into `package_dir`; existing
packages are reused without automatic updates or replacement. Set
`download_missing_packages: false` for an entirely preinstalled package set. Normal
atmosphere construction may require the ESO SkyCalc service or a compatible cached
spectrum. The real integration test uses a cached ESO spectrum through ScopeSim's
native local-file reader in a temporary copy of Armazones; production packages and
configuration are unchanged.

Direct HEEPS calls remain in `heeps_e2e.py`:

| Function | Source |
| --- | --- |
| `read_config`, `update_config` | `heeps/config/read_config.py`, `heeps/config/update_config.py` |
| `pupil` | `heeps/pupil/pupil.py` |
| `propagate_cube` | `heeps/wavefront/propagate_cube.py` |
| `resize_img`, `resize_cube` | `heeps/util/img_processing.py` |

These functions call HEEPS's internal optical routines in the usual way. The driver
does not invoke the old notebooks or `heeps.wavefront.load_errors`; input selection,
in-memory unit conversion and bounded resampling remain in the standalone HEEPS stage.

## Validation

```sh
python -m unittest discover -s tests -v
python tests/check_simulation_e2e.py \
  --package-dir /path/to/existing/inst_pkgs \
  --offline-skycalc /path/to/cached_eso_skycalc.fits \
  --report /tmp/simulation_validation.json
```

The real checker compares joint execution with separate HEEPS and ScopeSim calls
using the same fixed seed. It compares the PSF, science and calibrator cubes pixel
for pixel, verifies HEEPS resume, and checks that no `input_chunks` exists. Products
and copied instrument packages are temporary and removed after the check.
See [the recorded result](../docs/simulation_e2e_validation.json). This establishes
equivalence to manual execution in the tested environment, not an independent
physical calibration or validation of live ESO retrieval.

The [upstream compatibility report](../docs/heeps_upstream_validation.json) covers
eight one-frame monochromatic phase-only band/mask cases and an L-band five-wavelength,
two-frame phase+amplitude HEEPS -> ScopeSim comparison. The latter uses the existing
instrument packages and a cached atmosphere spectrum in temporary package copies.
No legacy HEEPS example downloads occur; the external HEEPS checkout remains clean.
This verifies the recorded revision, not every historical or future HEEPS release.
Reproduce it with:

```sh
python tests/check_upstream_heeps.py \
  --heeps-root /path/to/unchanged/HEEPS \
  --offline-skycalc /path/to/cached_eso_skycalc.fits \
  --report /tmp/heeps_upstream_validation.json
```
