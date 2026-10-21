# sim_heeps_pipe — End-to-end simulation pipeline

`sim_heeps_pipe` is the complete command-line workflow for converting optical path
difference (OPD) and optional amplitude FITS inputs into METIS detector simulations.
It runs HEEPS propagation, passes the generated PSF cube to ScopeSim automatically,
and writes the PSF, science and optional calibrator cubes in one result directory.

```text
OPD + optional amplitude + pupil + aperture mask
                      |
                    HEEPS
                      |
               PSF intensity cube
                      |
                   ScopeSim
                      |
          Science / calibrator detector cubes
```

**Entry point:** [simulation_e2e.py](simulation_e2e.py)

**Configuration:** [configs/simulation_l6.yaml](configs/simulation_l6.yaml)

Edit this one YAML and run the entry point. It handles input validation, HEEPS
propagation, PSF handoff, detector simulation and result bookkeeping. The two stage
scripts are runtime components called by the driver; users normally do not run
them separately. No notebook, `psf_workflow` folder or test suite is required.

## 1. Repository and data layout

The application folder belongs next to the official HEEPS source directory:

```text
Sim-heeps-pipeline/                  # Repository root
├── heeps/                          # Official HEEPS source
└── sim_heeps_pipe/
    ├── README.md                   # This guide
    ├── END_TO_END.md               # Additional technical details
    ├── simulation_e2e.py           # Run this entry point
    ├── heeps_e2e.py                # HEEPS stage used by the driver
    ├── scopesim_e2e_Lband.py        # ScopeSim stage used by the driver
    ├── configs/
    │   └── simulation_l6.yaml      # The configuration to edit
    ├── inputs/
    │   ├── calibration/
    │   │   ├── pupil/
    │   │   │   └── ELT_fullM1_feb2025_noseg.fits
    │   │   └── mask/
    │   │       ├── ls_SPM_LM6c.fits
    │   │       ├── ls_SPM_LM12c.fits
    │   │       ├── ls_SPM_N6c.fits
    │   │       └── ls_SPM_N12c.fits
    │   └── wavefront/
    │       ├── phase.fits
    │       └── amplitude.fits      # Optional
    ├── inst_pkgs/                  # Created automatically when packages download
    └── outputs/                    # Created by simulations
```

Only the mask needed for a simulation is required. FITS files are supplied
separately; cloning the repository does not provide your observation data or
calibration files. Filenames may differ from the example: update the YAML paths.

Paths inside YAML are relative to **the YAML file**, not the terminal directory.
For `configs/simulation_l6.yaml`, `../inputs/` means `sim_heeps_pipe/inputs/`.
Absolute paths and `~` are also supported. Command-line path overrides are relative
to the terminal working directory.

## 2. Environment setup

Use the following HEEPS reference versions, with ScopeSim 0.11.2 for the
pipeline's legacy ScopeSim interface:

| Package | Version |
| --- | --- |
| Python | 3.12.7 |
| PROPER | 3.3.3 |
| NumPy | 1.26.4 |
| SciPy | 1.15.0 |
| Matplotlib | 3.10.6 |
| VIP-HCI | 1.6.2 |
| Photutils | 2.1.0 |
| Astropy | 6.1.7 |
| scikit-image | 0.24.0 |
| ScopeSim | 0.11.2 |

The scientific-stack versions follow the original HEEPS reference environment.
ScopeSim 0.11.2 has been checked with this pipeline; its previous validation used
SciPy 1.17.0. The exact SciPy 1.15.0 combination listed here still needs a preview
on the recipient's machine. ScopeSim 0.12.0 requires newer NumPy/Astropy versions
and should not be substituted into this older scientific stack.

### Create the environment

```bash
git clone https://github.com/nickpan-afk/Sim-heeps-pipeline.git
cd Sim-heeps-pipeline

conda create -n heeps-sim python=3.12.7 pip
conda activate heeps-sim

python -m pip install "pip==25.2" "setuptools==80.9.0" wheel
python -m pip install \
  "numpy==1.26.4" \
  "scipy==1.15.0" \
  "matplotlib==3.10.6" \
  "vip_hci==1.6.2" \
  "photutils==2.1.0" \
  "astropy==6.1.7" \
  "scikit-image==0.24.0" \
  "scopesim==0.11.2" \
  "synphot==1.6.0" \
  "PyYAML==6.0.3" \
  requests beautifulsoup4
```

### Install Python PROPER

Download Python PROPER 3.3.3 from the
[official PROPER files page](https://sourceforge.net/projects/proper-library/files/),
extract it and install it in the same environment:

```bash
# Replace this path with your extracted PROPER directory.
python -m pip install /path/to/proper_v3.3.3_python
```

Follow PROPER's supplied installation instructions if a platform-specific C
compiler is required. The Python import is `proper`.

### Verify imports

From the repository root:

```bash
python -m pip check
python -c "import heeps, proper, scopesim; print('HEEPS:', heeps.__file__); print('ScopeSim:', scopesim.__version__)"
```

HEEPS should resolve to the repository's adjacent `heeps/` directory. The unchanged
official checkout at commit `c7dc898` is supported by the initialization adapter;
no local mask-generator extensions are needed when using explicit FITS masks.

For ScopeSim 0.11.2, obtain the matching instrument-package set from the project
maintainer. Use `download_missing_packages: false` once that set is installed.
The current stable METIS package targets ScopeSim 0.12, so a fresh stable download
is not the recommended package source for this legacy environment.

## 3. Configure an observation

Edit **`configs/simulation_L.yaml`** inside `sim_heeps_pipe/`. Start by setting
the optical input paths, OPD units, cadence, band and matching mask. This file
contains both the `heeps` and `scopesim` sections for the complete run.

This is a complete compact example for L-band, phase-only operation:

```yaml
name: example_L_LM6c
output_dir: ../outputs
run_mode: preview

heeps:
  pupil_file: ../inputs/calibration/pupil/ELT_fullM1_feb2025_noseg.fits
  mask_file: ../inputs/calibration/mask/ls_SPM_LM6c.fits
  phase_files:
    - path: ../inputs/wavefront/phase.fits
      global_start: 0
  phase_unit: m
  amplitude_files: []
  pure: false
  band: L
  mask_name: LM6c
  spectrum: auto
  pixel_scale_mas: auto
  frame_start: 0
  nframes: null
  cadence_s: 0.1
  hfov_arcsec: 1.1
  ngrid: 1024
  chunk_size: 100
  cpu_count: 1

scopesim:
  package_dir: ../inst_pkgs
  download_missing_packages: false
  package_release: stable
  source_kind: star
  mags: [7.0]
  nplanes: 10
  average_group: 1
  make_calibrator: true
  dit: 0.1
  ndit: 1
  mode: img_lm
  filter_name: Lp
  cold_stop_transmission: 0.07511
  flux_scale: 1.0
  seed: 20261021
  show_effects: true
  effects:
    psf: false
```

Omitted effect switches use the application's defaults, including the standard
detector/noise chain. Keep `psf: false` because HEEPS already supplies the PSF.
The example cold-stop transmission belongs to the supplied L6 setup; review it
when changing the mask or optical configuration.

### Input requirements

| Input | Expected content |
| --- | --- |
| Pupil | Finite, non-negative 2D real electric-field amplitude transmission. |
| Mask | Finite, non-negative 2D real aperture-mask/Lyot-stop amplitude transmission. |
| Phase | A 3D `(frame, y, x)` **OPD** cube, with explicit `phase_unit: m`, `um` or `nm`. Phase angles in radians are not accepted as OPD. |
| Amplitude | An optional real electric-field amplitude cube, static 2D map or one-frame cube. Intensity is not automatically converted to amplitude. |

Input grids must describe consistent pupil geometry and aligned time samples.
The program resamples images but cannot infer missing physical-coordinate
information or establish whether an input dataset is appropriate for a new band.

For phase + amplitude operation, replace the empty amplitude list with:

```yaml
  amplitude_files:
    - path: ../inputs/wavefront/amplitude.fits
      global_start: 0
```

For multiple input cubes, list each path and the global index of its first frame
using `global_start`. `frame_start` is the first requested global index. Indices
are zero-based; reported stop indices are excluded. Input files remain unchanged.

## 4. Run the complete pipeline

Run the following commands from **inside `sim_heeps_pipe/`**:

```bash
conda activate heeps-sim
cd /path/to/Sim-heeps-pipeline/sim_heeps_pipe

# Validate the complete execution plan; no downloads or simulation.
python simulation_e2e.py --config configs/simulation_l6.yaml --dry-run

# First test: one science plane and an optional calibrator plane.
python simulation_e2e.py --config configs/simulation_l6.yaml --preview

# Production: all configured detector planes and magnitudes.
python simulation_e2e.py --config configs/simulation_l6.yaml --full
```

Replace `/path/to/Sim-heeps-pipeline` with your checkout location. If you already
have a working environment named `heeps`, activate that environment instead of
`heeps-sim`.

A normal run performs these steps automatically:

1. Validate optical input files, frame coverage and detector frame requirements.
2. Compute or reuse valid HEEPS PSF chunks, then assemble `psf_cube.fits`.
3. Locate or download the required ScopeSim instrument packages.
4. Pass the PSF cube and its spatial sampling to ScopeSim.
5. Save science/calibrator detector cubes and a top-level `manifest.json`.

To run from the repository root instead, keep the application folder in both paths:

```bash
python sim_heeps_pipe/simulation_e2e.py \
  --config sim_heeps_pipe/configs/simulation_l6.yaml --preview
```

Override the output root from inside `sim_heeps_pipe/` without editing YAML:

```bash
python simulation_e2e.py --config configs/simulation_l6.yaml \
  --output-dir /path/to/results --full
```

The entry point accepts `--config`, `--output-dir`, `--dry-run`, `--preview` and
`--full`. Set scientific controls such as `nplanes`, `seed`, `band` and input paths
in YAML. Display the available command-line options with:

```bash
python simulation_e2e.py --help
```

`--preview` and `--full` override the top-level YAML `run_mode`. Without either
flag, the YAML value controls the run. Preview uses the first magnitude, one
science plane and an optional calibrator plane, with enough HEEPS frames for the
configured `average_group`.

`--dry-run` checks configuration, local input geometry/coverage and frame capacity.
It reports missing instrument-package directories but does not build METIS,
validate effect names against an instrument or validate every PSF pixel. Those
checks occur during execution. With downloads disabled, missing packages cause
an error rather than being scheduled for download.

## 5. Adjust HEEPS parameters

| Parameter | Meaning and practical use |
| --- | --- |
| `band` | Selects L, M, N1 or N2. |
| `mask_name`, `mask_file` | Set both the compatible mask label and the actual FITS file. The file is used directly; renaming the label does not change its contents. |
| `spectrum: auto` | Uses the selected band's spectral profile. An explicit spectrum mapping overrides it. |
| `pixel_scale_mas: auto` | Uses the selected band's PSF sampling; a numeric value overrides it. |
| `cadence_s` | Input time-sample interval, in seconds. Independent of ScopeSim DIT/NDIT. |
| `frame_start` | First input frame to use, in global indexing. |
| `nframes` | In the combined pipeline, `null` computes the minimum frames needed by ScopeSim. An integer requests an explicit count in full mode. |
| `chunk_size` | Number of PSF time frames computed per batch, not wavelength count or a temporal average. Smaller batches usually need less memory but create more chunk files. |
| `cpu_count` | Requested HEEPS worker limit. Start with 1 and increase after a successful preview. |
| `ngrid` | Propagation grid size. Larger grids increase memory and computation; incompatible pupil/detector dimensions are rejected. |
| `hfov_arcsec` | **Half** the PSF field of view. A value of 1.1 gives approximately a 2.2-arcsec full width. |
| `pure` | `true` computes a static no-error PSF and repeats it. Both error input lists must then be empty. Repeated frames are not independent wavefront realizations. |

Every frame completes all configured wavelength propagations before wavelength
weighting. For example, `chunk_size: 100` and `nlam: 5` still produce 100 temporal
PSF frames, not 500 output frames. Try a smaller chunk size if memory is limited.

### Per-band spectra

Both example YAML files include a `band_specs` table that can override the built-in
defaults. With `spectrum: auto` and `pixel_scale_mas: auto`, changing `band` selects
that band's row:

| Band | Center wavelength (microns) | PSF sampling (mas/pixel) | Default spectral sampling | Compatible masks |
| --- | --- | --- | --- | --- |
| L | 3.81 | 5.47 | Five equal-weight samples, fractional bandwidth 0.182 | LM6c, LM12c |
| M | 4.79 | 5.47 | Single wavelength | LM6c, LM12c |
| N1 | 8.70 | 6.79 | Single wavelength | N6c, N12c |
| N2 | 11.33 | 6.79 | Single wavelength | N6c, N12c |

Only L has an existing broadband prescription here. For another band, specify the
intended bandwidth and sample count instead of reusing L's bandwidth. For example,
this profile retains monochromatic N2 until you supply your intended broadband values:

```yaml
  band_specs:
    N2:
      lambda_center_m: 11.33e-6
      pixel_scale_mas: 6.79
      frac_bw: 0.0
      nlam: 1
```

These center/scale defaults follow the supported HEEPS configuration. Actual
propagation wavelengths may shift slightly because HEEPS rounds the pupil sampling;
requested and actual wavelengths are recorded in FITS headers.

For explicit wavelengths and weights, replace `spectrum: auto` with:

```yaml
  spectrum:
    wavelengths_m: [3.5e-6, 3.81e-6, 4.1e-6]
    weights: [0.2, 0.6, 0.2]
```

Weights are normalized automatically. Changing the HEEPS band does **not** select
the corresponding ScopeSim mode, filter, detector or throughput. Configure those
independently. The supplied ScopeSim example is L-band; N-band HEEPS propagation
checks do not establish N-band ScopeSim compatibility.

## 6. Adjust ScopeSim parameters

| Parameter | Meaning |
| --- | --- |
| `mags` | List of Vega-spectrum magnitudes, e.g. `[7.0]` or `[5.0, 7.0]`. Each magnitude gets its own detector products in full mode. |
| `nplanes` | Detector frames per magnitude and per target kind. Preview overrides this to 1. |
| `average_group` | Consecutive HEEPS intensity frames averaged to construct each input PSF. |
| `make_calibrator` | Adds a matching-magnitude calibrator cube using PSF groups from the second half. An optional companion is science-only. |
| `dit` | Seconds per detector integration. |
| `ndit` | Integrations per saved detector plane; exposure duration is `dit * ndit`. |
| `exposure_output_mode` | `average` or `sum` over NDIT through ScopeSim. These modes change output scaling. |
| `seed` | A fixed integer makes random effects repeatable under matching configuration, data, software and instrument packages. `null` leaves noise unfixed between runs. |
| `source_kind` | `star` or `empty_sky`. Empty sky has zero stellar flux but may retain atmosphere/optics/background emission. |
| `mode`, `filter_name` | ScopeSim instrument mode and filter. The supplied L setup uses `img_lm` and `Lp`. |
| `cold_stop_transmission` | Configured cold-stop transmission; review for the chosen mask/setup. |
| `flux_scale` | Sum of normalized PSF weights, normally 1. This is not an automatic measurement of HEEPS absolute throughput. |
| `companion_mag` | Optional science-only companion magnitude; separation and position angle have separate controls. |
| `effects` | Mapping of available effect names to `true`/`false`. Keep `psf: false` and the exposure-integration/output effects enabled. |

A fixed seed does not give every frame identical noise: the script derives separate
seeds for individual readouts/effects. Change the integer to generate another
repeatable noise realization.

For a noiseless linear comparison, change these effect switches while retaining
the same optical chain:

```yaml
  effects:
    psf: false
    shot_noise: false
    readout_noise: false
    dark_current: false
    detector_linearity: false
    ipc: false
```

Enabled atmosphere and optical components can still contribute mean background.
Disabling an optical effect can remove its transmission as well as its emission.

### Frame counts and exposure timing

For a star, the minimum number of HEEPS frames is:

```text
nplanes * average_group * (2 if make_calibrator else 1)
```

With the compact example above, full mode produces 10 science frames and 10
calibrator frames using 20 HEEPS frames. Preview produces one science/calibrator
pair using two HEEPS frames. To average 10 PSFs per output image, set
`average_group: 10`; the full example would then need 200 HEEPS frames.

`heeps.nframes: null` calculates this requirement; it does not process the entire
input cube. An explicit larger HEEPS frame count retains extra PSFs and keeps the
existing calibrator midpoint rule: calibrator groups start at
`number_of_complete_groups // 2`. Incomplete trailing groups are rejected unless
you explicitly select `remainder_policy: drop`.

The PSF-group duration is `average_group * cadence_s`; the detector exposure is
`dit * ndit`. They are independent settings. The driver reports a warning when
they differ and does not change either automatically. PSF averaging does not
perform a time-resolved detector integration within NDIT.

The combined driver manages stage output paths, ScopeSim `input_cube`, and input
PSF `pixel_scale_arcsec` automatically. Do not duplicate those managed settings
inside the stage sections.

## 7. Instrument packages and first-run downloads

```yaml
scopesim:
  package_dir: ../inst_pkgs
  download_missing_packages: true
  package_release: stable
```

`package_dir` is both the download destination and the later lookup directory.
On the first simulation, missing METIS, ELT and Armazones packages are downloaded
there. No post-download YAML edit is needed. Existing package directories are
reused without automatic upgrades or replacement. `package_release` only controls
new downloads; it does not update installed packages or the ScopeSim Python package.

For the ScopeSim 0.11.2 environment above, install its matching package bundle
and set `download_missing_packages: false`. The automatic-download settings shown
here describe the application capability; use a package release compatible with
your ScopeSim version. Current stable packages target ScopeSim 0.12.

First-time downloads need access to the ScopeSim instrument-package service.
Atmosphere initialization may separately require ESO SkyCalc or an appropriate
cached/local spectrum. `--dry-run` never downloads either resource.

### HTTP 503, timeout or server unavailable

A message such as `503 Service Unavailable` from
`https://scopesim.univie.ac.at/InstPkgSvr/` means the instrument-package request
failed. The current application reports the error and exits; it does not perform
automatic retries. Rerun after connectivity/service recovery.

To test while the service is unavailable, obtain a compatible, previously validated
instrument-package bundle from the project maintainer and arrange it as:

```text
sim_heeps_pipe/inst_pkgs/
├── METIS/
├── ELT/
└── Armazones/
```

For example, from the repository root:

```bash
rsync -av /path/to/validated/inst_pkgs/ sim_heeps_pipe/inst_pkgs/
```

Then rerun the same command. You may keep `download_missing_packages: true` to
reuse packages already present, or set it to `false` to require preinstalled
packages. Setting it to `false` by itself does not supply missing packages.
Use the same instrument-package versions when comparing runs or sharing a
reproducible example; the moving `stable` release can change over time.

## 8. Outputs and restarting a run

Combined runs use a flat directory for the main products:

```text
sim_heeps_pipe/outputs/
└── example_L_LM6c_preview_<hash>/
    ├── psf_cube.fits
    ├── science_mag7.fits
    ├── calibrator_mag7.fits         # If enabled
    ├── config_input.yaml
    ├── config_effective.yaml
    ├── manifest.json
    └── .work/                      # Hidden internal state; retain for resume
        ├── heeps/                  # PSF chunks and progress records
        └── scopesim_run01/         # Detector-stage configuration and provenance
```

| Product | Meaning |
| --- | --- |
| `psf_cube.fits` | HEEPS relative intensity, float32, `(time, y, x)`. Wavelength samples are already combined, not stored along the first axis. |
| `science_mag*.fits` | ScopeSim science detector cube, `(nplanes, y, x)`. Read `BUNIT` for its units. |
| `calibrator_mag*.fits` | Optional calibrator detector cube with the same shape convention. |
| `manifest.json` | Status, attempts and full paths to the PSF and detector products. |
| `config_input.yaml` | Source YAML used for the attempt. |
| `config_effective.yaml` | Resolved stage settings for inspection; this snapshot is not a replacement for the combined input schema. |

Even a single detector image is saved as a one-plane cube. For a 2D image:

```python
from astropy.io import fits

image = fits.getdata("/path/to/science_mag7.fits")[0]
```

For the same plan, rerunning checks and reuses valid HEEPS chunks and its completed
PSF cube. ScopeSim starts a new attempt rather than resuming individual detector
frames. Later detector products receive `_run02`, `_run03`, etc.; existing products
are retained. A changed input/configuration/code fingerprint can create a new run
directory and require recomputation. Old output directories are not migrated.

No `input_chunks/` is generated: input FITS slices and OPD unit conversion are
handled in memory. PSF chunks and the merged PSF still both occupy disk for resume.
Keep `.work/` and original input files until you no longer need continuation.

HEEPS PSF amplitudes are relative. ScopeSim averages and normalizes each temporal
PSF group, then assigns absolute source photometry from the configured spectrum
and magnitude. HEEPS throughput changes are not automatically inherited as an
absolute flux loss. This is a broadband spatial-PSF transfer, not a separate HEEPS
PSF at every ScopeSim wavelength.

## 9. Share the complete application

Keep `simulation_e2e.py`, `heeps_e2e.py` and `scopesim_e2e_Lband.py` together in
`sim_heeps_pipe/`, with `configs/simulation_l6.yaml` and this README. Place the folder
next to the official `heeps/` source as shown in the directory example. The driver
needs both stage files at runtime.

Tests and validation reports do not need to be uploaded for normal operation.
Distribute input/calibration FITS and, when needed, a validated instrument-package
bundle separately. Exclude observation data, `inst_pkgs/`, `outputs/` and Python
caches from GitHub. Record the environment and instrument-package versions used
for scientific comparisons.

A recipient should install the environment, supply the FITS inputs, edit the single
YAML, run `--dry-run`, then run `--preview` before starting `--full`.

The HEEPS stage directly calls `read_config`, `update_config`, `pupil`,
`propagate_cube`, `resize_img` and `resize_cube`. Their internal HEEPS propagation
routines remain in use. The application supplies bounded in-memory error arrays
instead of using `heeps.wavefront.load_errors` in production.

The official checkout needs no custom aperture-generator extensions because the
configured pupil and aperture mask are read from FITS. The initialization adapter
disables only HEEPS's unrelated legacy example download and leaves its installed
source and shared module bindings unchanged.

Additional technical details are in [END_TO_END.md](END_TO_END.md), if distributed
with the application. The commands and configuration guidance above are usable
without the repository's optional tests or validation-report files.
