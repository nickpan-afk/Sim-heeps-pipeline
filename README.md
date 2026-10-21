# HEEPS end-to-end PSF production

For the **combined OPD/amplitude -> HEEPS -> ScopeSim workflow**, use
[`simulation_e2e.py`](simulation_e2e.py) and
[`configs/simulation_l6.yaml`](configs/simulation_l6.yaml).
See [the combined workflow guide](END_TO_END.md) for commands, frame counts,
outputs, dependencies and validation. The guide below covers HEEPS-only operation.

`heeps_e2e.py` accepts local wavefront FITS files and a YAML configuration and
produces a non-coronagraphic, on-axis intensity PSF cube. It follows the configuration,
path-resolution, dry-run, and preview conventions of the supplied ScopeSim L-band
program. It does not execute either notebook or import the ScopeSim application.

The application is self-contained: `heeps_e2e.py` contains configuration validation,
chunking, unit conversion, wavelength averaging, HEEPS calls, FITS output, and resume.
Deploy just that script and a YAML file alongside your input/calibration files.
It does not import or read `psf_workflow`, execute notebooks, or require another
application helper module. The original helper directory remains for the notebook.

Use the existing **heeps** environment (NumPy, Astropy, PyYAML, PROPER, and HEEPS
dependencies). HEEPS must be installed/importable. The unchanged official checkout
at commit `c7dc898` is supported; a custom HEEPS download switch is not required.
Within this repository the adjacent
HEEPS checkout is used automatically; elsewhere use installed HEEPS or put the
HEEPS repository's parent-of-package directory on `PYTHONPATH`, for example:

```sh
PYTHONPATH=/path/to/HEEPS python /path/to/heeps_e2e.py --config /path/to/heeps_l6.yaml --preview
```

HEEPS is imported directly inside `heeps_backend()` when propagation starts.
Keeping that import inside the function lets `--dry-run` validate inputs without
loading HEEPS/PROPER. No service or instrument-package download is required here.

## Quick start

From the HEEPS repository root:

```sh
conda activate heeps

# Validate all requested input coverage and display the execution plan; no output.
python sim_heeps_pipe/heeps_e2e.py --config sim_heeps_pipe/configs/heeps_l6.yaml --dry-run

# Compute a two-frame broadband preview (the example defaults to preview).
python sim_heeps_pipe/heeps_e2e.py --config sim_heeps_pipe/configs/heeps_l6.yaml

# Process all configured/common input frames, with automatic chunk resume.
python sim_heeps_pipe/heeps_e2e.py --config sim_heeps_pipe/configs/heeps_l6.yaml --full
```

The single `configs/heeps_l6.yaml` supports both phase-only and phase + amplitude
runs. It uses paths for the local calibration dataset as examples.
Set `pupil_file` and `mask_file` once to your actual calibration files. For each new
observation, update the input paths, OPD unit, cadence, and requested frame range.
Paths inside YAML are resolved relative to the YAML file, including `~` and defined
environment variables. CLI path overrides are resolved relative to the shell cwd.
The script also works when invoked by its absolute path from a different directory.

For a newly received phase cube without amplitude data, replace the
`amplitude_files` list in that same YAML with:

```yaml
amplitude_files: []
```

Then use the same configuration file:

```sh
python sim_heeps_pipe/heeps_e2e.py \
  --config sim_heeps_pipe/configs/heeps_l6.yaml \
  --phase-file /path/to/new_phase.fits --preview

python sim_heeps_pipe/heeps_e2e.py \
  --config sim_heeps_pipe/configs/heeps_l6.yaml \
  --phase-file /path/to/new_phase.fits --full
```

To add amplitude data, set `amplitude_files` in YAML or pass
`--amplitude-file /path/to/amplitude.fits`. Additional overrides are `--nframes`,
`--cpu-count`, and `--output-dir`. Run `--help` for details.
No separate phase-only YAML is needed. Changing the amplitude input creates a new
configuration fingerprint, keeping previous results separate.

The application returns exit code 0 on success, 2 for configuration/input/run errors,
and 130 on keyboard interruption. A failed run leaves completed chunks available
for a subsequent invocation with the same inputs and configuration.

## Inputs and settings

| Setting | Meaning |
| --- | --- |
| `pupil_file` | Required finite, non-negative 2D entrance-pupil amplitude image. |
| `mask_file` | Required finite, non-negative 2D aperture-mask/Lyot-stop amplitude image. This is the exact file used; mask selection is not inferred from its filename. |
| `phase_files` | Optional list of 3D `(frame, y, x)` real optical-path-difference cubes. At least phase or amplitude is required unless `pure: true`. |
| `phase_unit` | Explicit `m`, `um`, or `nm` of OPD, converted to meters before HEEPS. Phase in radians is rejected rather than guessed. |
| `amplitude_files` | Optional real electric-field amplitude cube, a static 2D map, or a one-frame cube. An intensity map is not automatically square-rooted. `[]` disables amplitude errors. |
| `global_start` | Global index represented by a source file's first frame; defaults to 0. Applies independently to phase and amplitude sources. |
| `frame_start` | First global frame to produce, inclusive. |
| `nframes` | Positive output frame count. `null` infers the common available range. Pure/static-amplitude-only runs require an explicit count. |
| `cadence_s` | Required source time sampling; written to FITS WCS. Not inferred from filenames and not a ScopeSim detector exposure setting. |
| `pixel_scale_mas` | `auto` or omitted selects sampling from the band profile; a number explicitly overrides it, in mas/pixel. |
| `hfov_arcsec` | Half field of view, not the full width. |
| `spectrum` | `auto` or omitted selects the band profile; explicit center/bandwidth or wavelengths/weights override it. |
| `band_specs` | Optional per-band overrides for center wavelength, pixel scale, fractional bandwidth and sample count. |
| `chunk_size`, `cpu_count` | Memory/process controls; defaults are 100 frames and one worker. |
| `run_mode` | `preview` or `full`. Preview uses a separate name, at most `preview_frames` frames, and one worker. |
| `pure` | Disable all wavefront errors, compute one static PSF, then repeat it for the explicitly requested frame count. Input error lists must be empty. |
| `input_root` | Optional HEEPS path-resolution root; defaults to the pupil file's parent. Explicit input/calibration paths do not need the legacy HEEPS directory layout. Keep output outside this input root. |

Source examples:

```yaml
# Complete input cube:
phase_files:
  - path: incoming/phase.fits
    global_start: 0

# Alternatively, several pre-sliced cubes with explicit origins:
# phase_files:
#   - path: incoming/phase_0_10000.fits
#     global_start: 0
#   - path: incoming/phase_10000_20000.fits
#     global_start: 10000
```

All ranges are zero-based and stop-exclusive. The planner splits at input-source
boundaries and verifies coverage; no source file is modified. Non-finite input
samples are rejected. For production, confirm phase/amp temporal alignment and
calibration geometry before selecting a long frame range.

The default broadband calculation is an equal-weight mean of five requested
wavelengths:

```yaml
spectrum:
  lambda_center_m: 3.81e-6
  frac_bw: 0.182
  nlam: 5
```

For monochromatic propagation, use `frac_bw: 0` and `nlam: 1`. For explicit
transmission or other user-defined weights, replace that mapping with:

```yaml
spectrum:
  wavelengths_m: [3.5e-6, 3.81e-6, 4.1e-6]
  weights: [0.2, 0.6, 0.2]
```

Weights are normalized once; individual PSFs are not renormalized by their flux.
No source spectrum, detector conversion, or photon weighting is inferred. HEEPS
quantizes propagation wavelengths to its pupil sampling; requested and actual
wavelengths are both recorded. Band/mask combinations follow the existing HEEPS
rules: L/M with LM6c/LM12c and N1/N2 with N6c/N12c. Supply the corresponding mask
file, wavelength settings, and pixel scale explicitly for other bands.

## Output and resume

The terminal prints two machine-readable lines on success:

```text
PSF_CUBE=/absolute/output/name/fingerprint/full.fits
RESULT_JSON=/absolute/output/name/result.json
```

`full.fits` is a float32 intensity cube in `(time, y, x)` order. It includes spatial
and time WCS, cadence, global frame bounds, input error flags, wavelength samples,
weights, and FITS checksums. Intensities preserve the HEEPS entrance-pupil
normalization and mask/crop throughput; there is no per-frame unit-sum normalization.
This stage models neither stellar magnitudes nor thermal background, noise, or a
physical detector readout.

The name includes `_preview` for previews. Its directory contains:

```text
output_dir/
  experiment_name[_preview]/
    result.json                 # Pointer to the most recently completed run
    fingerprint/
      manifest.json             # Input identities, effective experiment, code/runtime digest
      chunks/                   # Completed PSF chunks with checksums
      full.fits                 # Final merged PSF cube
      progress.json
```

`result.json` records the absolute cube path, requested configuration, effective
experiment, execution plan, and completion time. It is published only after a
successful merge. There is no extra copy of the potentially multi-GB final cube.
For automation:

```python
import json
from pathlib import Path

result = json.loads(Path("/path/printed/as/RESULT_JSON").read_text())
psf_cube = Path(result["psf_cube"])
```

The same validated configuration reuses completed chunks and the final cube.
Different inputs, OPD units, calibrations, scientific settings, or source code produce
new fingerprints. Invalid existing results are retained and rejected rather than
silently overwritten. Merging uses bounded-memory memmap. The full merged output
and chunk files both occupy disk. Input FITS files are read directly with memory
mapping; only the current frame interval is sampled. OPD conversion happens in
memory, and no `input_chunks/` directory or converted-input FITS files are written.
The mappings remain open while the wavelengths for that chunk are propagated.
Reduce chunk size if memory is constrained. A force-killed process can leave a lock; check its
recorded PID before removing a stale lock.

Input identity uses path/size/mtime rather than a full hash of every large source
cube. Update the mtime if another tool modifies data while preserving both size and
timestamp. Source fingerprints cover this script and the actual importable HEEPS
package; no helper-directory source is consulted. Consolidating the application
changes that fingerprint, so the first run uses a new cache directory while old
outputs remain available. Version metadata may be unavailable for source-installed
PROPER; change the experiment name after upgrading such an installation.

## Passing the result to the supplied ScopeSim program

Use `PSF_CUBE` as the ScopeSim YAML's `input_cube`, or pass its absolute path to the
ScopeSim application's existing `--input` option. For L6, `pixel_scale_mas: 5.47`
corresponds to `pixel_scale_arcsec: 0.00547` in that program. ScopeSim's `dit`, `ndit`,
throughput parameters, source magnitudes, detector effects, and normalization remain
explicit downstream settings. They are not inferred from HEEPS cadence or mask flux.
The two supplied ScopeSim files are reference files and have not been edited or run.

## HEEPS scripts used

The single application script contains the orchestration and FITS routines. Its
`heeps_backend()` imports HEEPS directly and calls these existing HEEPS functions:

| HEEPS source | Function and role |
| --- | --- |
| `heeps/config/read_config.py` | `read_config`: HEEPS optical defaults and path initialization through the offline compatibility adapter. |
| `heeps/config/update_config.py` | `update_config`: pupil/detector sampling and actual propagation wavelength. |
| `heeps/pupil/pupil.py` | `pupil`: entrance wavefront and entrance-pupil normalization. |
| `heeps/util/img_processing.py` | `resize_img`, `resize_cube`: phase/amplitude sampling; `pad_img` is used internally. |
| `heeps/wavefront/propagate_cube.py` | `propagate_cube`: frame propagation using explicit error arrays. |

Indirectly executed through these functions:

- `heeps/wavefront/propagate_one.py`: per-frame propagation.
- `heeps/wavefront/add_errors.py`: OPD phase and amplitude multiplication.
- `heeps/optics/lyot_stop.py`: loads and applies the explicitly provided aperture mask.
- `heeps/optics/detector.py` and `heeps/optics/lens.py`: intensity extraction and detector-plane cropping.
- `heeps/util/multiCPU.py`: resampling/propagation workers.
- `heeps/util/round2.py`: HEEPS odd-grid sampling.
- `heeps/optics/apodizer.py` and `heeps/optics/fp_mask.py` are invoked by propagation,
  but their coronagraph branches do not apply in the fixed `ELT` mode.

`heeps/wavefront/load_errors.py` and `heeps/wavefront/propagate.py` are used only in
reference validation, not in production: their original loading path does not
forward the CPU limit to resampling. `heeps/util/save2fits.py` is not used for final
output; explicit FITS/WCS metadata is written inside `heeps_e2e.py`.

`_local_heeps_config` uses `download_inputs=False` when that parameter is explicitly
supported. For the official legacy reader, it executes the same function with a
private copy of its global bindings and a disabled example downloader. Other
initialization and optical defaults remain those of the installed HEEPS version.
The adapter does not edit HEEPS source or replace objects in its shared modules.
The original reader remains unchanged for other callers. Unrecognized readers
without a supported offline path fail before initialization.

Explicit pupil and aperture-mask FITS files mean this workflow does not require
the local mask-generation or automatic mask-selection extensions. Copy this
folder into the official HEEPS checkout to use the adjacent upstream package.
See [the stock-HEEPS validation](../docs/heeps_upstream_validation.json) and the
[combined guide](END_TO_END.md) for the tested revision and scope.

## Validation

```sh
python -m unittest discover -s tests -v
# Optional real-input CLI validation with isolated three-frame fixtures:
python tests/check_heeps_e2e.py
```

The unit tests include copying the script to a separate directory and blocking
both imports and filesystem access to `psf_workflow`. They also test the embedded
chunk/resume/merge routines directly, including checksums and corrupt-cache rejection.

The real-input check keeps its fixtures in the operating system's temporary
directory rather than the production output root. It uses the copied script with
the same isolation guards and exercises phase-only, nanometer OPD, phase+amplitude, pure, and static amplitude
from new local fixtures without the legacy reference dataset. It checks
FITS/WCS/checksums, a different cwd, resume without rewriting valid cubes, OPD
conversion, the previous implementation's outputs when available, and an independent
original HEEPS propagation reference. See
[the recorded results](../docs/heeps_e2e_validation.json).

### Single-file consolidation checks

- All 28 unit tests passed, including standalone deployment and embedded FITS/cache tests.
- Five real-input profiles passed with imports and filesystem access to the former
  helper blocked: phase-only, nanometer OPD, phase+amplitude, pure, and static amplitude.
- A saved pre-consolidation entry point was run in the same environment against the
  same fixtures. All five resulting cubes matched the standalone outputs exactly.
- A copied script also completed the five-wavelength, two-frame L6 preview with valid
  FITS checksums while the former helper remained blocked.
- No full 36000-frame production simulation was launched. The notebook helper and
  existing HEEPS source files were retained.

For an optional same-environment comparison with a saved previous entry point, use
`python tests/check_heeps_e2e.py --baseline-entry /path/to/previous_heeps_e2e.py`.
That previous file is used only as a validation reference, not by the application.

### Direct input reads

The application reads phase/amplitude intervals directly from the original FITS
files. Global source origins, static amplitude maps, non-finite-data rejection,
and OPD conversion retain their previous behavior, with no persistent input cache.
Only PSF chunks are stored for resume; their input files must remain available.
Existing input-cache folders from historical runs are not deleted automatically.

The direct-read change passed 30 unit tests and same-environment comparisons for
five real-input profiles. All five matched the earlier input-cache implementation
pixel for pixel, with no input-cache directory in the new run folders.
