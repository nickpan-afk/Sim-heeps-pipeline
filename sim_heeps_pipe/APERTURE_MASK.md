# Independent aperture-mask generation

Use the same Python environment as the simulation (HEEPS, PROPER, NumPy,
Astropy and PyYAML). No notebooks, wavefront files, instrument packages or
ScopeSim downloads are needed. This script is independent of `simulation_e2e.py`.

From `sim_heeps_pipe/`:

```sh
conda activate heeps
python create_aperture.py --config configs/aperture_mask.yaml --dry-run
python create_aperture.py --config configs/aperture_mask.yaml
```

The default example creates `inputs/calibration/mask/ls_SPM_LM6c_generated.fits`.
Existing outputs are protected; use another `output_file`, or pass `--overwrite`
explicitly. `--output-file /path/to/new_mask.fits` overrides the YAML output path.
YAML paths are relative to the YAML; CLI paths are relative to the current directory.
The dry run validates configuration only; pixel coverage is checked during generation.

## Configure the geometry

- `hole_centers_file`: path to a whitespace-separated TXT with exactly two columns,
  `x_mm y_mm`, in mask-plane millimeters. Paths are relative to the YAML file.
  Blank lines and `#` comments are supported. No radius or hole-index column is used.
  The origin is the pupil center; x increases along columns and y along rows
  (use `origin='lower'` when plotting). The number of data rows determines the hole count.
- `hole_radius_mm`: the common radius of every hole, in mask-plane millimeters.
- `mask_family`: `LM` or `N`. The fixed factors are 0.9686 and 0.9628 respectively;
  they apply to both centers and radius. New hole coordinates do not trigger
  recalculation. The existing `44.8/45` and `856.49/1000` conversions are preserved.
- `npupil`: optional, default 1385. It controls the square output sampling, not
  the physical geometry. Keep it consistent with your simulation's calibration
  sampling; resampling can change the represented edges of small holes.
- `name`: optional descriptive FITS metadata. It does not select geometry or a band.
- `base_stop`: explicit HEEPS geometry for the outer boundary, central obscuration
  and spiders. The supplied configuration reproduces the LM6c setup. Segments are
  disabled, matching the original notebook's mask creation.

For example, the YAML points to the independent coordinate input:

```yaml
hole_centers_file: ../inputs/calibration/mask_coordinates/LM6c.txt
hole_radius_mm: 2.355
```

The TXT contains:

```text
# x_mm y_mm
-4.82 -18.70
7.44 -13.21
18.60 5.21
7.74 13.12
-13.71 13.66
-15.09 0.21
```

Keep the original mask-plane coordinates in the TXT; do not pre-apply the fixed
scales. A future SAMpy integration can read the same coordinate file and apply
its own required units and conventions. This change does not add SAMpy execution.

The old notebook variants use these values in addition to their coordinate lists:

| Variant | Family | Radius (mm) | Base diam_ext (m) | Base diam_int (m) |
| --- | --- | --- | --- | --- |
| LM6c | LM | 2.355 | 36.067 | 10.8602932 |
| LM12c | LM | 0.865 | 36.067 | 11.213 |
| N6c | N | 2.200 | 35.844 | 11.213 |
| N12c | N | 0.985 | 35.844 | 11.213 |

Switching family alone does not change radii, coordinates or base-stop geometry.
Enter the geometry settings in YAML and the centers in TXT. The example TXT contains
the six-hole coordinates; for a twelve-hole mask use a TXT with twelve coordinate rows.

## Output and simulation use

One FITS file is produced: the primary image is a float32 amplitude-transmission
mask, computed as the base stop multiplied by the union of the circular holes.
There is no flux normalization or squaring. Fractional base-stop edges are retained;
circles use the original binary pixel-center criterion. Overlapping holes form a
union; holes partly behind the base stop are clipped. A fully blocked or unsampled
hole causes an error so missing apertures cannot go unnoticed.

FITS headers record the fixed factors, base-stop parameters, YAML hash and coordinate-TXT hash. A small
`HOLES` table inside the same FITS records the input centers and common radius.
HEEPS reads the primary image as usual. No sidecar files or plots are created.

Set the simulation configuration's HEEPS `mask_file` to the generated FITS manually.
Keep the simulation's supported band/mask selector consistent with its intended
setup; this utility does not change that selector or extend its accepted names.
Generating a mask never starts the simulation or changes its configuration.

## HEEPS code used

The script imports `heeps.pupil.create_stop` to generate the base stop. That HEEPS
routine uses PROPER, its polygon/hexagon helpers and `heeps.util.img_processing.resize_img`.
The circular-hole calculation is implemented in this standalone script, preserving
the old `heeps/pupil/aperture_mask.py` mathematics without importing that custom file.
An unmodified official HEEPS checkout can therefore supply the base-stop routine.
The official routine also applies a polygon outer-diameter correction that was
disabled in the local notebook version. When that correction is present, the
script compensates its input diameter with `sqrt(diam_ext**2 + seg_width**2)` so
the resulting polygon retains the YAML's requested diameter. HEEPS source files
are not changed. This adapter preserves the existing geometry across these two
implementations.
