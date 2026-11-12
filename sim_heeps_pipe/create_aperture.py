#!/usr/bin/env python3
"""Generate an aperture-mask FITS from YAML settings and a two-column coordinate TXT."""
from __future__ import annotations

import argparse
from hashlib import sha256
import inspect
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from astropy.io import fits
import yaml

REF_SCALE = 44.8 / 45.0
M_PER_MM = 856.49 / 1000.0
SCALE_FACTORS = {"LM": 0.9686, "N": 0.9628}
BASE_DEFAULTS = dict(
    ngrid=1024, pupil_img_size=39.9988, diam_nominal=38.542,
    diam_ext=36.067, diam_int=10.8602932, seg_width=1.417,
    spi_angles=[0, 60, 120, 180, 240, 300], spi_width=0.202,
    AP_angles=[180], AP_width=0.386, AP_center=1.0, AP_length=2.0,
    dRext=0.0, dRint=0.0, dRspi=0.0,
    force_sym=False, circ_ext=False, circ_int=False,
)


class UniqueLoader(yaml.SafeLoader):
    """Reject duplicate keys rather than silently changing geometry."""


def _mapping(loader, node):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node)
        if key in result:
            raise ValueError(f"Duplicate YAML key: {key}")
        result[key] = loader.construct_object(value_node)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def _number(value, key, minimum=None, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{key} must be a finite number")
    if not np.isfinite(value) or (minimum is not None and value < minimum) or (positive and value <= 0):
        raise ValueError(f"Invalid {key}: {value}")


def load_config(path):
    path = Path(path).expanduser().resolve()
    raw = path.read_bytes()
    config = yaml.load(raw, Loader=UniqueLoader)
    allowed = {"name", "mask_family", "hole_radius_mm", "hole_centers_file", "npupil", "base_stop", "output_file"}
    if not isinstance(config, dict) or set(config) - allowed:
        raise ValueError("Configuration must be a mapping with supported keys only")
    for key in ("mask_family", "hole_radius_mm", "hole_centers_file", "output_file"):
        if key not in config:
            raise ValueError(f"Missing {key}")
    if config["mask_family"] not in SCALE_FACTORS:
        raise ValueError("mask_family must be LM or N")
    if "name" in config and (not isinstance(config["name"], str) or not config["name"].strip()):
        raise ValueError("name must be a nonempty string")
    _number(config["hole_radius_mm"], "hole_radius_mm", positive=True)
    coordinate_file = config["hole_centers_file"]
    if not isinstance(coordinate_file, str) or not coordinate_file.strip():
        raise ValueError("hole_centers_file must be a nonempty TXT path")
    coordinate_file = (path.parent / Path(os.path.expandvars(coordinate_file)).expanduser()).resolve()
    coordinate_raw = coordinate_file.read_bytes()
    centers = []
    for line_number, line in enumerate(coordinate_raw.decode("utf-8-sig").splitlines(), start=1):
        columns = line.split("#", 1)[0].split()
        if not columns:
            continue
        try:
            if len(columns) != 2:
                raise ValueError("expected exactly two columns: x_mm y_mm")
            pair = [float(value) for value in columns]
            for value in pair:
                _number(value, "coordinate")
        except ValueError as exc:
            raise ValueError(f"{coordinate_file}:{line_number}: {exc}") from exc
        centers.append(pair)
    if not centers:
        raise ValueError(f"No hole coordinates in {coordinate_file}")
    config["hole_centers_file"] = coordinate_file
    config["hole_centers_mm"] = centers
    config["coordinates_sha256"] = sha256(coordinate_raw).hexdigest()
    npupil = config.setdefault("npupil", 1385)
    if type(npupil) is not int or npupil < 3:
        raise ValueError("npupil must be an integer >= 3")
    base = config.get("base_stop", {})
    if not isinstance(base, dict) or set(base) - BASE_DEFAULTS.keys():
        raise ValueError("base_stop contains unsupported keys")
    base = {**BASE_DEFAULTS, **base}
    for key, value in base.items():
        if key in ("circ_ext", "circ_int", "force_sym"):
            if type(value) is not bool:
                raise ValueError(f"{key} must be boolean")
        elif key in ("spi_angles", "AP_angles"):
            if not isinstance(value, list):
                raise ValueError(f"{key} must be a list")
            for angle in value:
                _number(angle, key)
        elif key == "ngrid":
            if type(value) is not int or value < 2:
                raise ValueError("ngrid must be an integer >= 2")
        else:
            _number(value, key, minimum=0, positive=key in ("pupil_img_size", "diam_nominal", "seg_width"))
    if base["diam_int"] >= base["diam_ext"] or base["diam_ext"] > base["pupil_img_size"]:
        raise ValueError("Require diam_int < diam_ext <= pupil_img_size")
    if base["diam_int"] and base["seg_width"] / 2 > base["diam_int"]:
        raise ValueError("seg_width is too large for the central obscuration")
    config["base_stop"] = base
    if not isinstance(config["output_file"], str) or not config["output_file"].strip():
        raise ValueError("output_file must be a nonempty path")
    output = Path(os.path.expandvars(config["output_file"])).expanduser()
    config["output_file"] = (path.parent / output).resolve()
    if config["output_file"] in (path, coordinate_file):
        raise ValueError("output_file must not replace the configuration or coordinate TXT")
    config["config_sha256"] = sha256(raw).hexdigest()
    return config


def normalized_holes(config):
    """Keep the exact operation order of the existing aperture_mask.py."""
    centers = np.asarray(config["hole_centers_mm"]) * REF_SCALE * M_PER_MM
    diameter = (2 * config["hole_radius_mm"]) * REF_SCALE * M_PER_MM
    scale = SCALE_FACTORS[config["mask_family"]]
    nominal_radius = config["base_stop"]["diam_nominal"] / 2
    return scale * centers / nominal_radius, scale * (diameter / 2) / nominal_radius


def render_mask(config):
    # An adjacent checkout is optional; a copied script can use installed HEEPS.
    repo = Path(__file__).resolve().parents[1]
    if (repo / "heeps" / "__init__.py").is_file() and str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    import proper
    from heeps.pupil import create_stop

    base = dict(config["base_stop"])
    # The official routine shrinks a polygon's diameter by cos(arcsin(w / D)).
    # The notebook's local version disabled that correction. Compensate the input
    # to the official routine so both backends represent the requested diameter.
    active_lines = {line.split("#", 1)[0].strip() for line in inspect.getsource(create_stop).splitlines()}
    if not base["circ_ext"] and "diam_ext *= np.cos(alpha)" in active_lines:
        base["diam_ext"] = float(np.hypot(base["diam_ext"], base["seg_width"]))
        print("HEEPS polygon diameter correction compensated to preserve the configured outer geometry.")
    previous_n = proper.n
    try:
        pupil = create_stop(**base, npupil=config["npupil"],
                            add_seg=False, mask_name=None, verbose=False)
    finally:
        proper.n = previous_n
    centers, radius = normalized_holes(config)
    n = config["npupil"]
    grid = np.linspace(-1 + 1 / n, 1 - 1 / n, n)
    x, y = (np.float32(a) for a in np.meshgrid(grid, grid))
    holes = np.zeros_like(pupil)
    for index, (x0, y0) in enumerate(centers, start=1):
        hole = (x - x0)**2 + (y - y0)**2 < radius**2
        if not hole.any():
            raise ValueError(f"Hole {index} falls outside the grid or is smaller than a sampled pixel")
        if not np.any(pupil[hole] > 0):
            raise ValueError(f"Hole {index} is fully blocked by the base stop")
        holes[hole] = 1
    mask = np.asarray(holes * pupil, dtype=np.float32)
    if not np.isfinite(mask).all() or not np.any(mask > 0):
        raise ValueError("The generated mask has no valid transmitting pixels")
    return mask


def write_mask(config, mask, overwrite=False):
    """Atomically publish one FITS; existing files are protected by default."""
    output = config["output_file"]
    output.parent.mkdir(parents=True, exist_ok=True)
    header = fits.Header()
    header["MASKFAM"] = config["mask_family"]
    if "name" in config:
        header["MASKNAME"] = config["name"]
    header["NHOLES"] = len(config["hole_centers_mm"])
    header["RADMM"] = (config["hole_radius_mm"], "Common input radius in mask-plane mm")
    header["SCALEFAC"] = SCALE_FACTORS[config["mask_family"]]
    header["REFSCALE"] = REF_SCALE
    header["MPERMM"] = M_PER_MM
    header["DIAMNOM"] = (config["base_stop"]["diam_nominal"], "Nominal pupil diameter [m]")
    header["PUPLSIZE"] = (config["base_stop"]["pupil_img_size"], "Pupil image width [m]")
    header["CFGSHA"] = config["config_sha256"]
    header["COORDSHA"] = config["coordinates_sha256"]
    header["BUNIT"] = "amplitude transmission"
    header.add_history("HEEPS create_stop (no segments) multiplied by the union of circular holes.")
    header.add_history("Legacy normalized-coordinate mapping preserved; no flux renormalization.")
    for key, value in config["base_stop"].items():
        header.add_history(f"base_stop.{key} = {value}")
    centers = np.asarray(config["hole_centers_mm"])
    table = fits.BinTableHDU.from_columns([
        fits.Column(name="X_MM", format="D", unit="mm", array=centers[:, 0]),
        fits.Column(name="Y_MM", format="D", unit="mm", array=centers[:, 1]),
        fits.Column(name="RADIUS_MM", format="D", unit="mm",
                    array=np.full(len(centers), config["hole_radius_mm"])),
    ], name="HOLES")
    hdus = fits.HDUList([fits.PrimaryHDU(mask, header=header), table])
    descriptor, temporary = tempfile.mkstemp(prefix=".mask-", suffix=".fits", dir=output.parent)
    os.close(descriptor)
    try:
        hdus.writeto(temporary, overwrite=True, checksum=True)
        if overwrite:
            os.replace(temporary, output)
        else:
            os.link(temporary, output)  # Exclusive creation, including concurrent runs.
    finally:
        Path(temporary).unlink(missing_ok=True)
        hdus.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true", help="Validate YAML without importing HEEPS or writing files")
    parser.add_argument("--output-file", type=Path, help="Override output path, relative to the current directory")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace an existing output")
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        if args.output_file:
            config["output_file"] = args.output_file.expanduser().resolve()
        if config["output_file"] in (args.config.expanduser().resolve(), config["hole_centers_file"]):
            raise ValueError("Output must not replace the configuration or coordinate TXT")
        print(f"Family: {config['mask_family']}; scale: {SCALE_FACTORS[config['mask_family']]}; "
              f"holes: {len(config['hole_centers_mm'])}; radius: {config['hole_radius_mm']} mm; "
              f"shape: {config['npupil']} x {config['npupil']}")
        print(f"Output: {config['output_file']}")
        if args.dry_run:
            return 0
        if config["output_file"].exists() and not args.overwrite:
            raise FileExistsError("Output already exists; choose a new path or explicitly pass --overwrite")
        mask = render_mask(config)
        write_mask(config, mask, overwrite=args.overwrite)
        print(f"Saved aperture mask: {config['output_file']}")
        return 0
    except (ValueError, OSError, yaml.YAMLError, ImportError) as exc:
        parser.exit(1, f"Mask generation error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
