#!/usr/bin/env python3
"""HEEPS PSF cube -> METIS detector FITS products, configured by YAML.

Standalone simulation entry point extracted from TryScopeL6clean.ipynb.
No imports from the notebook or project scripts, and no downstream analysis.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field, fields, replace
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import inspect
import json
import os
from numbers import Real
from pathlib import Path
import platform
import sys
import tempfile
from uuid import uuid4

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import scopesim as sim
import yaml


DEFAULT_EFFECTS = {
    "psf": False, "shot_noise": True, "readout_noise": True,
    "dark_current": True, "detector_linearity": True, "ipc": True,
    "cold_stop": True, "filter_wheel": True, "nd_filter_wheel": True,
    "skycalc_atmosphere": True, "telescope_reflection": True,
    "common_fore_optics": True, "img_lm_optics": True,
    "quantum_efficiency": True, "reference_pixel_mask": True,
    "ad_conversion": True, "auto_exposure": True,
    "exposure_integration": True, "exposure_output": True,
    "detector_readout_parameters": True,
}
REQUIRED_PACKAGES = ("METIS", "ELT", "Armazones")


@dataclass(frozen=True)
class SimulationConfig:
    input_cube: Path
    package_dir: Path
    output_dir: Path
    download_missing_packages: bool = True
    package_release: str = "stable"
    source_kind: str = "star"
    mags: tuple[float, ...] = (7.0,)
    nplanes: int = 1
    average_group: int = 1
    dit: float = 0.1
    ndit: int = 1
    mode: str = "img_lm"
    filter_name: str = "Lp"
    pixel_scale_arcsec: float = 0.00547
    cold_stop_transmission: float = 0.07511
    flux_scale: float = 1.0
    seed: int | None = None
    companion_mag: float | None = None
    companion_separation_arcsec: float = 0.07
    companion_pa_deg: float = -113.0
    nan_policy: str = "omit"
    remainder_policy: str = "error"
    run_mode: str = "full"
    make_calibrator: bool = True
    detector_readout_mode: str | None = None
    exposure_output_mode: str = "average"
    show_effects: bool = True
    effects: dict[str, bool] = field(default_factory=lambda: dict(DEFAULT_EFFECTS))

    def __post_init__(self):
        for name in ("input_cube", "package_dir", "output_dir"):
            object.__setattr__(self, name, Path(getattr(self, name)))
        values = self.mags
        if isinstance(values, Real) and not isinstance(values, bool):
            values = (values,)
        if not isinstance(values, (list, tuple)) or not values:
            raise ValueError("mags must be a nonempty list of finite numbers (or one number).")
        if any(isinstance(v, bool) or not isinstance(v, Real) or not np.isfinite(v) for v in values):
            raise ValueError("mags must contain finite numbers.")
        object.__setattr__(self, "mags", tuple(float(v) for v in values))
        if len({f"{mag:g}" for mag in self.mags}) != len(self.mags):
            raise ValueError("mags must produce distinct output filename tags.")
        for name in ("mode", "filter_name", "package_release"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"{name} must be a nonempty string.")
        if self.source_kind not in {"star", "empty_sky"}:
            raise ValueError("source_kind must be star or empty_sky.")
        if self.run_mode not in {"preview", "full"}:
            raise ValueError("run_mode must be preview or full.")
        if self.exposure_output_mode not in {"average", "sum"}:
            raise ValueError("exposure_output_mode must be average or sum.")
        if self.detector_readout_mode not in {None, "fast", "slow"}:
            raise ValueError("detector_readout_mode must be null, fast, or slow for LM imaging.")
        for name in ("nplanes", "average_group", "ndit"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")
        for name in ("dit", "pixel_scale_arcsec", "flux_scale"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive number.")
        for name in ("cold_stop_transmission", "companion_separation_arcsec", "companion_pa_deg"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite number.")
        if not 0 <= self.cold_stop_transmission <= 1 or self.companion_separation_arcsec < 0:
            raise ValueError("Transmission must be in [0,1], and companion separation nonnegative.")
        if self.seed is not None and (isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0):
            raise ValueError("seed must be a nonnegative integer or null.")
        if self.companion_mag is not None:
            if isinstance(self.companion_mag, bool) or not isinstance(self.companion_mag, Real) or not np.isfinite(self.companion_mag):
                raise ValueError("companion_mag must be finite or null.")
            if self.source_kind == "empty_sky":
                raise ValueError("empty_sky cannot include a companion.")
        if self.nan_policy not in {"omit", "error"} or self.remainder_policy not in {"drop", "error"}:
            raise ValueError("Invalid NaN or incomplete-group policy.")
        for name in ("make_calibrator", "show_effects", "download_missing_packages"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be true or false.")
        if not isinstance(self.effects, dict) or any(not isinstance(k, str) or not isinstance(v, bool) for k, v in self.effects.items()):
            raise ValueError("effects must be a mapping from effect names to true/false.")
        effects = dict(DEFAULT_EFFECTS)
        effects.update(self.effects)
        for name in ("exposure_integration", "exposure_output"):
            if not effects[name]:
                raise ValueError(f"Keep {name} enabled; use exposure_output_mode to select average/sum.")
        object.__setattr__(self, "effects", effects)


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject repeated YAML keys instead of silently using the last value."""


def _unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f"Duplicate YAML key: {key!r}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def load_config(path: Path) -> SimulationConfig:
    path = Path(path).expanduser().resolve()
    data = yaml.load(path.read_text(), Loader=UniqueKeyLoader)
    return config_from_mapping(data, path.parent)


def config_from_mapping(data, base) -> SimulationConfig:
    """Validate configuration data; resolve all paths against an explicit base."""
    base = Path(base).expanduser().resolve()
    if not isinstance(data, dict):
        raise ValueError("The YAML document must contain a configuration mapping.")
    data = dict(data)
    allowed = {f.name for f in fields(SimulationConfig)}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError("Unknown configuration keys: " + ", ".join(sorted(map(str, unknown))))
    for name in ("input_cube", "package_dir", "output_dir"):
        if name not in data or not isinstance(data[name], str) or not data[name].strip():
            raise ValueError(f"{name} must be a nonempty path string.")
        value = Path(data[name]).expanduser()
        data[name] = (value if value.is_absolute() else base / value).resolve()
    return SimulationConfig(**data)


def effective_config(config):
    if config.run_mode == "preview":
        return replace(config, run_mode="full", nplanes=1, mags=(config.mags[0],), output_dir=config.output_dir / "previews")
    return config


def config_dict(config):
    return json.loads(json.dumps(asdict(config), default=str))


def missing_instrument_packages(config):
    missing = []
    for name in REQUIRED_PACKAGES:
        path = config.package_dir / name
        if path.exists() and not path.is_dir():
            raise ValueError(f"Instrument package path is not a directory: {path}")
        if not path.is_dir():
            missing.append(name)
    return missing


def ensure_instrument_packages(config, *, progress=print):
    """Download missing packages to staging, then publish complete directories.

    Existing local packages are reused, never automatically upgraded or replaced.
    Download/unpack failures leave no partial package directories to reuse.
    """
    missing = missing_instrument_packages(config)
    if not missing:
        return
    if not config.download_missing_packages:
        raise FileNotFoundError(
            f"Missing instrument packages in {config.package_dir}: {', '.join(missing)}. "
            "Enable download_missing_packages or install them manually.")
    config.package_dir.mkdir(parents=True, exist_ok=True)
    if progress:
        progress(f"Downloading {', '.join(missing)} ({config.package_release}) to {config.package_dir}")
    try:
        with tempfile.TemporaryDirectory(prefix=".download_", dir=config.package_dir) as temporary:
            stage = Path(temporary)
            sim.download_packages(missing, release=config.package_release, save_dir=stage)
            for name in missing:
                path = stage / name
                if not path.is_dir() or not any(path.glob("*.yaml")):
                    raise RuntimeError(f"Download did not produce a complete {name} package.")
            for name in missing:
                target = config.package_dir / name
                if target.exists():
                    raise RuntimeError(f"Package appeared during download; refusing to replace {target}")
                (stage / name).rename(target)
    except Exception as exc:
        raise RuntimeError(
            f"Instrument package setup failed: {exc}. Existing local packages were retained. "
            "Check the network connection and package_release, then rerun.") from exc


def inspect_input(config, *, require_packages=True):
    if not config.input_cube.is_file():
        raise FileNotFoundError(f"Input cube not found: {config.input_cube}")
    missing = missing_instrument_packages(config)
    if require_packages and missing:
        raise FileNotFoundError(f"Missing instrument packages in {config.package_dir}: {', '.join(missing)}")
    with fits.open(config.input_cube, memmap=True) as hdul:
        shape = hdul[0].shape
    if len(shape) != 3 or min(shape) < 1:
        raise ValueError(f"HDU 0 must be a nonempty (frame,y,x) PSF cube, got {shape}.")
    groups, remainder = divmod(shape[0], config.average_group)
    sky = config.source_kind == "empty_sky"
    cal_start = groups // 2
    available = cal_start if config.make_calibrator else groups
    if not sky and remainder and config.remainder_policy == "error":
        raise ValueError(f"{remainder} input frames form an incomplete group; explicitly select remainder_policy: drop to discard them.")
    if not sky and config.nplanes > available:
        raise ValueError(f"Requested {config.nplanes} planes, but only {available} {'paired' if config.make_calibrator else 'science'} groups are available.")
    return {
        "missing_instrument_packages": missing,
        "download_missing_packages": config.download_missing_packages,
        "package_release": config.package_release,
        "shape": list(shape), "source_kind": config.source_kind,
        "input_usage": "shape/crop geometry only" if sky else "normalized PSF frame groups",
        "cal_start_group": cal_start,
        "science_frames": None if sky else [0, config.nplanes * config.average_group],
        "calibrator_frames": None if sky or not config.make_calibrator else [cal_start * config.average_group, (cal_start + config.nplanes) * config.average_group],
        "range_convention": "zero-based, stop excluded",
        "trailing_frames_dropped": remainder if not sky else 0,
        "output_shape": [config.nplanes, shape[1], shape[2]],
        "seconds_per_output_plane": config.dit * config.ndit,
        "seconds_per_target": config.nplanes * config.dit * config.ndit,
        "note": "DIT/NDIT are explicit; input cadence and wavelength sampling are not inferred.",
    }


def normalized_group(data, group, config):
    start, stop = group * config.average_group, (group + 1) * config.average_group
    if group < 0 or stop > data.shape[0]:
        raise IndexError("Input group lies outside the PSF cube.")
    total = np.zeros(data.shape[1:], dtype=np.float64)
    count = np.zeros(data.shape[1:], dtype=np.int64)
    for frame in data[start:stop]:
        finite = np.isfinite(frame)
        if np.isinf(frame).any() or (config.nan_policy == "error" and not finite.all()):
            raise ValueError(f"Invalid PSF samples in group {group}.")
        if np.any(frame[finite] < 0):
            raise ValueError(f"Negative PSF weights in group {group}.")
        total += np.where(finite, frame, 0)
        count += finite
    if np.any(count == 0):
        raise ValueError(f"No finite samples for one or more pixels in group {group}.")
    frame = total / count
    flux = frame.sum(dtype=np.float64)
    if not np.isfinite(flux) or flux <= 0:
        raise ValueError(f"Nonpositive or non-finite PSF total in group {group}.")
    return frame / flux * config.flux_scale


def image_wcs(header, nx, ny, scale):
    # Match the notebook's explicit centered LINEAR/first-two-axis WCS.
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [(nx + 1) / 2, (ny + 1) / 2]
    wcs.wcs.crval = [0, 0]
    wcs.wcs.cdelt = [scale, scale]
    wcs.wcs.cunit = ["arcsec", "arcsec"]
    wcs.wcs.ctype = [header.get("CTYPE1", "LINEAR"), header.get("CTYPE2", "LINEAR")]
    return wcs


def frame_source(frame, wcs, mag):
    hdu = fits.ImageHDU(frame, header=wcs.to_header())
    spectrum = sim.source.source_templates.vega_spectrum(mag=mag)
    return sim.Source(image_hdu=hdu, spectra=[spectrum])


def input_source(data, group, wcs, mag, config, science=True):
    if config.source_kind == "empty_sky":
        return sim.source.source_templates.empty_sky(flux=0)
    frame = normalized_group(data, group, config)
    source = frame_source(frame, wcs, mag)
    if science and config.companion_mag is not None:
        companion = frame_source(frame.copy(), wcs, config.companion_mag)
        angle = np.radians(config.companion_pa_deg)
        companion.shift(dx=config.companion_separation_arcsec * np.sin(angle),
                        dy=config.companion_separation_arcsec * np.cos(angle))
        source = source + companion
    return source


def build_instrument(config):
    sim.link_irdb(str(config.package_dir))
    properties = {"!OBS.filter_name": config.filter_name}
    if config.detector_readout_mode is not None:
        properties["!OBS.detector_readout_mode"] = config.detector_readout_mode
    cmd = sim.UserCommands(use_instrument="METIS", set_modes=[config.mode], properties=properties)
    if config.seed is not None:
        cmd["!SIM.random.seed"] = config.seed
    metis = sim.OpticalTrain(cmd)
    metis["psf"].include = False
    metis["cold_stop"].include = True
    cold_stop = metis["cold_stop"]
    cold_stop.update_transmission(config.cold_stop_transmission)
    # ScopeSim 0.11.2 update_transmission makes unitless emission=[0,0].
    # Units on an identically-zero emission repair Synphot without changing flux.
    emission = cold_stop.surface.meta.get("emission")
    if emission is not None and not cold_stop.surface.meta.get("emission_unit") and np.all(np.asarray(emission) == 0):
        unit = "ph s-1 cm-2 Angstrom-1 arcsec-2"
        cold_stop.surface.meta["emission_unit"] = unit
        cold_stop.meta["emission_unit"] = unit
    metis["detector_readout_parameters"].include = True
    if "chop_nod" in metis.effects["name"]:
        metis["chop_nod"].include = False
    metis.update()
    for name, enabled in config.effects.items():
        try:
            metis[name].include = enabled
        except (KeyError, ValueError) as exc:
            raise ValueError(f"Cannot configure effect {name!r}; check the instrument package and spelling.") from exc
    metis["exposure_output"].set_mode(config.exposure_output_mode)
    metis.update()
    return metis


def detector_frame(metis, source, config, readout_id, nx, ny):
    # Preserve the seed scheme used by the notebook helper, including 0.11 caches.
    if config.seed is not None:
        for index, effect in enumerate(metis.optics_manager.all_effects):
            if "random_seed" in effect.meta:
                effect.meta["random_seed"] = int(np.random.SeedSequence([config.seed, readout_id, index]).generate_state(1)[0])
    metis.observe(source, update=True)
    kwargs = {"dit": config.dit, "ndit": config.ndit}
    if "roid" in inspect.signature(metis.readout).parameters:
        kwargs["roid"] = readout_id
    outputs = metis.readout(**kwargs)
    try:
        if len(outputs) != 1 or len(outputs[0]) < 2:
            raise ValueError("Expected one detector array with image data in HDU 1.")
        hdul = outputs[0]
        data = hdul[1].data
        if data.ndim != 2 or ny > data.shape[0] or nx > data.shape[1]:
            raise ValueError(f"Cannot center-crop detector shape {data.shape} to {(ny,nx)}.")
        x0, y0 = (data.shape[1] - nx) // 2, (data.shape[0] - ny) // 2
        image = np.array(data[y0:y0+ny, x0:x0+nx], dtype=np.float32, copy=True)
        if not np.isfinite(image).all():
            raise ValueError("Detector output contains non-finite pixels.")
        header = hdul[0].header.copy()
        header.extend(hdul[1].header, update=True)
        for key in list(header):
            if key.startswith("CRPIX1"):
                header[key] -= x0
            elif key.startswith("CRPIX2"):
                header[key] -= y0
        header["PIXSCALE"] = (float(metis.cmds["!INST.pixel_scale"]), "Detector arcsec/pixel")
        header["DETOMODE"] = metis["exposure_output"].current_mode
        if not metis["ad_conversion"].include:
            header["BUNIT"] = "electron"
        return image, header
    finally:
        for hdul in outputs:
            hdul.close()


def output_header(header, config, plan, mag, kind, readout_id):
    header = header.copy()
    header["SRCTYPE"] = config.source_kind
    if mag is not None:
        header["SRCMAG"] = float(mag)
    else:
        header.remove("SRCMAG", ignore_missing=True)
    for key, value in {"AVGGRP": config.average_group, "NPLANES": config.nplanes,
                       "DIT": config.dit, "NDIT": config.ndit,
                       "INTTIME": config.dit * config.ndit,
                       "TOTINT": config.nplanes * config.dit * config.ndit,
                       "FILTER": config.filter_name, "SIMMODE": config.mode,
                       "COLDSTOP": config.cold_stop_transmission,
                       "FLXSCALE": config.flux_scale, "INPSCALE": config.pixel_scale_arcsec,
                       "OBJKIND": kind, "RNGSEED": str(config.seed),
                       "ROID0": readout_id, "ROIDSTEP": 2,
                       "SCOPSIM": importlib.metadata.version("scopesim"),
                       "MANIFEST": "manifest.json"}.items():
        header[key] = value
    companion = kind == "science" and config.companion_mag is not None
    header["COMPAN"] = companion
    if companion:
        header["COMPMAG"] = config.companion_mag
        header["COMPSEP"] = config.companion_separation_arcsec
        header["COMPPA"] = config.companion_pa_deg
    interval = plan["science_frames" if kind == "science" else "calibrator_frames"]
    if interval is not None:
        header["INSTART"], header["INSTOP"] = interval
    if config.source_kind == "empty_sky":
        header.add_history("Zero-flux empty sky; input PSF pixels unused (crop geometry only).")
    else:
        header.add_history("HEEPS PSF groups averaged, normalized, then passed through ScopeSim.")
    header.add_history("No downstream analysis performed. See config_effective.yaml and manifest.json.")
    return header


def json_write(path, value):
    temporary = path.with_suffix(".json.partial")
    temporary.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + "\n")
    temporary.replace(path)


def provenance(config):
    versions = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "astropy", "synphot", "scopesim", "PyYAML"):
        versions[name] = importlib.metadata.version(name)
    packages = {}
    for name in ("METIS", "ELT", "Armazones"):
        root = config.package_dir / name
        packages[name] = {
            "version_file": (root / "version.yaml").read_text() if (root / "version.yaml").is_file() else None,
            "yaml_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(root.rglob("*")) if p.suffix in {".yaml", ".yml"}},
        }
    stat = config.input_cube.stat()
    return {"versions": versions, "instrument_packages": packages,
            "input": {"path": str(config.input_cube.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns},
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "seed_scheme": "SeedSequence(root_seed, readout_id, effect_index); roid if supported"}


def run_simulation(config, *, config_path=None, instrument_factory=None, progress=print,
                   run_dir=None, metadata_dir=None, product_suffix=''):
    """Simulate exposures; optional explicit directories support the joint driver."""
    if not isinstance(product_suffix, str) or any(
            c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in product_suffix):
        raise ValueError('product_suffix may contain only letters, digits, underscores and hyphens')
    config = effective_config(config)
    inspect_input(config, require_packages=False)  # Validate input before any download.
    ensure_instrument_packages(config, progress=progress)
    plan = inspect_input(config)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid4().hex[:8]
    explicit_dir = run_dir is not None
    run_dir = Path(run_dir).resolve() if explicit_dir else config.output_dir / run_id
    metadata_dir = Path(metadata_dir).resolve() if metadata_dir is not None else run_dir
    run_dir.mkdir(parents=True, exist_ok=explicit_dir)
    if metadata_dir != run_dir:
        metadata_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = metadata_dir / "manifest.json"
    if manifest_path.exists():
        raise FileExistsError(f'Refusing to overwrite existing stage manifest: {manifest_path}')
    manifest = {"schema_version": 1, "status": "running", "run_id": run_id,
                "config": config_dict(config), "plan": plan, "products": []}
    json_write(manifest_path, manifest)
    try:
        manifest["provenance"] = provenance(config)
        if config_path is not None:
            source = Path(config_path)
            (metadata_dir / "config_input.yaml").write_text(source.read_text())
            manifest["configuration_source"] = {"path": str(source.resolve()),
                "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
        (metadata_dir / "config_effective.yaml").write_text(yaml.safe_dump(config_dict(config), sort_keys=False))
        metis = (instrument_factory or build_instrument)(config)
        if config.show_effects:
            print("Detector mode:", metis.cmds["!DET.mode"])
            metis.effects.pprint_all()
        manifest["effective_effects"] = [{"name": str(row["name"]), "class": str(row["class"]),
                                           "include": bool(row["included"])} for row in metis.effects]
        json_write(manifest_path, manifest)
        with fits.open(config.input_cube, memmap=True) as hdul:
            data = hdul[0].data
            _, ny, nx = data.shape
            wcs = image_wcs(hdul[0].header, nx, ny, config.pixel_scale_arcsec)
            mags = (None,) if config.source_kind == "empty_sky" else config.mags
            kinds = ("science", "calibrator") if config.make_calibrator else ("science",)
            for mag_index, mag in enumerate(mags):
                tag = "empty_sky" if mag is None else f"mag{mag:g}"
                with tempfile.TemporaryDirectory(prefix=".work_", dir=run_dir) as scratch:
                    arrays = {kind: np.lib.format.open_memmap(Path(scratch) / f"{kind}.npy", mode="w+",
                               dtype=np.float32, shape=(config.nplanes, ny, nx)) for kind in kinds}
                    headers = {}
                    try:
                        for plane in range(config.nplanes):
                            for kind_index, kind in enumerate(kinds):
                                group = plane if kind == "science" else plan["cal_start_group"] + plane
                                source = input_source(data, group, wcs, mag, config, science=kind == "science")
                                rid = 2 * (mag_index * config.nplanes + plane) + kind_index
                                image, header = detector_frame(metis, source, config, rid, nx, ny)
                                arrays[kind][plane] = image
                                if kind not in headers:
                                    headers[kind] = header
                            if progress:
                                progress(f"{tag}: plane {plane + 1}/{config.nplanes}")
                        product = {"mag": mag, "source_kind": config.source_kind}
                        for kind_index, kind in enumerate(kinds):
                            path = run_dir / f"{kind}_{tag}{product_suffix}.fits"
                            if path.exists():
                                raise FileExistsError(f'Refusing to overwrite existing detector product: {path}')
                            header = output_header(headers[kind], config, plan, mag, kind,
                                                   2 * mag_index * config.nplanes + kind_index)
                            partial = path.with_suffix(".fits.partial")
                            fits.writeto(partial, arrays[kind], header=header, overwrite=False)
                            partial.replace(path)
                            product[kind] = os.path.relpath(path, metadata_dir)
                        manifest["products"].append(product)
                        json_write(manifest_path, manifest)
                    finally:
                        for array in arrays.values():
                            array.flush()
                            array._mmap.close()
        manifest["status"] = "complete"
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        json_write(manifest_path, manifest)
    if progress:
        progress(f"Completed. Manifest: {manifest_path}")
        for product in manifest["products"]:
            for kind in ("science", "calibrator"):
                if kind in product:
                    progress(f"  {kind}: {(metadata_dir / product[kind]).resolve()}")
    return manifest_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="YAML configuration; paths inside it are relative to the YAML file")
    parser.add_argument("--input", type=Path, help="Override input cube (relative to the current working directory)")
    parser.add_argument("--output-dir", type=Path, help="Override output directory (relative to the current working directory)")
    parser.add_argument("--nplanes", type=int, help="Override configured plane count (preview still limits it to 1)")
    parser.add_argument("--seed", type=int, help="Override root random seed")
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument("--dry-run", action="store_true", help="Validate configuration and input geometry; do not build an instrument or simulate")
    actions.add_argument("--inspect-effects", action="store_true", help="Build and print the configured instrument effects, then exit without simulating")
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        changes = {}
        for option, name in ((args.input, "input_cube"), (args.output_dir, "output_dir")):
            if option is not None:
                changes[name] = option.expanduser().resolve()
        for name in ("nplanes", "seed"):
            if getattr(args, name) is not None:
                changes[name] = getattr(args, name)
        config = effective_config(replace(config, **changes))
        plan = inspect_input(config, require_packages=not config.download_missing_packages)
        print("Effective configuration:")
        print(yaml.safe_dump(config_dict(config), sort_keys=False))
        print("Execution plan:")
        print(json.dumps(plan, indent=2))
        if args.dry_run:
            return 0
        if args.inspect_effects:
            ensure_instrument_packages(config)
            metis = build_instrument(config)
            print("Detector mode:", metis.cmds["!DET.mode"])
            metis.effects.pprint_all()
            return 0
        run_simulation(config, config_path=args.config)
        return 0
    except (ValueError, TypeError, FileNotFoundError, yaml.YAMLError) as exc:
        print(f"Configuration/input error: {exc}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(f"Simulation/setup error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
