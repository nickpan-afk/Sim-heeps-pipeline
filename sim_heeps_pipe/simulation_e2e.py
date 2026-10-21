#!/usr/bin/env python3
"""Run OPD/amplitude -> HEEPS PSFs -> ScopeSim detector products from one YAML."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np
from astropy.io import fits
import yaml

if __package__:
    from . import heeps_e2e as heeps_stage
    from . import scopesim_e2e_Lband as scope_stage
else:
    import heeps_e2e as heeps_stage
    import scopesim_e2e_Lband as scope_stage


@dataclass(frozen=True)
class Pipeline:
    config_path: Path
    run_dir: Path
    heeps: heeps_stage.PipelineConfig
    scopesim: scope_stage.SimulationConfig
    summary: dict


def _path(value, base):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('output_dir must be a nonempty path string')
    path = Path(value).expanduser()
    return (base / path).resolve()


def load_config(path, *, run_mode=None, output_dir=None):
    """Build a read-only joint plan, including checks possible before PSFs exist."""
    path = Path(path).expanduser().resolve()
    data = yaml.load(path.read_text(), Loader=heeps_stage.UniqueKeyLoader)
    if not isinstance(data, dict):
        raise ValueError('The YAML document must be a mapping')
    allowed = {'name', 'output_dir', 'run_mode', 'heeps', 'scopesim'}
    if set(data) - allowed:
        raise ValueError(f'Unknown pipeline keys: {sorted(set(data) - allowed)}')
    if allowed - set(data):
        raise ValueError(f'Missing pipeline keys: {sorted(allowed - set(data))}')
    name = data['name']
    if not isinstance(name, str) or not name or name in {'.', '..'} or Path(name).name != name:
        raise ValueError('name must be a nonempty single directory name')
    mode = run_mode or data['run_mode']
    if mode not in {'preview', 'full'}:
        raise ValueError('run_mode must be preview or full')
    root = Path(output_dir).expanduser().resolve() if output_dir is not None else _path(data['output_dir'], path.parent)
    for section in ('heeps', 'scopesim'):
        if not isinstance(data[section], dict):
            raise ValueError(f'{section} must be a mapping')
    hdata, sdata = dict(data['heeps']), dict(data['scopesim'])
    reserved = {'heeps': {'name', 'output_dir', 'run_mode', 'preview_frames'},
                'scopesim': {'input_cube', 'output_dir', 'run_mode', 'pixel_scale_arcsec'}}
    for section, values in (('heeps', hdata), ('scopesim', sdata)):
        if set(values) & reserved[section]:
            raise ValueError(f'{section}: pipeline manages {sorted(set(values) & reserved[section])}')
    # Use existing validators without serializing intermediate configuration files.
    hdata.update(name=name, output_dir=str(root / '_planning'), run_mode='full')
    hconfig = heeps_stage.config_from_mapping(hdata, path.parent)
    sdata.update(input_cube=str(root / '_pending.fits'), output_dir=str(root / '_planning'),
                 run_mode='full', pixel_scale_arcsec=hconfig.pixel_scale_mas * 1e-3)
    sconfig = scope_stage.config_from_mapping(sdata, path.parent)
    if sconfig.effects.get('psf'):
        raise ValueError('scopesim.effects.psf must be false: HEEPS supplies the PSF')
    if mode == 'preview':
        sconfig = replace(sconfig, nplanes=1, mags=sconfig.mags[:1])
    minimum = sconfig.nplanes * sconfig.average_group * (2 if sconfig.make_calibrator else 1)
    # Empty sky still passes through the requested HEEPS stage, but uses only geometry.
    if sconfig.source_kind == 'empty_sky':
        minimum = 1
    if hconfig.nframes is not None and hconfig.nframes < minimum:
        raise ValueError(f'HEEPS nframes={hconfig.nframes} is insufficient; ScopeSim needs at least {minimum}')
    if mode == 'preview' or hconfig.nframes is None:
        hconfig = replace(hconfig, nframes=minimum)
    _, hsummary = heeps_stage.inspect_inputs(hconfig)
    frames = hsummary['effective_frames']
    groups, remainder = divmod(frames, sconfig.average_group)
    if sconfig.source_kind != 'empty_sky' and remainder and sconfig.remainder_policy == 'error':
        raise ValueError('HEEPS frame count is not divisible by scopesim.average_group; select remainder_policy: drop explicitly')
    missing = scope_stage.missing_instrument_packages(sconfig)
    if missing and not sconfig.download_missing_packages:
        raise FileNotFoundError(f'Missing instrument packages in {sconfig.package_dir}: {", ".join(missing)}')
    identity = dict(heeps=hsummary['fingerprint'], scopesim=scope_stage.config_dict(sconfig),
                    run_mode=mode,
                    code=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                          for p in (__file__, scope_stage.__file__)])
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    run_dir = root / f'{name}_{mode}_{digest}'
    hconfig = replace(hconfig, output_dir=run_dir / '.work' / 'heeps')
    plan, hsummary = heeps_stage.inspect_inputs(
        hconfig, work_dir=hconfig.output_dir, psf_file=run_dir / 'psf_cube.fits')
    sconfig = replace(sconfig, input_cube=plan.full_file, output_dir=run_dir)
    warnings = []
    sampled = sconfig.average_group * hconfig.cadence_s
    exposure = sconfig.dit * sconfig.ndit
    if not np.isclose(sampled, exposure, rtol=1e-9, atol=0):
        warnings.append('PSF group duration differs from detector exposure; DIT/NDIT remain explicit and are not inferred.')
    summary = dict(run_mode=mode, run_dir=str(run_dir), heeps=hsummary,
                   scopesim=dict(nplanes=sconfig.nplanes, average_group=sconfig.average_group,
                                 make_calibrator=sconfig.make_calibrator, source_kind=sconfig.source_kind,
                                 science_local_frames=None if sconfig.source_kind == 'empty_sky' else [0, sconfig.nplanes*sconfig.average_group],
                                 calibrator_local_frames=[groups//2*sconfig.average_group, (groups//2+sconfig.nplanes)*sconfig.average_group]
                                 if sconfig.make_calibrator and sconfig.source_kind != 'empty_sky' else None,
                                 pixel_scale_arcsec=sconfig.pixel_scale_arcsec,
                                 missing_instrument_packages=missing),
                   timing=dict(heeps_cadence_s=hconfig.cadence_s, psf_group_duration_s=sampled,
                               detector_exposure_s=exposure), warnings=warnings,
                   method='Broadband HEEPS PSF, normalized by ScopeSim; absolute photometry from source spectrum/magnitude.')
    return Pipeline(path, run_dir, hconfig, sconfig, summary)


def _write_json(path, data):
    temporary = path.with_name(f'.{path.name}.{uuid4().hex}.partial')
    try:
        temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _lock(run_dir):
    run_dir.mkdir(parents=True, exist_ok=True)
    lock = run_dir / '.pipeline.lock'
    try:
        handle = lock.open('x')
    except FileExistsError:
        raise RuntimeError(f'Run is locked: {lock}. If a previous process was killed, verify it has stopped before removing this lock.') from None
    try:
        with handle:
            handle.write(f'pid={os.getpid()}\n')
        yield
    finally:
        lock.unlink(missing_ok=True)


def _input_sampling(cube, pipeline):
    with fits.open(cube, memmap=True) as hdul:
        header = hdul[0].header
        scale = float(header['PSCALE']) * 1e-3
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('HEEPS output has invalid PSCALE')
        for axis in (1, 2):
            if header.get(f'CUNIT{axis}') != 'arcsec' or not np.isclose(float(header[f'CDELT{axis}']), scale, rtol=1e-9, atol=0):
                raise ValueError('HEEPS output has inconsistent spatial sampling metadata')
        if not np.isclose(scale, pipeline.heeps.pixel_scale_mas * 1e-3, rtol=1e-9, atol=0):
            raise ValueError('HEEPS output sampling differs from the execution plan')
        if hdul[0].shape[0] != pipeline.heeps.nframes:
            raise ValueError('HEEPS output frame count differs from the execution plan')
        if not np.isclose(float(header['CDELT3']), pipeline.heeps.cadence_s, rtol=1e-9, atol=0):
            raise ValueError('HEEPS output cadence differs from the execution plan')
        return scale


def run_pipeline(pipeline, *, progress=print):
    """Reuse HEEPS validated chunks; each attempt starts a fresh ScopeSim run."""
    with _lock(pipeline.run_dir):
        manifest_path = pipeline.run_dir / 'manifest.json'
        previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        attempt = dict(started_utc=datetime.now(timezone.utc).isoformat(), status='running', stage='heeps')
        manifest = dict(schema_version=1, status='running', plan=pipeline.summary,
                        attempts=previous.get('attempts', []) + [attempt], products={})
        effective = dict(name=pipeline.heeps.name, run_mode=pipeline.summary['run_mode'],
                         heeps=heeps_stage.config_dict(pipeline.heeps),
                         scopesim=scope_stage.config_dict(pipeline.scopesim))
        (pipeline.run_dir / 'config_input.yaml').write_text(pipeline.config_path.read_text())
        (pipeline.run_dir / 'config_effective.yaml').write_text(yaml.safe_dump(effective, sort_keys=False))
        _write_json(manifest_path, manifest)
        try:
            if progress:
                progress('Stage 1/2: HEEPS PSF propagation (valid existing chunks are reused).')
            cube, result = heeps_stage.run_pipeline(
                pipeline.heeps, pipeline.config_path, work_dir=pipeline.heeps.output_dir,
                psf_file=pipeline.run_dir / 'psf_cube.fits')
            scale = _input_sampling(cube, pipeline)
            scope_config = replace(pipeline.scopesim, input_cube=cube, pixel_scale_arcsec=scale)
            manifest['products'].update(psf_cube=str(cube), heeps_result=str(result))
            attempt['stage'] = 'scopesim'
            _write_json(manifest_path, manifest)
            if progress:
                progress('Stage 2/2: ScopeSim detector simulation.')
            attempt_number = len(manifest['attempts'])
            scope_manifest = scope_stage.run_simulation(
                scope_config, progress=progress, run_dir=pipeline.run_dir,
                metadata_dir=pipeline.run_dir / '.work' / f'scopesim_run{attempt_number:02d}',
                product_suffix='' if attempt_number == 1 else f'_run{attempt_number:02d}')
            report = json.loads(scope_manifest.read_text())
            if report['status'] != 'complete':
                raise RuntimeError('ScopeSim did not report a completed run')
            products = [{key: str((scope_manifest.parent / value).resolve()) if key in {'science', 'calibrator'} else value
                         for key, value in product.items()} for product in report['products']]
            manifest['products'].update(scopesim_manifest=str(scope_manifest), detector_products=products)
            attempt.update(status='complete', stage='complete', scopesim_manifest=str(scope_manifest))
            manifest['status'] = 'complete'
        except BaseException as exc:
            attempt.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed',
                           error=f'{type(exc).__name__}: {exc}')
            manifest['status'] = attempt['status']
            raise
        finally:
            attempt['finished_utc'] = datetime.now(timezone.utc).isoformat()
            _write_json(manifest_path, manifest)
        return manifest_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output-dir', type=Path, help='Override output root, relative to shell cwd')
    parser.add_argument('--dry-run', action='store_true', help='Check inputs and joint plan; no simulation, downloads, or outputs')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--preview', action='store_true')
    modes.add_argument('--full', action='store_true')
    args = parser.parse_args(argv)
    try:
        mode = 'preview' if args.preview else 'full' if args.full else None
        pipeline = load_config(args.config, run_mode=mode, output_dir=args.output_dir)
        print(json.dumps(pipeline.summary, indent=2))
        if args.dry_run:
            return 0
        path = run_pipeline(pipeline)
        print(f'PIPELINE_MANIFEST={path}')
        return 0
    except KeyboardInterrupt:
        print('Interrupted. Rerun the same command to reuse valid HEEPS outputs.', file=sys.stderr)
        return 130
    except (ValueError, TypeError, OSError, RuntimeError, KeyError, yaml.YAMLError) as exc:
        print(f'Pipeline error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
