#!/usr/bin/env python3
"""Local phase/amplitude FITS inputs -> non-coronagraphic HEEPS PSF cube.

Self-contained application: this script + YAML + importable HEEPS and dependencies.
YAML paths are relative to the YAML file. HEEPS is imported inside heeps_backend
only when propagation starts; --dry-run remains read-only and does not import it.
No notebook execution, implicit downloads, detector simulation, or analysis.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields, replace
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata
import importlib.util
import inspect
import json
import logging
import math
from numbers import Real
import os
from pathlib import Path
import re
import shutil
import sys
from typing import Callable
from types import FunctionType
from uuid import uuid4

# Use the adjacent HEEPS checkout when present; copied scripts use installed HEEPS
# or PYTHONPATH. No other project modules are imported.
_repo_candidate = Path(__file__).resolve().parents[1]
if (_repo_candidate / 'heeps' / '__init__.py').is_file() and str(_repo_candidate) not in sys.path:
    sys.path.insert(0, str(_repo_candidate))

import numpy as np
from astropy.io import fits
import yaml


# -----------------------------------------------------------------------------
# Configuration, planning, HEEPS propagation, and FITS output
# -----------------------------------------------------------------------------
LOG = logging.getLogger(__name__)
SCHEMA = 1

# Centers and sampling follow heeps.config.read_config. Preserve the existing L
# broadband prescription; the other profiles start monochromatic.
BAND_SPECS = {
    'L':   {'lambda_center_m': 3.81e-6,  'pixel_scale_mas': 5.47,  'frac_bw': 0.182, 'nlam': 5},
    'M':   {'lambda_center_m': 4.79e-6,  'pixel_scale_mas': 5.47,  'frac_bw': 0., 'nlam': 1},
    'N1':  {'lambda_center_m': 8.70e-6,  'pixel_scale_mas': 6.79,  'frac_bw': 0., 'nlam': 1},
    'N2':  {'lambda_center_m': 11.33e-6, 'pixel_scale_mas': 6.79,  'frac_bw': 0., 'nlam': 1},
}


@dataclass(frozen=True)
class Spectrum:
    wavelengths_m: tuple[float, ...]
    weights: tuple[float, ...]

    def __post_init__(self):
        wavelengths = np.asarray(self.wavelengths_m, dtype=float)
        weights = np.asarray(self.weights, dtype=float)
        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError('At least one wavelength is required')
        if weights.shape != wavelengths.shape:
            raise ValueError('Wavelengths and weights must have the same length')
        if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
            raise ValueError('Wavelengths must be finite and positive, in meters')
        if np.any(np.diff(wavelengths) <= 0):
            raise ValueError('Wavelengths must be strictly increasing')
        if not np.all(np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
            raise ValueError('Weights must be finite, non-negative, with positive sum')
        object.__setattr__(self, 'wavelengths_m', tuple(map(float, wavelengths)))
        object.__setattr__(self, 'weights', tuple(map(float, weights / weights.sum())))

    @classmethod
    def uniform(cls, lambda_center_m: float, frac_bw: float, nlam: int = 5):
        if isinstance(nlam, bool) or not isinstance(nlam, int) or nlam < 1:
            raise ValueError('nlam must be a positive integer')
        if not np.isfinite(frac_bw) or not 0 <= frac_bw < 2:
            raise ValueError('frac_bw must be in [0, 2)')
        if not np.isfinite(lambda_center_m) or lambda_center_m <= 0:
            raise ValueError('lambda_center_m must be finite and positive')
        if nlam == 1 or frac_bw == 0:
            return cls((lambda_center_m,), (1.,))
        return cls(tuple(np.linspace(lambda_center_m * (1 - frac_bw / 2),
                                     lambda_center_m * (1 + frac_bw / 2), nlam)),
                   (1.,) * nlam)


@dataclass(frozen=True)
class Source:
    """A source file whose first frame corresponds to global_start.

    Static amplitude maps can be 2D or a one-frame 3D cube.
    """
    filename: str
    global_start: int = 0

    def __post_init__(self):
        if not isinstance(self.filename, str) or not self.filename:
            raise ValueError('Source filename must be a non-empty string')
        if (isinstance(self.global_start, bool) or not isinstance(self.global_start, int)
                or self.global_start < 0):
            raise ValueError('Source global_start must be a non-negative integer')


@dataclass(frozen=True)
class Experiment:
    name: str
    band: str
    mask_name: str
    spectrum: Spectrum
    phase: tuple[Source, ...] = ()
    amplitude: tuple[Source, ...] = ()
    pupil: str = 'pupil/ELT_fullM1_feb2025_noseg.fits'
    start: int = 0
    stop: int = 36000
    chunk_size: int = 1000
    hfov_arcsec: float = 1.1
    ngrid: int = 1024
    cpu_count: int = 1
    dit_s: float = 0.3
    pscale_mas: float | None = None
    pure: bool = False
    lyot_stop: str | None = None
    phase_unit: str = 'm'

    def __post_init__(self):
        object.__setattr__(self, 'phase', tuple(self.phase))
        object.__setattr__(self, 'amplitude', tuple(self.amplitude))
        if not isinstance(self.spectrum, Spectrum) or not isinstance(self.pure, bool):
            raise ValueError('Require a Spectrum and an explicit boolean pure flag')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', self.name):
            raise ValueError('Experiment name must be a safe filename')
        valid_masks = {'L': ('LM6c', 'LM12c'), 'M': ('LM6c', 'LM12c'),
                       **{band: ('N6c', 'N12c') for band in ('N1', 'N2')}}
        if self.mask_name not in valid_masks.get(self.band, ()):
            raise ValueError('Incompatible band and mask_name')
        for key in ('start', 'stop', 'chunk_size', 'ngrid', 'cpu_count'):
            value = getattr(self, key)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f'{key} must be an integer')
        if self.start < 0 or self.stop <= self.start:
            raise ValueError('Require 0 <= start < stop')
        if self.chunk_size < 1 or self.cpu_count < 1 or self.ngrid < 3:
            raise ValueError('chunk_size, cpu_count and ngrid must be positive')
        for key in ('hfov_arcsec', 'dit_s', 'pscale_mas'):
            value = getattr(self, key)
            if value is not None and (not np.isfinite(value) or value <= 0):
                raise ValueError(f'{key} must be finite and positive')
        if self.pure and (self.phase or self.amplitude):
            raise ValueError('Pure experiments must not specify error inputs')
        if not self.pure and not (self.phase or self.amplitude):
            raise ValueError('Non-pure experiments require phase and/or amplitude inputs')
        if self.phase_unit not in ('m', 'um', 'nm'):
            raise ValueError('phase_unit must be m, um, or nm of optical path difference')
        for sources in (self.phase, self.amplitude):
            starts = [s.global_start for s in sources]
            if starts != sorted(set(starts)) or any(s < 0 for s in starts):
                raise ValueError('Source global_start values must be unique, sorted and non-negative')

    def preview(self, nframes: int = 2):
        """An isolated small run; the production experiment remains unchanged."""
        if not isinstance(nframes, int) or nframes < 1:
            raise ValueError('nframes must be positive')
        return replace(self, name=self.name + '_preview', stop=min(self.start + nframes, self.stop),
                       chunk_size=nframes, cpu_count=1)


def _resolve(root: Path, filename: str) -> Path:
    path = Path(os.path.expandvars(filename)).expanduser()
    return (root / path).resolve()


def _shape(path: Path) -> tuple[int, ...]:
    with fits.open(path, memmap=True) as hdul:
        h = hdul[0].header
        return tuple(h[f'NAXIS{i}'] for i in range(h['NAXIS'], 0, -1))


def _identity(path: Path):
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def _digest(value) -> str:
    return sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _runtime_versions():
    versions = {'python': sys.version, 'numpy': np.__version__}
    for name in ('astropy', 'proper', 'scipy', 'scikit-image'):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = 'unavailable'
    return versions


def _code_digest():
    """Hash this script and the importable HEEPS sources, without importing HEEPS.

    This remains valid when the script is copied outside the repository. No other
    application helper directory participates in execution or cache provenance.
    """
    spec = importlib.util.find_spec('heeps')
    if spec is None or not spec.submodule_search_locations:
        raise ModuleNotFoundError('HEEPS is not importable; install HEEPS or add its parent directory to PYTHONPATH')
    digest = sha256(json.dumps(_runtime_versions(), sort_keys=True).encode())
    digest.update(b'heeps_e2e.py')
    digest.update(Path(__file__).read_bytes())
    for location in spec.submodule_search_locations:
        package_root = Path(location)
        for path in sorted(package_root.rglob('*.py')):
            digest.update(('heeps/' + path.relative_to(package_root).as_posix()).encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def _source_at(sources, index):
    candidates = [s for s in sources if s.global_start <= index]
    if not candidates:
        raise ValueError(f'No source covers global frame {index}')
    return candidates[-1]


@dataclass(frozen=True)
class Plan:
    experiment: Experiment
    input_root: Path
    run_dir: Path
    fingerprint: str
    chunks: tuple[tuple[int, int], ...]
    inputs: tuple[dict, ...]
    code_hash: str
    output_file: Path | None = None

    @property
    def full_file(self):
        return self.output_file if self.output_file is not None else self.run_dir / 'full.fits'

    def summary(self):
        e = self.experiment
        return dict(experiment=e.name, band=e.band, mask=e.mask_name, pure=e.pure,
                    frame_range=[e.start, e.stop], chunk_count=len(self.chunks),
                    wavelengths_m=e.spectrum.wavelengths_m, weights=e.spectrum.weights,
                    hfov_arcsec=e.hfov_arcsec, cpu_count=e.cpu_count,
                    fingerprint=self.fingerprint, output=str(self.run_dir))


def plan_experiment(experiment: Experiment, input_root, output_root) -> Plan:
    """Read headers/stat only; never creates output or imports HEEPS."""
    root = Path(input_root).expanduser().resolve()
    e = experiment
    files = {_resolve(root, e.pupil)}
    if e.lyot_stop is not None:
        mask_path = _resolve(root, e.lyot_stop)
        if not mask_path.is_file():
            raise FileNotFoundError(mask_path)
        if len(_shape(mask_path)) != 2:
            raise ValueError('Explicit Lyot stop must be a 2D FITS image')
        files.add(mask_path)
    # Include optical calibration files, including the auto-selected Lyot stops.
    files.update(p for p in (root / 'pupil').rglob('*') if p.is_file())
    for sources, kind in ((e.phase, 'phase'), (e.amplitude, 'amplitude')):
        for i, source in enumerate(sources):
            path = _resolve(root, source.filename)
            if not path.is_file():
                raise FileNotFoundError(path)
            shape = _shape(path)
            if not shape or min(shape) < 1:
                raise ValueError(f'{kind} must be non-empty: {path}: {shape}')
            static = kind == 'amplitude' and (len(shape) == 2 or (len(shape) == 3 and shape[0] == 1))
            if not static and len(shape) != 3:
                raise ValueError(f'{kind} must be a 3D cube: {path}: {shape}')
            begin = max(e.start, source.global_start)
            end = min(e.stop, sources[i + 1].global_start if i + 1 < len(sources) else e.stop)
            if end > begin and not static and end - source.global_start > shape[0]:
                raise ValueError(f'{path} does not cover global frames {begin}:{end}')
            files.add(path)
        if sources and sources[0].global_start > e.start:
            raise ValueError(f'{kind} inputs do not cover frame {e.start}')
    if not _resolve(root, e.pupil).is_file():
        raise FileNotFoundError(_resolve(root, e.pupil))
    if len(_shape(_resolve(root, e.pupil))) != 2:
        raise ValueError('Pupil must be a 2D FITS image')
    inputs = tuple(_identity(p) for p in sorted(files))
    code_hash = _code_digest()
    fingerprint = _digest(dict(schema=SCHEMA, experiment=asdict(e), inputs=inputs, code_hash=code_hash))
    boundaries = {e.start, e.stop}
    boundaries.update(range(e.start, e.stop, e.chunk_size))
    boundaries.update(s.global_start for s in (*e.phase, *e.amplitude) if e.start < s.global_start < e.stop)
    edges = sorted(boundaries)
    chunks = tuple(zip(edges[:-1], edges[1:]))
    run_dir = Path(output_root).expanduser().resolve() / e.name / fingerprint[:16]
    if run_dir == root or root in run_dir.parents:
        raise ValueError('Output must be outside input_root')
    return Plan(e, root, run_dir, fingerprint, chunks, inputs, code_hash)


def _assert_current(plan):
    for item in plan.inputs:
        if _identity(Path(item['path'])) != item:
            raise RuntimeError(f'Input changed after planning: {item["path"]}; create a new plan')
    if _code_digest() != plan.code_hash:
        raise RuntimeError('Workflow or HEEPS code changed after planning; create a new plan')


@contextmanager
def _run_lock(plan):
    plan.run_dir.mkdir(parents=True, exist_ok=True)
    lock = plan.run_dir / '.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise RuntimeError(f'Run is locked: {lock}. Check the recorded PID before removing a stale lock.') from None
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(str(os.getpid()))
        yield
    finally:
        lock.unlink()


def _atomic_json(path, payload):
    tmp = path.with_suffix(path.suffix + '.partial')
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + '\n')
    os.replace(tmp, path)


def _header(plan, start, stop, pscale, actual_lams):
    e = plan.experiment
    h = fits.Header()
    for key, value in dict(RUNHASH=plan.fingerprint, COMPLETE=True, FRAME0=start, FRAME1=stop,
                           BAND=e.band, MASK=e.mask_name, PURE=e.pure, DIT=e.dit_s,
                           PSCALE=pscale, NLAM=len(e.spectrum.weights), NORM='WEIGHTED_MEAN',
                           BUNIT='relative intensity', REPEAT=e.pure,
                           ADDPHASE=bool(e.phase), ADDAMP=bool(e.amplitude), WFEUNIT='m').items():
        h[key] = value
    for i, (requested, actual, weight) in enumerate(zip(e.spectrum.wavelengths_m, actual_lams, e.spectrum.weights)):
        h[f'HIERARCH WAVE{i} REQ'] = requested
        h[f'HIERARCH WAVE{i} ACT'] = actual
        h[f'HIERARCH WAVE{i} WEIGHT'] = weight
    h['CTYPE1'] = 'LINEAR'
    h['CTYPE2'] = 'LINEAR'
    h['CTYPE3'] = 'TIME'
    h['CUNIT1'] = h['CUNIT2'] = 'arcsec'
    h['CUNIT3'] = 's'
    h['CDELT1'] = h['CDELT2'] = pscale * 1e-3
    h['CDELT3'] = e.dit_s
    h['CRVAL1'] = h['CRVAL2'] = 0.
    h['CRVAL3'] = start * e.dit_s
    h['CRPIX3'] = 1.
    h.add_history('Temporal frame ranges are half-open [FRAME0, FRAME1).')
    h.add_history('Spectral weights are normalized; no per-image flux renormalization.')
    if e.pure:
        h.add_history('Static PSF repeated in time; frames are not independent realizations.')
    return h


def _write_atomic(path, data, header):
    partial = path.with_suffix('.partial.fits')
    fits.writeto(partial, np.asarray(data, dtype=np.float32), header, overwrite=True, checksum=True)
    os.replace(partial, path)


def _valid_output(path, plan, start, stop, expected_spatial=None):
    if not path.exists():
        return False
    try:
        with fits.open(path, memmap=True) as hdul:
            hdul.verify('exception')
            h = hdul[0].header
            shape = hdul[0].data.shape
            return (h.get('RUNHASH') == plan.fingerprint and h.get('COMPLETE') is True
                    and h.get('FRAME0') == start and h.get('FRAME1') == stop
                    and len(shape) == 3 and shape[0] == stop - start
                    and h.get('BITPIX') == -32
                    and (expected_spatial is None or shape[1:] == expected_spatial)
                    and 'CHECKSUM' in h and 'DATASUM' in h
                    and hdul[0].verify_checksum() == 1 and hdul[0].verify_datasum() == 1)
    except (OSError, ValueError, TypeError, KeyError):
        LOG.warning('Invalid FITS cache: %s', path)
        return False


@contextmanager
def _source_slice(plan, source, kind, start, stop):
    """Yield a read-only source interval without creating an input cache.

    The FITS mapping stays open through propagation. Only explicitly requested OPD
    unit conversion allocates a raw chunk; meter inputs and amplitude maps use views.
    """
    if source is None:
        yield None
        return
    src = _resolve(plan.input_root, source.filename)
    with fits.open(src, mode='readonly', memmap=True) as hdul:
        data = hdul[0].data
        static = kind == 'amplitude' and (data.ndim == 2 or data.shape[0] == 1)
        local_start, local_stop = start-source.global_start, stop-source.global_start
        if (stop <= start or local_start < 0 or
                (not static and (data.ndim != 3 or local_stop > data.shape[0]))):
            raise ValueError(f'Invalid {kind} frame range {start}:{stop} for {src}')
        cut = data if static else data[local_start:local_stop]
        if kind == 'phase' and plan.experiment.phase_unit != 'm':
            cut = np.asarray(cut, dtype=np.float64) * {'um': 1e-6, 'nm': 1e-9}[plan.experiment.phase_unit]
        if not np.all(np.isfinite(cut)):
            raise ValueError(f'Non-finite {kind} input in frames {start}:{stop}: {src}')
        yield cut


def _bounded_errors(conf, phase_data, amp_data, nframes, cpu_count):
    """Resample in-memory input slices using an explicit worker limit.

    HEEPS load_errors does not forward cpu_count to resize_cube. The single-core
    path below also avoids multiCPU's repeated stacking of growing arrays.
    """
    from heeps.util.img_processing import resize_cube, resize_img

    def resample(data):
        if data.ndim == 2:
            data = data[None]
        if cpu_count == 1 or len(data) == 1:
            first = resize_img(data[0], conf['npupil'])
            result = np.empty((len(data), *first.shape), dtype=first.dtype)
            result[0] = first
            for i in range(1, len(data)):
                result[i] = resize_img(data[i], conf['npupil'])
            return result
        result = resize_cube(data, conf['npupil'], cpu_count=min(cpu_count, len(data)))
        return np.array(result, ndmin=3)

    phases = resample(phase_data) if phase_data is not None else [None] * nframes
    amplitudes = resample(amp_data) if amp_data is not None else np.array([None] * nframes)
    tiptilts = np.array([None] * nframes)
    apo = np.array([conf['apo_misalign']] * nframes)
    stops = np.array([conf['ls_misalign']] * nframes)
    return phases, amplitudes, tiptilts, apo, stops


def _local_heeps_config(reader, **overrides):
    """Read installed HEEPS defaults without fetching its legacy example dataset.

    New readers can expose an explicit download switch. The published reader at
    c7dc898 instead calls a global extract_zip unconditionally when its example
    is absent. Run that reader with a private globals dictionary whose downloader
    is disabled; leave the installed module and all other callers untouched.
    """
    if 'download_inputs' in inspect.signature(reader).parameters:
        return reader(download_inputs=False, **overrides)
    if not isinstance(reader, FunctionType) or 'extract_zip' not in reader.__code__.co_names:
        raise RuntimeError('Unsupported HEEPS read_config: no download switch or recognized legacy downloader')

    def skip_example_download(*args, **kwargs):
        LOG.debug('Skipped HEEPS legacy example download; using configured local FITS inputs')

    def config_print(*args, **kwargs):
        # The original reader announces the download before calling extract_zip.
        if args and isinstance(args[0], str) and args[0].startswith('Downloading input files from Google Drive'):
            return
        print(*args, **kwargs)

    namespace = dict(reader.__globals__)
    namespace.update(extract_zip=skip_example_download, print=config_print)
    isolated_reader = FunctionType(reader.__code__, namespace, reader.__name__,
                                   reader.__defaults__, reader.__closure__)
    isolated_reader.__kwdefaults__ = reader.__kwdefaults__
    return isolated_reader(**overrides)


def heeps_backend(plan: Plan, start: int, stop: int):
    """Return (cube, pixel_scale_mas, actual_wavelengths); use one static frame for pure."""
    import heeps

    e = plan.experiment
    phase_source = _source_at(e.phase, start) if e.phase else None
    amp_source = _source_at(e.amplitude, start) if e.amplitude else None
    phase_path = str(_resolve(plan.input_root, phase_source.filename)) if phase_source else ''
    amp_path = str(_resolve(plan.input_root, amp_source.filename)) if amp_source else ''
    total = None
    actual_lams = []
    spatial = None
    pscale = None
    with _source_slice(plan, phase_source, 'phase', start, stop) as phase, \
            _source_slice(plan, amp_source, 'amplitude', start, stop) as amp:
        for lam, weight in zip(e.spectrum.wavelengths_m, e.spectrum.weights):
            conf = _local_heeps_config(heeps.config.read_config,
                dir_current=str(plan.run_dir), dir_input=str(plan.input_root),
                dir_output=str(plan.run_dir / 'heeps'), dir_temp=str(plan.run_dir / 'temp'),
                f_pupil=e.pupil, f_phase=phase_path, f_amp=amp_path,
                mode='ELT', ngrid=e.ngrid, nframes=1 if e.pure else stop-start,
                cpu_count=e.cpu_count, dit=e.dit_s, nstep=1, send_to=None,
                **({'select_lyot': '', 'f_lyot_stop': str(_resolve(plan.input_root, e.lyot_stop))}
                   if e.lyot_stop is not None else {}))
            band_spec = {**conf['band_specs'][e.band], 'lam': lam}
            if e.pscale_mas is not None:
                band_spec['pscale'] = e.pscale_mas
            conf.update(band=e.band, mask_name=e.mask_name, hfov=e.hfov_arcsec,
                        band_specs={**conf['band_specs'], e.band: band_spec},
                        add_phase=bool(e.phase), add_amp=bool(e.amplitude),
                        add_point_err=False, add_apo_drift=False, add_ls_drift=False,
                        add_disp=False, add_cl_det=False, add_cl_vort=False)
            conf = heeps.config.update_config(saveconf=False, verbose=False, **conf)
            if not 0 < conf['npupil'] <= e.ngrid or not 0 < conf['ndet'] <= e.ngrid:
                raise ValueError(f'Incompatible grid: npupil={conf["npupil"]}, ndet={conf["ndet"]}, ngrid={e.ngrid}')
            if not Path(conf['f_lyot_stop']).is_file():
                raise FileNotFoundError(conf['f_lyot_stop'])
            wf = heeps.pupil.pupil(savefits=False, verbose=False, **conf)
            errors = _bounded_errors(conf, phase, amp, 1 if e.pure else stop-start, e.cpu_count)
            cube = np.asarray(heeps.wavefront.propagate_cube(
                wf, *errors, onaxis=True, avg=False, savefits=False, verbose=False, **conf), dtype=np.float32)
            del errors
            if cube.ndim == 2:
                cube = cube[np.newaxis]
            wanted_frames = 1 if e.pure else stop-start
            if cube.shape != (wanted_frames, conf['ndet'], conf['ndet']):
                raise ValueError(f'Unexpected HEEPS output shape: {cube.shape}')
            if not np.all(np.isfinite(cube)) or np.any(cube < 0):
                raise ValueError('PSF must be finite and non-negative')
            if spatial is not None and (spatial != cube.shape[1:] or pscale != conf['pscale']):
                raise ValueError('Wavelengths produced incompatible spatial sampling')
            spatial, pscale = cube.shape[1:], float(conf['pscale'])
            actual_lams.append(float(conf['lam']))
            cube *= weight
            if total is None:
                total = cube
            else:
                total += cube
            del cube, wf
    return total, pscale, tuple(actual_lams)


def _chunk_path(plan, start, stop):
    return plan.run_dir / 'chunks' / f'psf_{start:08d}_{stop:08d}.fits'


def _manifest(plan):
    return dict(schema=SCHEMA, fingerprint=plan.fingerprint, experiment=asdict(plan.experiment),
                inputs=plan.inputs, code_hash=plan.code_hash, runtime=_runtime_versions(), chunks=plan.chunks)


def run_experiment(plan: Plan, *, backend: Callable = heeps_backend) -> list[Path]:
    """Compute missing chunks; refuse invalid existing results instead of overwriting.

    A backend argument permits small, deterministic integration tests without PROPER.
    Production calls should use the default backend.
    """
    _assert_current(plan)
    with _run_lock(plan):
        manifest = plan.run_dir / 'manifest.json'
        content = _manifest(plan)
        if manifest.exists() and json.loads(manifest.read_text()) != json.loads(json.dumps(content)):
            raise RuntimeError('Existing manifest differs from this plan')
        _atomic_json(manifest, content)
        (plan.run_dir / 'chunks').mkdir(exist_ok=True)
        results = []
        static_result = None
        for start, stop in plan.chunks:
            _assert_current(plan)
            path = _chunk_path(plan, start, stop)
            if _valid_output(path, plan, start, stop):
                LOG.info('Reuse frames %s:%s', start, stop)
                results.append(path)
                continue
            if path.exists():
                raise RuntimeError(f'Existing chunk failed validation; inspect or remove it: {path}')
            LOG.info('Compute %s frames %s:%s', plan.experiment.name, start, stop)
            if plan.experiment.pure and static_result is not None:
                cube, pscale, actual_lams = static_result
            else:
                cube, pscale, actual_lams = backend(plan, start, stop)
                cube = np.asarray(cube, dtype=np.float32)
                if cube.ndim == 2:
                    cube = cube[None]
                if plan.experiment.pure:
                    static_result = (cube, pscale, actual_lams)
            nframes = 1 if plan.experiment.pure else stop-start
            if cube.ndim != 3 or cube.shape[0] != nframes or min(cube.shape[1:]) < 1:
                raise ValueError(f'Unexpected backend cube shape: {cube.shape}')
            if not np.all(np.isfinite(cube)) or np.any(cube < 0):
                raise ValueError('Backend returned invalid PSF intensity')
            if (not np.isfinite(pscale) or pscale <= 0 or len(actual_lams) != len(plan.experiment.spectrum.weights)
                    or not np.all(np.isfinite(actual_lams)) or np.any(np.asarray(actual_lams) <= 0)):
                raise ValueError('Backend returned invalid sampling metadata')
            header = _header(plan, start, stop, pscale, actual_lams)
            header['CRPIX1'], header['CRPIX2'] = (cube.shape[2]+1)/2, (cube.shape[1]+1)/2
            shape = (stop-start, *cube.shape[1:])
            needed = math.prod(shape) * 4 + 28800
            if shutil.disk_usage(plan.run_dir).free < needed:
                raise OSError(f'Insufficient free space for chunk: need at least {needed} bytes')
            # A broadcast view avoids materializing a repeated pure cube in RAM.
            data = np.broadcast_to(cube, shape) if plan.experiment.pure else cube
            _assert_current(plan)
            _write_atomic(path, data, header)
            if not _valid_output(path, plan, start, stop):
                raise RuntimeError(f'Written chunk failed verification: {path}')
            results.append(path)
            _atomic_json(plan.run_dir / 'progress.json', dict(completed=[p.name for p in results],
                                                           last_stop=stop, status='running'))
            del data, cube
        _atomic_json(plan.run_dir / 'progress.json', dict(completed=[p.name for p in results],
                                                       last_stop=plan.experiment.stop, status='chunks_complete'))
        return results


def _empty_fits(path, shape, template):
    h = fits.PrimaryHDU(data=np.empty((1, 1, 1), dtype=np.float32)).header
    h['NAXIS1'], h['NAXIS2'], h['NAXIS3'] = shape[2], shape[1], shape[0]
    h.extend(template, update=True, strip=True)
    for key in ('CHECKSUM', 'DATASUM', 'BSCALE', 'BZERO', 'BLANK'):
        h.remove(key, ignore_missing=True)
    header_bytes = h.tostring().encode('ascii')
    data_bytes = math.prod(shape) * 4
    with path.open('wb') as stream:
        stream.write(header_bytes)
        stream.seek(len(header_bytes) + ((data_bytes + 2879)//2880)*2880 - 1)
        stream.write(b'\0')
    return np.memmap(path, dtype='>f4', mode='r+', offset=len(header_bytes), shape=shape)


def _sampling(header):
    keys = ['PSCALE', 'DIT', 'BAND', 'MASK', 'PURE', 'NLAM', 'NORM']
    keys += [f'WAVE{i} {field}' for i in range(header['NLAM']) for field in ('REQ', 'ACT', 'WEIGHT')]
    return tuple(header[k] for k in keys)


def merge_chunks(plan: Plan) -> Path:
    """Validate the full expected chunk set, then merge using bounded memory."""
    _assert_current(plan)
    with _run_lock(plan):
        e = plan.experiment
        if _valid_output(plan.full_file, plan, e.start, e.stop):
            _atomic_json(plan.run_dir / 'progress.json', dict(status='complete', output=str(plan.full_file)))
            return plan.full_file
        if plan.full_file.exists():
            raise RuntimeError(f'Existing full cube failed validation: {plan.full_file}')
        template = spatial = sampling = None
        paths = []
        for start, stop in plan.chunks:
            path = _chunk_path(plan, start, stop)
            if not _valid_output(path, plan, start, stop, spatial):
                raise ValueError(f'Missing or invalid chunk: {path}')
            with fits.open(path, memmap=True) as hdul:
                header = hdul[0].header.copy()
                if template is None:
                    spatial = hdul[0].data.shape[1:]
                    sampling = _sampling(header)
                    template = header
                elif _sampling(header) != sampling:
                    raise ValueError(f'Incompatible chunk sampling: {path}')
            paths.append((start, stop, path))
        shape = (e.stop-e.start, *spatial)
        needed = math.prod(shape) * 4 + 28800
        if shutil.disk_usage(plan.run_dir).free < needed:
            raise OSError(f'Insufficient free space for merge: need at least {needed} bytes')
        template['FRAME0'], template['FRAME1'] = e.start, e.stop
        template['CRVAL3'] = e.start * e.dit_s
        partial = plan.run_dir / 'full.partial.fits'
        out = _empty_fits(partial, shape, template)
        try:
            for start, stop, path in paths:
                with fits.open(path, memmap=True) as hdul:
                    # Copy in small blocks even if a caller chooses large chunks.
                    for offset in range(0, stop-start, 16):
                        block = hdul[0].data[offset:offset+16]
                        dest = start-e.start+offset
                        out[dest:dest+len(block)] = block
            out.flush()
        finally:
            del out
        with fits.open(partial, mode='update', memmap=True) as hdul:
            hdul[0].add_checksum()
            hdul.flush(output_verify='exception')
        _assert_current(plan)
        if not _valid_output(partial, plan, e.start, e.stop, spatial):
            raise RuntimeError(f'Merged FITS verification failed: {partial}')
        os.replace(partial, plan.full_file)
        _atomic_json(plan.run_dir / 'progress.json', dict(status='complete', output=str(plan.full_file)))
        return plan.full_file


# -----------------------------------------------------------------------------
# YAML configuration and command-line entry point
# -----------------------------------------------------------------------------
class UniqueKeyLoader(yaml.SafeLoader):
    """Reject duplicate keys at every YAML nesting level."""


def _unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str):
            raise ValueError('YAML keys must be strings')
        if key in result:
            raise ValueError(f'Duplicate YAML key: {key!r}')
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def _path(value, base):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Paths must be non-empty strings')
    expanded = os.path.expandvars(value)
    if '$' in expanded:
        raise ValueError(f'Unresolved environment variable in path: {value}')
    return (base / Path(expanded).expanduser()).resolve()


def _positive(value, name, *, integer=False, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, int if integer else Real):
        raise ValueError(f'{name} must be a finite {"integer" if integer else "number"}')
    if not np.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f'{name} must be {"non-negative" if allow_zero else "positive"}')


def _sources(value, base, name):
    if not isinstance(value, list):
        raise ValueError(f'{name} must be a list; use [] to disable this input')
    result = []
    for item in value:
        if isinstance(item, str):
            item = {'path': item}
        if not isinstance(item, dict) or set(item) - {'path', 'global_start'} or 'path' not in item:
            raise ValueError(f'{name} entries require path and optional global_start')
        result.append(Source(str(_path(item['path'], base)), item.get('global_start', 0)))
    return tuple(result)


@dataclass(frozen=True)
class PipelineConfig:
    name: str
    pupil_file: Path
    mask_file: Path
    output_dir: Path
    spectrum: Spectrum
    band: str
    mask_name: str
    cadence_s: float
    phase_files: tuple[Source, ...] = ()
    amplitude_files: tuple[Source, ...] = ()
    phase_unit: str = 'm'
    input_root: Path | None = None
    frame_start: int = 0
    nframes: int | None = None
    chunk_size: int = 100
    cpu_count: int = 1
    ngrid: int = 1024
    hfov_arcsec: float = 1.1
    pixel_scale_mas: float = 5.47
    pure: bool = False
    run_mode: str = 'preview'
    preview_frames: int = 2

    def __post_init__(self):
        for name in ('pupil_file', 'mask_file', 'output_dir'):
            object.__setattr__(self, name, Path(getattr(self, name)).expanduser().resolve())
        root = self.pupil_file.parent if self.input_root is None else Path(self.input_root).expanduser().resolve()
        object.__setattr__(self, 'input_root', root)
        for name in ('frame_start', 'chunk_size', 'cpu_count', 'ngrid', 'preview_frames'):
            _positive(getattr(self, name), name, integer=True, allow_zero=name == 'frame_start')
        if self.nframes is not None:
            _positive(self.nframes, 'nframes', integer=True)
        for name in ('cadence_s', 'hfov_arcsec', 'pixel_scale_mas'):
            _positive(getattr(self, name), name)
        if self.run_mode not in ('preview', 'full'):
            raise ValueError('run_mode must be preview or full')
        # Validate source ordering, error flags, units, band/mask and name immediately.
        self.experiment(stop=self.frame_start + (self.nframes or 1))

    def experiment(self, *, stop):
        return Experiment(
            name=self.name, band=self.band, mask_name=self.mask_name, spectrum=self.spectrum,
            phase=self.phase_files, amplitude=self.amplitude_files, pupil=str(self.pupil_file),
            lyot_stop=str(self.mask_file), phase_unit=self.phase_unit,
            start=self.frame_start, stop=stop, chunk_size=self.chunk_size,
            cpu_count=self.cpu_count, ngrid=self.ngrid, hfov_arcsec=self.hfov_arcsec,
            pscale_mas=self.pixel_scale_mas, dit_s=self.cadence_s, pure=self.pure)


def _spectrum(data):
    if not isinstance(data, dict):
        raise ValueError('spectrum must be a mapping')
    explicit = 'wavelengths_m' in data
    allowed = {'wavelengths_m', 'weights'} if explicit else {'lambda_center_m', 'frac_bw', 'nlam'}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f'Unknown or mixed spectrum keys: {sorted(unknown)}')
    if explicit:
        wavelengths = data['wavelengths_m']
        weights = data.get('weights')
        if not isinstance(wavelengths, list) or not wavelengths:
            raise ValueError('wavelengths_m must be a non-empty list')
        if weights is None:
            weights = [1.] * len(wavelengths)
        if not isinstance(weights, list):
            raise ValueError('weights must be a list')
        for value in wavelengths:
            _positive(value, 'wavelengths_m')
        for value in weights:
            _positive(value, 'weights', allow_zero=True)
        return Spectrum(tuple(wavelengths), tuple(weights))
    if 'lambda_center_m' not in data:
        raise ValueError('spectrum requires lambda_center_m or wavelengths_m')
    center, bw, nlam = data['lambda_center_m'], data.get('frac_bw', 0.), data.get('nlam', 1)
    _positive(center, 'lambda_center_m')
    _positive(bw, 'frac_bw', allow_zero=True)
    _positive(nlam, 'nlam', integer=True)
    return Spectrum.uniform(center, bw, nlam)


def load_config(path):
    path = Path(path).expanduser().resolve()
    data = yaml.load(path.read_text(), Loader=UniqueKeyLoader)
    return config_from_mapping(data, path.parent)


def _band_specs(value):
    """Validate optional per-band overrides, without importing HEEPS."""
    if not isinstance(value, dict):
        raise ValueError('band_specs must be a mapping')
    if set(value) - set(BAND_SPECS):
        raise ValueError(f'Unknown band_specs bands: {sorted(set(value) - set(BAND_SPECS))}')
    result = {band: dict(spec) for band, spec in BAND_SPECS.items()}
    for band, override in value.items():
        if not isinstance(override, dict):
            raise ValueError(f'band_specs.{band} must be a mapping')
        unknown = set(override) - {'lambda_center_m', 'pixel_scale_mas', 'frac_bw', 'nlam'}
        if unknown:
            raise ValueError(f'Unknown band_specs.{band} keys: {sorted(unknown)}')
        result[band].update(override)
    for band, spec in result.items():
        _positive(spec['pixel_scale_mas'], f'band_specs.{band}.pixel_scale_mas')
        _spectrum({key: value for key, value in spec.items() if key != 'pixel_scale_mas'})
    return result


def config_from_mapping(data, base):
    """Validate configuration data; resolve all paths against an explicit base."""
    base = Path(base).expanduser().resolve()
    if not isinstance(data, dict):
        raise ValueError('The YAML document must contain a configuration mapping')
    data = dict(data)
    allowed = {f.name for f in fields(PipelineConfig)} | {'band_specs'}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f'Unknown configuration keys: {sorted(unknown)}')
    required = {'name', 'pupil_file', 'mask_file', 'output_dir',
                'band', 'mask_name', 'cadence_s'}
    missing = required - set(data)
    if missing:
        raise ValueError(f'Missing required configuration keys: {sorted(missing)}')
    profiles = _band_specs(data.pop('band_specs', {}))
    if not isinstance(data['band'], str) or data['band'] not in profiles:
        raise ValueError(f'Unsupported band: {data["band"]!r}; choose from {list(profiles)}')
    profile = profiles[data['band']]
    if 'spectrum' not in data or data['spectrum'] == 'auto':
        data['spectrum'] = {key: value for key, value in profile.items() if key != 'pixel_scale_mas'}
    if 'pixel_scale_mas' not in data or data['pixel_scale_mas'] == 'auto':
        data['pixel_scale_mas'] = profile['pixel_scale_mas']
    for key in ('pupil_file', 'mask_file', 'output_dir', 'input_root'):
        if key in data and data[key] is not None:
            data[key] = _path(data[key], base)
    for key in ('pupil_file', 'mask_file', 'output_dir'):
        if data[key] is None:
            raise ValueError(f'{key} must be a path string')
    for key in ('phase_files', 'amplitude_files'):
        data[key] = _sources(data.get(key, []), base, key)
    if data['phase_files'] and 'phase_unit' not in data:
        raise ValueError('Declare phase_unit explicitly: m, um, or nm of optical path difference')
    data['spectrum'] = _spectrum(data['spectrum'])
    return PipelineConfig(**data)


def _available_stop(config):
    """Infer the common final frame, leaving gap/overlap checks to the planner."""
    ends = []
    for sources, kind in ((config.phase_files, 'phase'), (config.amplitude_files, 'amplitude')):
        if not sources:
            continue
        for source in sources:
            with fits.open(source.filename, memmap=True) as hdul:
                shape = hdul[0].shape
            static = kind == 'amplitude' and (len(shape) == 2 or (len(shape) == 3 and shape[0] == 1))
            if not static and (len(shape) != 3 or min(shape) < 1):
                raise ValueError(f'{kind} must be a non-empty 3D FITS cube: {source.filename}')
        if not static:
            ends.append(sources[-1].global_start + shape[0])
    if not ends:
        raise ValueError('Set nframes explicitly for pure or static-amplitude-only runs')
    stop = min(ends)
    if stop <= config.frame_start:
        raise ValueError('No input frames remain after frame_start')
    return stop


def inspect_inputs(config, *, work_dir=None, psf_file=None):
    """Validate local assets and build a read-only plan without importing HEEPS."""
    if not config.input_root.is_dir():
        raise FileNotFoundError(f'input_root directory not found: {config.input_root}')
    for path, name in ((config.pupil_file, 'pupil_file'), (config.mask_file, 'mask_file')):
        if not path.is_file():
            raise FileNotFoundError(f'{name} not found: {path}')
        with fits.open(path, memmap=True) as hdul:
            data = hdul[0].data
            if data.ndim != 2 or min(data.shape) < 1 or not np.all(np.isfinite(data)):
                raise ValueError(f'{name} must be a finite, non-empty 2D FITS image')
            if np.any(data < 0) or not np.any(data > 0):
                raise ValueError(f'{name} must contain non-negative amplitude transmission with positive support')
    # Equivalent to HEEPS round2odd(2 * hfov / pixel_scale * 1000) for positive sizes.
    ndet = 2 * int(config.hfov_arcsec / config.pixel_scale_mas * 1000) + 1
    if ndet > config.ngrid:
        raise ValueError(f'Detector size {ndet} exceeds ngrid={config.ngrid}; reduce hfov or increase ngrid')
    stop = config.frame_start + config.nframes if config.nframes is not None else _available_stop(config)
    experiment = config.experiment(stop=stop)
    # Validate the complete requested range even when executing only a preview.
    full_plan = plan_experiment(experiment, config.input_root, config.output_dir)
    plan = full_plan
    if config.run_mode == 'preview':
        plan = plan_experiment(experiment.preview(config.preview_frames), config.input_root, config.output_dir)
    if work_dir is not None or psf_file is not None:
        work = Path(work_dir).expanduser().resolve() if work_dir is not None else plan.run_dir
        output = Path(psf_file).expanduser().resolve() if psf_file is not None else None
        for destination in (work, output):
            if destination is not None and (destination == plan.input_root or plan.input_root in destination.parents):
                raise ValueError('Output must be outside input_root')
        plan = replace(plan, run_dir=work, output_file=output)
    return plan, dict(plan.summary(), requested_frames=stop-config.frame_start,
                      effective_frames=plan.experiment.stop-plan.experiment.start,
                      run_mode=config.run_mode, cadence_s=config.cadence_s,
                      phase_unit=config.phase_unit, phase_enabled=bool(config.phase_files),
                      amplitude_enabled=bool(config.amplitude_files),
                      pixel_scale_mas=config.pixel_scale_mas, psf_cube=str(plan.full_file),
                      output_shape=[plan.experiment.stop-plan.experiment.start, ndet, ndet],
                      cube_data_bytes=(plan.experiment.stop-plan.experiment.start)*ndet*ndet*4)


def config_dict(config):
    return json.loads(json.dumps(asdict(config), default=str))


def run_pipeline(config, config_path=None, *, work_dir=None, psf_file=None):
    plan, summary = inspect_inputs(config, work_dir=work_dir, psf_file=psf_file)
    run_experiment(plan)
    plan.full_file.parent.mkdir(parents=True, exist_ok=True)
    cube = merge_chunks(plan)
    # A small stable result pointer avoids duplicating multi-GB cubes or symlinking data.
    result = dict(status='complete', psf_cube=str(cube), fingerprint=plan.fingerprint,
                  config_path=str(Path(config_path).resolve()) if config_path else None,
                  completed_utc=datetime.now(timezone.utc).isoformat(),
                  requested_configuration=config_dict(config),
                  effective_experiment=asdict(plan.experiment), plan=summary)
    result_path = (plan.run_dir if work_dir is not None else plan.run_dir.parent) / 'result.json'
    tmp = result_path.with_name(f'.result.{uuid4().hex}.partial.json')
    try:
        tmp.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        os.replace(tmp, result_path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return cube, result_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path, help='YAML configuration; paths are relative to this file')
    parser.add_argument('--phase-file', type=Path, help='Override with one phase cube, global origin 0; path relative to the shell cwd')
    parser.add_argument('--amplitude-file', type=Path, help='Override with one amplitude input, global origin 0; path relative to the shell cwd')
    parser.add_argument('--output-dir', type=Path, help='Override output root; path relative to the shell cwd')
    parser.add_argument('--nframes', type=int, help='Override requested frame count; preview still limits execution')
    parser.add_argument('--cpu-count', type=int, help='Override workers; previews always use one worker')
    parser.add_argument('--dry-run', action='store_true', help='Validate configuration and local inputs, print plan; create no output or simulation')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--preview', action='store_true', help='Force preview mode')
    modes.add_argument('--full', action='store_true', help='Force full mode')
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        changes = {}
        for option, key in ((args.phase_file, 'phase_files'), (args.amplitude_file, 'amplitude_files')):
            if option is not None:
                changes[key] = (Source(str(option.expanduser().resolve())),)
        if args.output_dir is not None:
            changes['output_dir'] = args.output_dir.expanduser().resolve()
        for key in ('nframes', 'cpu_count'):
            if getattr(args, key) is not None:
                changes[key] = getattr(args, key)
        if args.preview or args.full:
            changes['run_mode'] = 'preview' if args.preview else 'full'
        config = replace(config, **changes)
        _, summary = inspect_inputs(config)
        print('Requested configuration:')
        print(yaml.safe_dump(config_dict(config), sort_keys=False))
        print('Execution plan:')
        print(json.dumps(summary, indent=2))
        if args.dry_run:
            return 0
        cube, result_path = run_pipeline(config, args.config)
        print(f'PSF_CUBE={cube}')
        print(f'RESULT_JSON={result_path}')
        return 0
    except KeyboardInterrupt:
        print('Interrupted; completed valid chunks can be resumed by rerunning the command.', file=sys.stderr)
        return 130
    except (ValueError, TypeError, OSError, RuntimeError, yaml.YAMLError) as exc:
        print(f'Configuration/input/run error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    raise SystemExit(main())
