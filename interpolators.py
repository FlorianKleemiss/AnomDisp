"""
interpolators.py
================
All interpolation infrastructure for the theta-integral I_integral(qr, l, m, typ, part).

Exports
-------
M_VALUES_BY_TYPE
get_m_values_for_type(typ)
get_parts_for_typ_m(typ, m)
_interpolator_key(n0, l, m, typ, part)
_build_interpolator(...)
_build_adaptive_interpolator(...)
I_integral(args)
I_interpolator
I_adaptive_interpolator
generate_interpolators_file(...)
load_interpolators_file(...)
enable_precomputed_interpolators(...)
generate_adaptive_interpolators_file(...)
load_adaptive_interpolators_file(...)
enable_precomputed_adaptive_interpolators(...)
init_precomputed_interpolators(...)
_pool_worker_init(...)
test_precomputed_interpolators(...)
test_precomputed_adaptive_interpolators(...)
resolve_interpolators_file(...)
PRECOMPUTED_INTERPOLATORS
PRECOMPUTED_INTERPOLATORS_ADAPTIVE
ACTIVE_INTERPOLATORS_FILE
DEFAULT_INTERPOLATOR_FILES
"""

import bisect
import json
import os
import math
import re

import mpmath
import numpy as np
from scipy.special import lpmv

digits = 14
mpmath.mp.dps = digits

# ---------------------------------------------------------------------------
# Constants / type mappings
# ---------------------------------------------------------------------------

M_VALUES_BY_TYPE = {
    "s": (1,),
    "p": (0, 1, 2),
    "d": (0, 1, 2, 3, 4),
}

DEFAULT_INTERPOLATOR_FILES = ("interpolator.json", "interpolators.json")
DEFAULT_CUTOFF_FILE = "interpolator_cutoffs_reasonable.json"

# ---------------------------------------------------------------------------
# Module-level cache globals (set by enable_* / init_* functions)
# ---------------------------------------------------------------------------

PRECOMPUTED_INTERPOLATORS = None
PRECOMPUTED_INTERPOLATORS_ADAPTIVE = None
ACTIVE_INTERPOLATORS_FILE = None

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def resolve_interpolators_file(input_file=None):
    if input_file is not None and os.path.exists(input_file):
        return input_file
    for candidate in DEFAULT_INTERPOLATOR_FILES:
        if os.path.exists(candidate):
            return candidate
    return input_file if input_file is not None else DEFAULT_INTERPOLATOR_FILES[0]


def get_m_values_for_type(typ):
    if typ not in M_VALUES_BY_TYPE:
        raise ValueError("Unknown typ={}".format(typ))
    return M_VALUES_BY_TYPE[typ]


def get_m_values_for_type_l(typ, l):
    """Return only physically meaningful m values for a given l (m <= l)."""
    l_int = int(l)
    return tuple(m for m in get_m_values_for_type(typ) if m <= l_int)


def get_parts_for_typ_m(typ, m):
    if typ == "s":
        return [1]
    if typ == "p":
        if m == 0:
            return [1, 2]
        if m in (1, 2):
            return [1]
    if typ == "d":
        if m in (0, 1, 4):
            return [1, 2]
        if m in (2, 3):
            return [1]
    raise ValueError("Unsupported typ={}, m={} combination".format(typ, m))


def _interpolator_key(n0, l, m, typ, part):
    return "n0={}|l={}|m={}|typ={}|part={}".format(n0, l, m, typ, part)


def _integral_cache_key(n0, l, m, typ, part):
    return "n0={}|l={}|m={}|typ={}|part={}".format(int(n0), int(l), int(m), str(typ), int(part))


def _integral_cache_file_path(cache_dir, n0, l, m, typ, part):
    key = _integral_cache_key(n0, l, m, typ, part)
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
    return os.path.join(str(cache_dir), safe + ".json")


class _ChannelIntegralPointStore:
    """Persistent on-disk cache for per-channel I_integral(qr) samples.

    The cache is append/merge style: existing points are loaded and preserved,
    new points are added, and flush writes the merged set back to disk.
    """

    def __init__(self, cache_dir, n0, l, m, typ, part, enabled=True):
        self.enabled = bool(enabled)
        self.n0 = int(n0)
        self.l = int(l)
        self.m = int(m)
        self.typ = str(typ)
        self.part = int(part)
        self.cache_dir = str(cache_dir)
        self.file_path = _integral_cache_file_path(self.cache_dir, self.n0, self.l, self.m, self.typ, self.part)
        self.points = {}
        self._dirty = False
        if self.enabled:
            os.makedirs(self.cache_dir, exist_ok=True)
            self._load()

    @staticmethod
    def _qkey(q):
        qf = float(q)
        if qf == 0.0:
            qf = 0.0
        return format(qf, ".17g")

    def _load(self):
        if not os.path.exists(self.file_path):
            return
        try:
            with open(self.file_path, "r") as inp:
                raw = json.load(inp)
            pts = raw.get("points", []) if isinstance(raw, dict) else []
            for row in pts:
                if not isinstance(row, (list, tuple)) or len(row) < 3:
                    continue
                q = float(row[0])
                r = float(row[1])
                im = float(row[2])
                self.points[self._qkey(q)] = (q, complex(r, im))
        except Exception:
            # Keep generation robust even if cache file is malformed.
            self.points = {}

    def evaluate(self, q):
        qf = float(q)
        k = self._qkey(qf)
        if k in self.points:
            return self.points[k][1]
        v = _I_integral_direct((qf, self.l, self.m, self.typ, self.part))
        if self.enabled:
            self.points[k] = (qf, complex(v))
            self._dirty = True
        return v

    def get_points_in_range(self, q_min, q_max):
        q0 = float(min(q_min, q_max))
        q1 = float(max(q_min, q_max))
        out = {}
        for q, v in self.points.values():
            if q0 <= float(q) <= q1:
                out[float(q)] = complex(v)
        return out

    def flush(self):
        if not self.enabled or (not self._dirty):
            return

        merged = dict(self.points)
        if os.path.exists(self.file_path):
            try:
                with open(self.file_path, "r") as inp:
                    raw = json.load(inp)
                pts = raw.get("points", []) if isinstance(raw, dict) else []
                for row in pts:
                    if not isinstance(row, (list, tuple)) or len(row) < 3:
                        continue
                    q = float(row[0])
                    r = float(row[1])
                    im = float(row[2])
                    merged[self._qkey(q)] = (q, complex(r, im))
            except Exception:
                pass

        rows = []
        for _, (q, v) in sorted(merged.items(), key=lambda item: item[1][0]):
            rows.append([float(q), float(complex(v).real), float(complex(v).imag)])

        payload = {
            "channel_key": _integral_cache_key(self.n0, self.l, self.m, self.typ, self.part),
            "points": rows,
        }
        tmp_path = self.file_path + ".tmp"
        with open(tmp_path, "w") as out:
            json.dump(payload, out)
        os.replace(tmp_path, self.file_path)

        self.points = merged
        self._dirty = False


def _load_cutoff_recommendations_file(input_file=DEFAULT_CUTOFF_FILE):
    if input_file is None or (not os.path.exists(input_file)):
        return {}
    with open(input_file, "r") as inp:
        raw = json.load(inp)
    return raw.get("cutoffs", {})


def _resolve_channel_qr_max(
    key,
    global_qr_max,
    cutoff_data=None,
    cutoff_field="cutoffs_relaxed_5pct",
    cutoff_threshold="0.2",
    cutoff_pad_factor=1.02,
    min_qr_max=1.0,
):
    qmax = float(global_qr_max)
    if not cutoff_data or key not in cutoff_data:
        return qmax

    entry = cutoff_data.get(key, {})
    field_data = entry.get(cutoff_field, {}) if isinstance(entry, dict) else {}
    if not isinstance(field_data, dict):
        return qmax

    if cutoff_threshold in field_data:
        raw_cut = field_data.get(cutoff_threshold)
    else:
        # Fallback for formatting mismatches like "0.2" vs "0.200000".
        raw_cut = field_data.get(str(float(cutoff_threshold)), qmax)

    try:
        cut = float(raw_cut)
    except Exception:
        return qmax

    if cut <= 0.0:
        return max(float(min_qr_max), min(float(global_qr_max), float(min_qr_max)))

    padded = cut * float(cutoff_pad_factor)
    return max(float(min_qr_max), min(float(global_qr_max), float(padded)))


def _n0_typ_valid(n0, typ):
    """Return True if this (n0, typ) combination is physically available."""
    if n0 == 1 and typ != "s":
        return False
    if n0 == 2 and typ == "d":
        return False
    if n0 == 4 and typ == "d":
        return False
    if n0 == 5 and typ != "s":
        return False
    return True


def _extrema_indices_from_series(y):
    """Return indices of local extrema for a 1D sampled series."""
    arr = np.asarray(y, dtype=float)
    if arr.size < 5:
        return []
    dy = np.diff(arr)
    # Replace exact zeros with tiny signed values from nearest non-zero slope.
    for i in range(dy.size):
        if dy[i] == 0.0:
            left = i - 1
            right = i + 1
            repl = 0.0
            while left >= 0:
                if dy[left] != 0.0:
                    repl = dy[left]
                    break
                left -= 1
            if repl == 0.0:
                while right < dy.size:
                    if dy[right] != 0.0:
                        repl = dy[right]
                        break
                    right += 1
            dy[i] = repl
    idx = []
    for i in range(1, dy.size):
        if (dy[i - 1] > 0.0 and dy[i] < 0.0) or (dy[i - 1] < 0.0 and dy[i] > 0.0):
            idx.append(i)
    return idx


def _probe_first_cycles_range(
    l,
    m,
    typ,
    part,
    qr_max=4500.0,
    n_probe=120,
    fit_num_cycles=3,
    scan_qr_max=120.0,
    scan_points=240,
    eval_fn=None,
):
    """Probe only the first oscillation cycles and return a compact effective range.

    This is a fast path for fit-based models where we train on early minima/maxima
    cycles and extrapolate the tail analytically.
    """
    fit_num_cycles = max(int(fit_num_cycles), 1)
    q_scan_max = max(10.0, min(float(qr_max), float(scan_qr_max)))
    n_scan = max(int(scan_points), 80)

    eval_q = eval_fn if eval_fn is not None else (lambda q: I_integral((q, l, m, typ, part)))
    scan_q = [float(v) for v in np.linspace(0.0, q_scan_max, num=n_scan)]
    scan_map = {q: eval_q(q) for q in scan_q}

    probe_max = min(10.0, q_scan_max)
    re_sq = 0.0
    im_sq = 0.0
    for q, v in scan_map.items():
        if q <= probe_max:
            c = complex(v)
            re_sq += c.real ** 2
            im_sq += c.imag ** 2
    use_real = re_sq >= im_sq
    y = [complex(scan_map[q]).real if use_real else complex(scan_map[q]).imag for q in scan_q]
    ext_idx = _extrema_indices_from_series(y)

    needed_extrema = 2 * fit_num_cycles
    if len(ext_idx) >= needed_extrema:
        end_idx = ext_idx[needed_extrema - 1]
    elif len(ext_idx) > 0:
        end_idx = ext_idx[-1]
    else:
        end_idx = len(scan_q) - 1

    end_idx = min(max(end_idx, 8), len(scan_q) - 1)
    effective_max = float(scan_q[end_idx])
    effective_max = max(effective_max, min(q_scan_max, 8.0))

    n_probe = max(int(n_probe), 80)
    log_start = max(effective_max * 1e-5, 1e-3)
    log_pts = list(np.geomspace(log_start, effective_max, num=n_probe))
    lin_pts = list(np.linspace(0.0, effective_max, num=n_probe))
    probe_pts = sorted(set([0.0] + [float(v) for v in log_pts + lin_pts]))

    value_map = dict(scan_map)
    for q in probe_pts:
        if q not in value_map:
            value_map[q] = eval_q(q)

    filtered = {q: v for q, v in value_map.items() if q <= effective_max * 1.001}
    return effective_max, filtered


def _theta_nodes_for_qr(qr, qr_split_start=40.0, panels_per_oscillation=3.0, max_panels=2400):
    """Return theta breakpoints for robust oscillatory integration.

    For high |qr| the phase exp(-i*qr*cos(theta)) oscillates rapidly.  We split
    [0, pi] into qr-dependent panels so mpmath.quad can resolve cancellation
    reliably instead of integrating the whole interval in one block.
    """
    q_abs = abs(float(qr))
    if q_abs <= float(qr_split_start):
        return [mpmath.mpf(0.0), mpmath.pi]

    # Around theta=pi/2, the local oscillation count over [0, pi] is ~|qr|/2.
    # Use multiple panels per oscillation for stable high-q cancellation.
    est_osc = 0.5 * q_abs
    n_panels = int(math.ceil(est_osc * float(panels_per_oscillation)))
    n_panels = max(8, min(int(max_panels), n_panels))

    step = mpmath.pi / n_panels
    return [step * i for i in range(n_panels + 1)]


# ---------------------------------------------------------------------------
# Core theta-integral
# ---------------------------------------------------------------------------

def _I_integral_direct(args):
    qr, l, m, typ, part = args

    def integrant(theta):
        ct = mpmath.cos(theta)
        oszillating_stuff = mpmath.exp(-1j * qr * ct)

        if typ == "s":
            st = mpmath.sin(theta) ** 2
            p_lm = float(lpmv(1, int(l), float(ct)))
            return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
        elif typ == "p":
            if m == 0:
                if part == 1:
                    st = mpmath.sin(theta)
                    p_lm = float(lpmv(0, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
                elif part == 2:
                    st = mpmath.sin(theta) ** 3
                    p_lm = float(lpmv(0, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 1:
                st = mpmath.sin(theta) ** 2 * ct
                p_lm = float(lpmv(1, int(l), float(ct)))
                return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 2:
                st = mpmath.sin(theta) ** 3
                p_lm = float(lpmv(2, int(l), float(ct)))
                return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
        elif typ == "d":
            if m == 0:
                if part == 1:
                    st = mpmath.sin(theta) * ct
                    p_lm = float(lpmv(0, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
                elif part == 2:
                    st = mpmath.sin(theta) ** 3 * ct
                    p_lm = float(lpmv(0, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 1:
                if part == 1:
                    st = mpmath.sin(theta) ** 2
                    p_lm = float(lpmv(1, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
                elif part == 2:
                    st = mpmath.sin(theta) ** 4
                    p_lm = float(lpmv(1, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 2:
                st = mpmath.sin(theta) ** 3 * ct
                p_lm = float(lpmv(2, int(l), float(ct)))
                return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 3:
                st = mpmath.sin(theta) ** 4
                p_lm = float(lpmv(3, int(l), float(ct)))
                return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
            elif m == 4:
                if part == 1:
                    st = mpmath.sin(theta) ** 2
                    p_lm = float(lpmv(1, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff
                elif part == 2:
                    st = mpmath.sin(theta) ** 2 * ct ** 2
                    p_lm = float(lpmv(1, int(l), float(ct)))
                    return st * mpmath.mpf(p_lm) * (-1.0) ** m * oszillating_stuff

    theta_nodes = _theta_nodes_for_qr(qr)
    result = mpmath.quad(integrant, theta_nodes)
    return result


def I_integral(args):
    return _I_integral_direct(args)


# ---------------------------------------------------------------------------
# Linear (uniform-grid) interpolator
# ---------------------------------------------------------------------------

class I_interpolator:
    def __init__(self, n0, l, m, typ, part, precomputed_cache=None, qr_values=None):
        self.n0 = n0
        self.l = l
        self.m = m
        self.typ = typ
        self.part = part
        self.distance_cache_mpf = []
        self.value_cache_real_mpf = []
        self.value_cache_imag_mpf = []
        self.interpolator_real_mpf = None
        self.interpolator_imag_mpf = None
        if precomputed_cache is None:
            self.make_cache(qr_values=qr_values)
        else:
            self.load_cache(precomputed_cache)

    def __call__(self, qr):
        return self.interpolate(qr)

    def interpolate(self, qr):
        if self.interpolator_real_mpf is None or self.interpolator_imag_mpf is None:
            raise ValueError("Interpolator not initialized. Call make_cache() first.")
        real_part = float(self.interpolator_real_mpf(float(qr)))
        imag_part = float(self.interpolator_imag_mpf(float(qr)))
        return mpmath.mpc(real_part, imag_part)

    def interpolate_mpf(self, qr):
        if self.interpolator_real_mpf is None or self.interpolator_imag_mpf is None:
            raise ValueError("MPF interpolator not initialized. Call make_cache() first.")
        real_part = self.interpolator_real_mpf(qr)
        imag_part = self.interpolator_imag_mpf(qr)
        return mpmath.mpc(real_part, imag_part)

    def make_linear_mpf_interpolator(self, x_vals, y_vals):
        return lambda x: (
            y_vals[0]
            if x <= x_vals[0]
            else y_vals[-1]
            if x >= x_vals[-1]
            else (
                y_vals[bisect.bisect_left(x_vals, x) - 1]
                + (mpmath.mpf(x) - x_vals[bisect.bisect_left(x_vals, x) - 1])
                / (x_vals[bisect.bisect_left(x_vals, x)] - x_vals[bisect.bisect_left(x_vals, x) - 1])
                * (y_vals[bisect.bisect_left(x_vals, x)] - y_vals[bisect.bisect_left(x_vals, x) - 1])
            )
        )

    def set_interpolators_from_cache(self):
        if all(v == mpmath.mpf(0) for v in self.value_cache_real_mpf):
            self.interpolator_real_mpf = lambda x: mpmath.mpf("0")
        else:
            self.interpolator_real_mpf = self.make_linear_mpf_interpolator(
                self.distance_cache_mpf, self.value_cache_real_mpf
            )
        if all(v == mpmath.mpf(0) for v in self.value_cache_imag_mpf):
            self.interpolator_imag_mpf = lambda x: mpmath.mpf("0")
        else:
            self.interpolator_imag_mpf = self.make_linear_mpf_interpolator(
                self.distance_cache_mpf, self.value_cache_imag_mpf
            )

    def load_cache(self, precomputed_cache):
        self.distance_cache_mpf = [mpmath.mpf(v) for v in precomputed_cache["distance_cache"]]
        self.value_cache_real_mpf = [mpmath.mpf(v) for v in precomputed_cache["value_cache_real"]]
        self.value_cache_imag_mpf = [mpmath.mpf(v) for v in precomputed_cache["value_cache_imag"]]
        self.set_interpolators_from_cache()

    def make_cache(self, qr_values=None):
        if qr_values is None:
            qr_values = np.linspace(0, 4500, num=4500)
        value_cache = np.array([I_integral((qr, self.l, self.m, self.typ, self.part)) for qr in qr_values])
        self.distance_cache_mpf = [mpmath.mpf(qr) for qr in qr_values]
        self.value_cache_real_mpf = [mpmath.mpf(v.real) for v in value_cache]
        self.value_cache_imag_mpf = [mpmath.mpf(v.imag) for v in value_cache]
        self.set_interpolators_from_cache()


# ---------------------------------------------------------------------------
# Adaptive PCHIP interpolator
# ---------------------------------------------------------------------------

class I_adaptive_interpolator:
    """Adaptive-grid real/imag PCHIP interpolator for I_integral.

    The adaptive grid is built by successively bisecting intervals whose
    midpoint linear-interpolation error exceeds rel_tol.  PCHIP (piecewise
    cubic Hermite, from scipy) is used for reconstruction, which is C1-
    continuous and avoids the zero-crossing / phase-jump artefacts of
    log-amp/phase representations.
    """

    def __init__(
        self,
        n0,
        l,
        m,
        typ,
        part,
        precomputed_cache=None,
        qr_min=0.0,
        qr_max=4500.0,
        rel_tol=1e-3,
        initial_points=24,
        max_points=8000,
        auto_qr_max=True,
        auto_threshold=1e-10,
        abs_tol=1e-10,
        low_q_switch=0.1,
        low_q_anchor_points=32,
        auto_low_q_switch=True,
        low_q_change_threshold=1e-4,
        phase_tol=0.35,
        use_q_scaling=True,
        phase_amp_rel_floor=1e-3,
        min_interval_rel=1e-5,
        min_interval_abs=1e-6,
        max_new_points_per_iter=128,
        zero_component_threshold=1e-2,
        use_tail_fit=False,
        tail_fit_start_q=25.0,
        fit_num_cycles=0,
        fit_scan_qr_max=120.0,
        fit_scan_points=240,
        use_point_cache=True,
        point_cache_dir="integral_point_cache",
    ):
        from scipy.interpolate import PchipInterpolator

        self._PchipInterpolator = PchipInterpolator
        self.n0 = n0
        self.l = l
        self.m = m
        self.typ = typ
        self.part = part
        self.distance_cache_mpf = []
        self.value_cache_real_mpf = []
        self.value_cache_imag_mpf = []
        self._interp_real = None
        self._interp_imag = None
        self._interp_low_q_real = None
        self._interp_low_q_imag = None
        self._distance_cache_float = []
        self.low_q_switch_effective = 0.0
        self.qr_max_effective = float(qr_max)
        self.build_stats = {}
        self.use_q_scaling = bool(use_q_scaling)
        self._dominant_component = "both"
        self._tail_fit = None
        if precomputed_cache is None:
            self.make_cache(
                qr_min=qr_min,
                qr_max=qr_max,
                rel_tol=rel_tol,
                initial_points=initial_points,
                max_points=max_points,
                auto_qr_max=auto_qr_max,
                auto_threshold=auto_threshold,
                abs_tol=abs_tol,
                low_q_switch=low_q_switch,
                low_q_anchor_points=low_q_anchor_points,
                auto_low_q_switch=auto_low_q_switch,
                low_q_change_threshold=low_q_change_threshold,
                phase_tol=phase_tol,
                use_q_scaling=use_q_scaling,
                phase_amp_rel_floor=phase_amp_rel_floor,
                min_interval_rel=min_interval_rel,
                min_interval_abs=min_interval_abs,
                max_new_points_per_iter=max_new_points_per_iter,
                zero_component_threshold=zero_component_threshold,
                use_tail_fit=use_tail_fit,
                tail_fit_start_q=tail_fit_start_q,
                fit_num_cycles=fit_num_cycles,
                fit_scan_qr_max=fit_scan_qr_max,
                fit_scan_points=fit_scan_points,
                use_point_cache=use_point_cache,
                point_cache_dir=point_cache_dir,
            )
        else:
            self.load_cache(precomputed_cache)

    def __call__(self, qr):
        return self.interpolate(qr)

    def _detect_dominant_component(self, value_map, probe_max=10.0, threshold=1e-2):
        """Check which component (real/imag) is significant in [0, probe_max].

        Returns 'real', 'imag', or 'both'.
        A component is treated as negligible when its RMS over the probe range
        is less than *threshold* times the total RMS.
        """
        re_sq = 0.0
        im_sq = 0.0
        for q, v in value_map.items():
            if q <= probe_max:
                c = complex(v)
                re_sq += c.real ** 2
                im_sq += c.imag ** 2
        total = re_sq + im_sq
        if total < 1e-60:
            return "both"
        thr_sq = float(threshold) ** 2
        if im_sq / total < thr_sq:
            return "real"
        if re_sq / total < thr_sq:
            return "imag"
        return "both"

    def _build_pchip(self):
        # Keep x strictly increasing for scipy PCHIP; cached merge points can
        # occasionally contain duplicate/unsorted q due to float round-trips.
        triples = []
        for q_mpf, r_mpf, i_mpf in zip(self.distance_cache_mpf, self.value_cache_real_mpf, self.value_cache_imag_mpf):
            qf = float(q_mpf)
            if not np.isfinite(qf):
                continue
            triples.append((qf, float(r_mpf), float(i_mpf)))
        triples.sort(key=lambda t: t[0])

        xs_list = []
        rs_list = []
        ims_list = []
        for qf, rv, iv in triples:
            if not xs_list or qf > xs_list[-1]:
                xs_list.append(qf)
                rs_list.append(rv)
                ims_list.append(iv)
            else:
                # Replace duplicate/non-increasing point with the latest value.
                xs_list[-1] = qf
                rs_list[-1] = rv
                ims_list[-1] = iv

        if len(xs_list) < 2:
            raise ValueError("Adaptive cache has fewer than 2 valid distinct q points.")

        # Persist the cleaned cache so all downstream logic sees a monotone grid.
        self.distance_cache_mpf = [mpmath.mpf(q) for q in xs_list]
        self.value_cache_real_mpf = [mpmath.mpf(v) for v in rs_list]
        self.value_cache_imag_mpf = [mpmath.mpf(v) for v in ims_list]

        xs = np.array(xs_list, dtype=float)
        self._distance_cache_float = list(xs_list)
        if self._dominant_component != "imag":
            rs = np.array([self._encode_value(v, x) for x, v in zip(xs_list, rs_list)], dtype=float)
            self._interp_real = self._PchipInterpolator(xs, rs, extrapolate=False)
        else:
            self._interp_real = None
        if self._dominant_component != "real":
            ims = np.array([self._encode_value(v, x) for x, v in zip(xs_list, ims_list)], dtype=float)
            self._interp_imag = self._PchipInterpolator(xs, ims, extrapolate=False)
        else:
            self._interp_imag = None

    @staticmethod
    def _damped_oscillator(q, alpha, omega, c_cos, c_sin):
        return np.exp(-alpha * q) * (c_cos * np.cos(omega * q) + c_sin * np.sin(omega * q))

    @staticmethod
    def _damped_oscillator_2h(q, alpha, omega, c1_cos, c1_sin, c2_cos, c2_sin):
        return np.exp(-alpha * q) * (
            c1_cos * np.cos(omega * q)
            + c1_sin * np.sin(omega * q)
            + c2_cos * np.cos(2.0 * omega * q)
            + c2_sin * np.sin(2.0 * omega * q)
        )

    @staticmethod
    def _damped_oscillator_ap(q, alpha, omega, A, phi):
        """Amplitude+phase form: exp(-alpha*q) * A * cos(omega*q + phi)."""
        return np.exp(-alpha * q) * A * np.cos(omega * q + phi)

    @staticmethod
    def _damped_oscillator_2h_ap(q, alpha, omega, A1, phi1, A2, phi2):
        """Two-harmonic amplitude+phase form."""
        return np.exp(-alpha * q) * (
            A1 * np.cos(omega * q + phi1)
            + A2 * np.cos(2.0 * omega * q + phi2)
        )

    def _fit_dominant_tail(self, tail_fit_start_q=25.0, fit_num_cycles=0):
        """Fit dominant component tail to a damped oscillatory model.

        Base model: exp(-alpha*q) * (c_cos*cos(omega*q) + c_sin*sin(omega*q))
        Enhanced model (optional): add second harmonic terms with 2*omega.
        """
        if self._dominant_component not in ("real", "imag"):
            self._tail_fit = None
            return

        xs = np.array([float(v) for v in self.distance_cache_mpf], dtype=float)
        if xs.size < 16:
            self._tail_fit = None
            return

        q_start = float(max(0.0, tail_fit_start_q))
        q_stop = float(self.qr_max_effective)
        mask = (xs >= q_start) & (xs <= q_stop)
        if np.count_nonzero(mask) < 16:
            self._tail_fit = None
            return

        if self._dominant_component == "real":
            ys_all = np.array([float(v) for v in self.value_cache_real_mpf], dtype=float)
        else:
            ys_all = np.array([float(v) for v in self.value_cache_imag_mpf], dtype=float)

        x = xs[mask]
        y = ys_all[mask]
        if y.size < 16:
            self._tail_fit = None
            return

        fit_num_cycles = max(int(fit_num_cycles), 0)
        if fit_num_cycles > 0:
            ext = _extrema_indices_from_series(y)
            needed = 2 * fit_num_cycles
            if len(ext) >= needed:
                stop = ext[needed - 1] + 1
                stop = min(max(stop, 12), len(x))
                x = x[:stop]
                y = y[:stop]
            elif len(ext) > 0:
                stop = ext[-1] + 1
                stop = min(max(stop, 12), len(x))
                x = x[:stop]
                y = y[:stop]

        amp_peak = float(np.max(np.abs(y)))
        if amp_peak < 1e-14:
            self._tail_fit = None
            return

        eps = 1e-18
        yabs = np.abs(y)
        valid = yabs > (amp_peak * 1e-6 + eps)
        if np.count_nonzero(valid) < 6:
            alpha0 = 0.02
            p0 = 1.0
        else:
            xv = x[valid]
            yv = np.log(yabs[valid] + eps)
            # Exponential initial guess: log|y| ~ -alpha*x
            slope_exp, _ = np.polyfit(xv, yv, 1)
            alpha0 = max(0.0, min(2.0, float(-slope_exp)))
            # Power-law initial guess: log|y| ~ -p*log(1+x)
            lx = np.log1p(xv)
            if lx.max() > lx.min() + 1e-8:
                slope_pow, _ = np.polyfit(lx, yv, 1)
                p0 = max(0.01, min(8.0, float(-slope_pow)))
            else:
                p0 = 1.0

        sgn = np.sign(y)
        crossing = np.where(sgn[:-1] * sgn[1:] < 0.0)[0]
        if crossing.size >= 2:
            dx = np.diff(x[crossing])
            dx = dx[dx > 1e-12]
            if dx.size > 0:
                omega0 = float(np.clip(np.pi / np.mean(dx), 1e-4, 40.0))
            else:
                omega0 = 1.0
        else:
            omega0 = 1.0

        try:
            from scipy.optimize import minimize as _sp_minimize

            sigma = np.maximum(np.abs(y), 0.05 * amp_peak)
            best = None

            def _maybe_take(candidate):
                nonlocal best
                if candidate is None:
                    return
                if (best is None) or (float(candidate.get("rel_rmse", 1e99)) < float(best.get("rel_rmse", 1e99))):
                    best = candidate

            def _decay_envelope(param, decay):
                """Compute decay envelope vector for the given model."""
                if decay == "exp":
                    return np.exp(-param * x)
                else:  # "power": (1+x)^(-param)
                    return (1.0 + x) ** (-param)

            def _solve_linear(param, omega, n_harmonics, weights=None, decay="exp"):
                """Given nonlinear params (param, omega), solve linear coefficients exactly.
                param is alpha for exp decay, power exponent p for power-law decay.
                Returns ([c_cos, c_sin, c2_cos, c2_sin], rel_rmse) or None on failure.
                """
                e = _decay_envelope(param, decay)
                cols = [e * np.cos(omega * x), e * np.sin(omega * x)]
                if n_harmonics >= 2:
                    cols += [e * np.cos(2.0 * omega * x), e * np.sin(2.0 * omega * x)]
                D = np.column_stack(cols)
                try:
                    if weights is not None:
                        solved, *_ = np.linalg.lstsq(D * weights[:, None], y * weights, rcond=None)
                    else:
                        solved, *_ = np.linalg.lstsq(D, y, rcond=None)
                    yhat = D @ solved
                    rel_rmse = float(np.sqrt(np.mean((yhat - y) ** 2)) / max(amp_peak, 1e-14))
                    coeff = [float(c) for c in solved]
                    # Pad to 4 elements.
                    while len(coeff) < 4:
                        coeff.append(0.0)
                    return coeff, rel_rmse
                except Exception:
                    return None

            def _varpro_obj(ao, n_harmonics, weights, decay):
                """Varpro objective: RMSE after solving linear params exactly."""
                param, omega = float(ao[0]), float(ao[1])
                if param < 0.0 or omega < 1e-6:
                    return 1e20
                result = _solve_linear(param, omega, n_harmonics, weights, decay)
                if result is None:
                    return 1e20
                return result[1]

            # Multiple (param, omega) starting points to escape local minima.
            omega_starts = sorted({
                omega0,
                float(np.clip(omega0 * 0.5, 1e-4, 79.0)),
                float(np.clip(omega0 * 2.0, 1e-4, 79.0)),
            })

            # Candidates: (n_harmonics, use_weights, decay_model, param0)
            # decay="exp": envelope = exp(-alpha*q),  param = alpha in [1e-8, 5]
            # decay="power": envelope = (1+q)^(-p), param = p in [0.01, 8]
            candidates = [
                (1, False, "exp",   alpha0),
                (1, True,  "exp",   alpha0),
                (2, False, "exp",   alpha0),
                (2, True,  "exp",   alpha0),
                (1, False, "power", p0),
                (1, True,  "power", p0),
                (2, False, "power", p0),
                (2, True,  "power", p0),
            ]
            bounds_by_decay = {
                "exp":   [(1e-8, 5.0), (1e-4, 80.0)],
                "power": [(0.01, 8.0), (1e-4, 80.0)],
            }
            for n_harm, use_w, decay, param0 in candidates:
                w = (1.0 / sigma) if use_w else None
                bounds = bounds_by_decay[decay]
                best_local_rmse = 1e20
                best_local_ao = None
                for omega_start in omega_starts:
                    try:
                        res = _sp_minimize(
                            _varpro_obj,
                            x0=[param0, omega_start],
                            args=(n_harm, w, decay),
                            method="L-BFGS-B",
                            bounds=bounds,
                            options={"maxiter": 2000, "ftol": 1e-14, "gtol": 1e-10},
                        )
                        if res.fun < best_local_rmse:
                            best_local_rmse = float(res.fun)
                            best_local_ao = res.x
                    except Exception:
                        pass
                if best_local_ao is None:
                    continue
                param_fit = float(best_local_ao[0])
                omega_fit = float(best_local_ao[1])
                lin = _solve_linear(param_fit, omega_fit, n_harm, w, decay)
                if lin is None:
                    continue
                coeff, rel_rmse = lin
                entry = {
                    "component": self._dominant_component,
                    "q_start": float(q_start),
                    "alpha": param_fit,
                    "omega": omega_fit,
                    "c_cos": coeff[0],
                    "c_sin": coeff[1],
                    "c2_cos": coeff[2] if n_harm >= 2 else 0.0,
                    "c2_sin": coeff[3] if n_harm >= 2 else 0.0,
                    "harmonics": n_harm,
                    "decay_model": decay,
                    "rel_rmse": rel_rmse,
                    "fit_num_cycles": int(fit_num_cycles),
                }
                _maybe_take(entry)

            self._tail_fit = best
        except Exception:
            self._tail_fit = None

    def _eval_tail_fit(self, q):
        if not isinstance(self._tail_fit, dict):
            return None
        qf = float(q)
        if qf < float(self._tail_fit.get("q_start", 0.0)):
            return None
        alpha = float(self._tail_fit.get("alpha", 0.0))
        omega = float(self._tail_fit.get("omega", 0.0))
        c_cos = float(self._tail_fit.get("c_cos", 0.0))
        c_sin = float(self._tail_fit.get("c_sin", 0.0))
        c2_cos = float(self._tail_fit.get("c2_cos", 0.0))
        c2_sin = float(self._tail_fit.get("c2_sin", 0.0))
        harmonics = int(self._tail_fit.get("harmonics", 1))
        decay_model = self._tail_fit.get("decay_model", "exp")

        if decay_model == "power":
            decay_factor = (1.0 + qf) ** (-alpha)
        else:
            decay_factor = np.exp(-alpha * qf)

        val = decay_factor * (c_cos * np.cos(omega * qf) + c_sin * np.sin(omega * qf))
        if harmonics >= 2:
            val += decay_factor * (c2_cos * np.cos(2.0 * omega * qf) + c2_sin * np.sin(2.0 * omega * qf))
        valf = float(val)

        if self._tail_fit.get("component", "both") == "real":
            return complex(valf, 0.0)
        if self._tail_fit.get("component", "both") == "imag":
            return complex(0.0, valf)
        return None

    def _encode_value(self, value, q):
        if self.use_q_scaling:
            return float(value) * (1.0 + float(q))
        return float(value)

    def _decode_value(self, encoded, q):
        if self.use_q_scaling:
            return float(encoded) / (1.0 + float(q))
        return float(encoded)

    def _interp_complex_from(self, interp_real, interp_imag, q):
        er = float(interp_real(q)) if interp_real is not None else 0.0
        ei = float(interp_imag(q)) if interp_imag is not None else 0.0
        return complex(self._decode_value(er, q), self._decode_value(ei, q))

    def _phase_corrected_predict(self, q0, q1, v0, v1, q):
        if q1 <= q0:
            return complex(v0)
        z0 = complex(v0)
        z1 = complex(v1)
        a0 = abs(z0)
        a1 = abs(z1)
        if a0 < 1e-30 or a1 < 1e-30:
            t = (q - q0) / (q1 - q0)
            return z0 + t * (z1 - z0)

        dphi = math.atan2((z1 * z0.conjugate()).imag, (z1 * z0.conjugate()).real)
        omega = dphi / (q1 - q0)
        u0 = z0 * complex(math.cos(-omega * q0), math.sin(-omega * q0))
        u1 = z1 * complex(math.cos(-omega * q1), math.sin(-omega * q1))
        t = (q - q0) / (q1 - q0)
        uq = u0 + t * (u1 - u0)
        return uq * complex(math.cos(omega * q), math.sin(omega * q))

    def _interpolate_phase_corrected(self, qf):
        xs = self._distance_cache_float
        if len(xs) < 2:
            return None
        i = bisect.bisect_right(xs, qf) - 1
        if i < 0:
            i = 0
        if i >= len(xs) - 1:
            i = len(xs) - 2
        q0 = xs[i]
        q1 = xs[i + 1]
        v0 = complex(float(self.value_cache_real_mpf[i]), float(self.value_cache_imag_mpf[i]))
        v1 = complex(float(self.value_cache_real_mpf[i + 1]), float(self.value_cache_imag_mpf[i + 1]))
        return self._phase_corrected_predict(q0, q1, v0, v1, qf)

    def _build_low_q_patch(self):
        q_hi = float(self.low_q_switch_effective)
        self._interp_low_q_real = None
        self._interp_low_q_imag = None
        if q_hi <= 0.0:
            return

        xs = []
        rs = []
        ims = []
        for q_mpf, r_mpf, i_mpf in zip(self.distance_cache_mpf, self.value_cache_real_mpf, self.value_cache_imag_mpf):
            qf = float(q_mpf)
            if qf <= q_hi * 1.0000000001:
                xs.append(qf)
                rs.append(self._encode_value(float(r_mpf), qf))
                ims.append(self._encode_value(float(i_mpf), qf))

        if len(xs) < 4:
            return

        if self._dominant_component != "imag":
            self._interp_low_q_real = self._PchipInterpolator(np.array(xs), np.array(rs), extrapolate=False)
        else:
            self._interp_low_q_real = None
        if self._dominant_component != "real":
            self._interp_low_q_imag = self._PchipInterpolator(np.array(xs), np.array(ims), extrapolate=False)
        else:
            self._interp_low_q_imag = None

    def _ensure_zero_tail_anchor(self):
        """Force terminal cache point at 10*qr_max_effective with zero amplitude.

        Skipped when a tail fit is active: the fit handles extrapolation beyond
        qr_max_effective, so a forced-zero anchor would only distort the PCHIP.
        """
        if not self.distance_cache_mpf:
            return

        # When a fit is present the fit extrapolates beyond qr_max_effective;
        # adding a forced-zero anchor would pull the PCHIP down prematurely.
        if isinstance(self._tail_fit, dict):
            return

        base_qr_max = max(float(self.qr_max_effective), 0.0)
        if base_qr_max <= 0.0:
            return

        target_q = 10.0 * base_qr_max
        last_q = float(self.distance_cache_mpf[-1])
        tol = max(1e-12, 1e-12 * target_q)

        if abs(last_q - target_q) <= tol:
            self.value_cache_real_mpf[-1] = mpmath.mpf(0.0)
            self.value_cache_imag_mpf[-1] = mpmath.mpf(0.0)
            return

        if last_q < target_q:
            self.distance_cache_mpf.append(mpmath.mpf(target_q))
            self.value_cache_real_mpf.append(mpmath.mpf(0.0))
            self.value_cache_imag_mpf.append(mpmath.mpf(0.0))

    def interpolate(self, qr):
        if self._interp_real is None and self._interp_imag is None:
            raise ValueError("Adaptive interpolator not initialized.")
        qf = float(qr)
        x_min = float(self.distance_cache_mpf[0])
        x_max = float(self.qr_max_effective)
        if qf <= x_min:
            real_part = float(self.value_cache_real_mpf[0])
            imag_part = float(self.value_cache_imag_mpf[0])
        elif qf >= x_max:
            c = self._eval_tail_fit(qf)
            if c is not None:
                real_part = c.real
                imag_part = c.imag
            elif isinstance(self._tail_fit, dict):
                # Fit present but q below q_start -- clamp to boundary fit value.
                c2 = self._eval_tail_fit(x_max)
                if c2 is not None:
                    real_part = c2.real
                    imag_part = c2.imag
                else:
                    real_part = float(self.value_cache_real_mpf[-1])
                    imag_part = float(self.value_cache_imag_mpf[-1])
            else:
                real_part = float(self.value_cache_real_mpf[-1])
                imag_part = float(self.value_cache_imag_mpf[-1])
        elif self._interp_low_q_real is not None or self._interp_low_q_imag is not None:
            if qf <= float(self.low_q_switch_effective):
                c = self._interp_complex_from(self._interp_low_q_real, self._interp_low_q_imag, qf)
                real_part = c.real
                imag_part = c.imag
            else:
                c = self._eval_tail_fit(qf)
                if c is None:
                    c = self._interpolate_phase_corrected(qf)
                if c is None:
                    c = self._interp_complex_from(self._interp_real, self._interp_imag, qf)
                real_part = c.real
                imag_part = c.imag
        else:
            c = self._eval_tail_fit(qf)
            if c is None:
                c = self._interpolate_phase_corrected(qf)
            if c is None:
                c = self._interp_complex_from(self._interp_real, self._interp_imag, qf)
            real_part = c.real
            imag_part = c.imag
        if self._dominant_component == "real":
            imag_part = 0.0
        elif self._dominant_component == "imag":
            real_part = 0.0
        return mpmath.mpc(real_part, imag_part)

    def load_cache(self, precomputed_cache):
        self.distance_cache_mpf = [mpmath.mpf(v) for v in precomputed_cache["distance_cache"]]
        self.value_cache_real_mpf = [mpmath.mpf(v) for v in precomputed_cache["value_cache_real"]]
        self.value_cache_imag_mpf = [mpmath.mpf(v) for v in precomputed_cache["value_cache_imag"]]
        self.qr_max_effective = float(precomputed_cache.get("qr_max_effective", float(self.distance_cache_mpf[-1])))
        self.low_q_switch_effective = float(precomputed_cache.get("low_q_switch_effective", 0.0))
        self.use_q_scaling = bool(precomputed_cache.get("use_q_scaling", True))
        self._dominant_component = precomputed_cache.get("dominant_component", "both")
        self._tail_fit = precomputed_cache.get("tail_fit", None)
        if self._tail_fit is not None:
            # In fit mode we intentionally avoid a separate low-q interpolator split.
            self.low_q_switch_effective = 0.0
        self._ensure_zero_tail_anchor()
        self._build_pchip()
        self._build_low_q_patch()

    def make_cache(
        self,
        qr_min=0.0,
        qr_max=4500.0,
        rel_tol=1e-3,
        initial_points=24,
        max_points=8000,
        auto_qr_max=True,
        auto_threshold=1e-8,
        abs_tol=1e-10,
        low_q_switch=0.1,
        low_q_anchor_points=32,
        auto_low_q_switch=True,
        low_q_change_threshold=1e-4,
        phase_tol=0.35,
        use_q_scaling=True,
        phase_amp_rel_floor=1e-3,
        min_interval_rel=1e-5,
        min_interval_abs=1e-6,
        max_new_points_per_iter=128,
        zero_component_threshold=1e-2,
        use_tail_fit=False,
        tail_fit_start_q=25.0,
        fit_num_cycles=0,
        fit_scan_qr_max=120.0,
        fit_scan_points=240,
        use_point_cache=True,
        point_cache_dir="integral_point_cache",
    ):
        self.use_q_scaling = bool(use_q_scaling)
        point_store = _ChannelIntegralPointStore(
            cache_dir=point_cache_dir,
            n0=self.n0,
            l=self.l,
            m=self.m,
            typ=self.typ,
            part=self.part,
            enabled=bool(use_point_cache),
        )
        eval_q = point_store.evaluate

        def _merge_known_points_in_range(vmap, qmin, qmax):
            known = point_store.get_points_in_range(qmin, qmax)
            if not known:
                return vmap
            merged = dict(vmap)
            for qk, vk in known.items():
                if qk not in merged:
                    merged[qk] = vk
            return merged

        if bool(use_tail_fit) and int(fit_num_cycles) > 0:
            effective_qr_max, value_map = _probe_first_cycles_range(
                self.l,
                self.m,
                self.typ,
                self.part,
                qr_max=qr_max,
                n_probe=max(int(initial_points), 80),
                fit_num_cycles=int(fit_num_cycles),
                scan_qr_max=float(fit_scan_qr_max),
                scan_points=int(fit_scan_points),
                eval_fn=eval_q,
            )
            value_map = _merge_known_points_in_range(value_map, qr_min, effective_qr_max)
            qr_values = sorted(q for q in value_map if q >= qr_min)
            if len(qr_values) < initial_points:
                for q in np.linspace(qr_min, effective_qr_max, num=initial_points):
                    qf = float(q)
                    if qf not in value_map:
                        value_map[qf] = eval_q(qf)
                value_map = _merge_known_points_in_range(value_map, qr_min, effective_qr_max)
                qr_values = sorted(q for q in value_map if q >= qr_min)
        elif auto_qr_max:
            # Probe log-spaced to find where amplitude drops below threshold*peak.
            # Reuse these evaluations as the initial grid seed — they are already
            # log-spaced which gives better coverage near qr=0 than uniform points.
            effective_qr_max, value_map = _probe_significant_range(
                self.l, self.m, self.typ, self.part,
                qr_max=qr_max, threshold=auto_threshold,
                eval_fn=eval_q,
            )
            value_map = _merge_known_points_in_range(value_map, qr_min, effective_qr_max)
            qr_values = sorted(q for q in value_map if q >= qr_min)
            # Ensure a minimum coverage count with uniform support points in the
            # effective range, so all cache points are spent where the signal matters.
            if len(qr_values) < initial_points:
                for q in np.linspace(qr_min, effective_qr_max, num=initial_points):
                    qf = float(q)
                    if qf not in value_map:
                        value_map[qf] = eval_q(qf)
                value_map = _merge_known_points_in_range(value_map, qr_min, effective_qr_max)
                qr_values = sorted(q for q in value_map if q >= qr_min)
        else:
            effective_qr_max = qr_max
            qr_values = [float(v) for v in np.linspace(qr_min, qr_max, num=initial_points)]
            value_map = {q: eval_q(q) for q in qr_values}
            value_map = _merge_known_points_in_range(value_map, qr_min, effective_qr_max)
            qr_values = sorted(q for q in value_map if q >= qr_min)

        if auto_low_q_switch and qr_min <= 0.0:
            detected = _auto_detect_low_q_switch(
                value_map,
                effective_qr_max,
                fallback_switch=low_q_switch,
                change_threshold=low_q_change_threshold,
            )
            self.low_q_switch_effective = max(0.0, min(float(detected), float(effective_qr_max)))
        else:
            self.low_q_switch_effective = max(0.0, min(float(low_q_switch), float(effective_qr_max)))
        if qr_min <= 0.0 and self.low_q_switch_effective > 0.0:
            q_hi = self.low_q_switch_effective
            n_low = max(int(low_q_anchor_points), 8)
            q_start = min(1e-4, q_hi)
            if q_hi > 0.0:
                if q_start > 0.0:
                    low_pts = list(np.geomspace(q_start, q_hi, num=n_low))
                else:
                    low_pts = [0.0]
                forced_pts = [0.0, q_hi]
                if q_hi >= 1e-4:
                    forced_pts.append(1e-4)
                for q in sorted(set([float(v) for v in low_pts + forced_pts])):
                    if q >= qr_min and q not in value_map:
                        value_map[q] = eval_q(q)
            qr_values = sorted(q for q in value_map if q >= qr_min)

        self.qr_max_effective = effective_qr_max
        min_interval = max(float(min_interval_abs), float(min_interval_rel) * max(float(effective_qr_max), 1.0))
        eps = 1e-20
        peak_amp = max((abs(complex(v)) for v in value_map.values()), default=0.0)
        # Floor the denominator to avoid over-refining tiny-amplitude tail oscillations.
        scale_floor = max(peak_amp * 1e-12, 1e-14)
        phase_amp_floor = max(float(phase_amp_rel_floor) * max(peak_amp, 1e-30), 1e-14)
        iterations = 0
        midpoint_evals = 0
        last_iter_max_abs_err = 0.0
        last_iter_max_rel_err = 0.0
        last_iter_worst_q = None
        last_iter_max_phase_err = 0.0
        tolerance_met = False

        while len(qr_values) < max_points:
            iterations += 1
            qr_values = sorted(qr_values)
            values = [value_map[q] for q in qr_values]
            self.distance_cache_mpf = [mpmath.mpf(q) for q in qr_values]
            self.value_cache_real_mpf = [mpmath.mpf(complex(v).real) for v in values]
            self.value_cache_imag_mpf = [mpmath.mpf(complex(v).imag) for v in values]
            self._build_pchip()

            to_add = []
            iter_max_abs_err = 0.0
            iter_max_rel_err = 0.0
            iter_worst_q = None
            iter_max_phase_err = 0.0
            for i in range(len(qr_values) - 1):
                q0 = qr_values[i]
                q1 = qr_values[i + 1]
                if (q1 - q0) <= min_interval:
                    continue
                v0 = value_map[q0]
                v1 = value_map[q1]
                qm = (q0 + q1) / 2.0
                if qm in value_map:
                    continue
                vm = eval_q(qm)
                value_map[qm] = vm
                midpoint_evals += 1
                amp = abs(complex(vm))
                if amp > peak_amp:
                    peak_amp = amp
                    scale_floor = max(peak_amp * 1e-12, 1e-14)
                    phase_amp_floor = max(float(phase_amp_rel_floor) * max(peak_amp, 1e-30), 1e-14)

                vpred = self._phase_corrected_predict(q0, q1, v0, v1, qm)
                abs_err = abs(complex(vm) - vpred)
                scale = max(abs(complex(v0)), abs(complex(v1)), amp, scale_floor, eps)
                rel_err = abs_err / scale
                vp_amp = abs(vpred)
                if min(amp, vp_amp) >= phase_amp_floor:
                    phase_err = abs(
                        math.atan2(
                            (complex(vm) * complex(vpred).conjugate()).imag,
                            (complex(vm) * complex(vpred).conjugate()).real,
                        )
                    )
                else:
                    phase_err = 0.0
                if rel_err > iter_max_rel_err:
                    iter_max_rel_err = rel_err
                    iter_max_abs_err = abs_err
                    iter_worst_q = qm
                if phase_err > iter_max_phase_err:
                    iter_max_phase_err = phase_err
                if (abs_err > (abs_tol + rel_tol * scale)) or (phase_err > float(phase_tol)):
                    to_add.append((rel_err, qm))

            last_iter_max_abs_err = iter_max_abs_err
            last_iter_max_rel_err = iter_max_rel_err
            last_iter_worst_q = iter_worst_q
            last_iter_max_phase_err = iter_max_phase_err
            if not to_add:
                tolerance_met = True
                break
            capacity = max_points - len(qr_values)
            add_limit = capacity
            if max_new_points_per_iter is not None and int(max_new_points_per_iter) > 0:
                add_limit = min(add_limit, int(max_new_points_per_iter))
            if len(to_add) > add_limit:
                to_add = sorted(to_add, key=lambda x: x[0], reverse=True)[:add_limit]
            new_points = [qm for _, qm in to_add]
            qr_values = sorted(set(qr_values + new_points))

        qr_values = sorted(qr_values)
        values = [value_map[q] for q in qr_values]
        self.distance_cache_mpf = [mpmath.mpf(q) for q in qr_values]
        self.value_cache_real_mpf = [mpmath.mpf(complex(v).real) for v in values]
        self.value_cache_imag_mpf = [mpmath.mpf(complex(v).imag) for v in values]
        # Detect which component is physically significant (the other is a numerical artifact).
        probe_max = min(10.0, float(self.qr_max_effective))
        self._dominant_component = self._detect_dominant_component(
            value_map, probe_max=probe_max, threshold=float(zero_component_threshold)
        )
        # Zero out the non-dominant component in the cache so it is stored as exact 0.
        if self._dominant_component == "real":
            self.value_cache_imag_mpf = [mpmath.mpf(0.0)] * len(self.value_cache_imag_mpf)
        elif self._dominant_component == "imag":
            self.value_cache_real_mpf = [mpmath.mpf(0.0)] * len(self.value_cache_real_mpf)

        self._tail_fit = None
        fit_start_effective = float(tail_fit_start_q)
        if int(fit_num_cycles) > 0 and fit_start_effective >= 0.8 * float(self.qr_max_effective):
            fit_start_effective = max(0.0, float(qr_min))
        if bool(use_tail_fit):
            self._fit_dominant_tail(tail_fit_start_q=fit_start_effective, fit_num_cycles=fit_num_cycles)
            # In fit mode we intentionally avoid a separate low-q interpolator split.
            self.low_q_switch_effective = 0.0
        self._ensure_zero_tail_anchor()
        self._build_pchip()
        self._build_low_q_patch()
        self.build_stats = {
            "iterations": int(iterations),
            "midpoint_evaluations": int(midpoint_evals),
            "grid_points": int(len(qr_values)),
            "tolerance_met": bool(tolerance_met),
            "hit_point_cap": bool((len(qr_values) >= int(max_points)) and (not tolerance_met)),
            "last_iter_max_abs_err": float(last_iter_max_abs_err),
            "last_iter_max_rel_err": float(last_iter_max_rel_err),
            "last_iter_worst_q": (None if last_iter_worst_q is None else float(last_iter_worst_q)),
            "last_iter_max_phase_err": float(last_iter_max_phase_err),
            "low_q_switch_effective": float(self.low_q_switch_effective),
            "qr_max_effective": float(self.qr_max_effective),
            "phase_tol": float(phase_tol),
            "phase_amp_rel_floor": float(phase_amp_rel_floor),
            "min_interval": float(min_interval),
            "max_new_points_per_iter": int(max_new_points_per_iter) if max_new_points_per_iter is not None else None,
            "use_q_scaling": bool(self.use_q_scaling),
            "dominant_component": self._dominant_component,
            "use_tail_fit": bool(self._tail_fit is not None),
            "tail_fit_start_q": float(tail_fit_start_q),
            "tail_fit_start_q_effective": float(fit_start_effective),
            "fit_num_cycles": int(fit_num_cycles),
            "use_point_cache": bool(use_point_cache),
            "point_cache_dir": str(point_cache_dir),
        }
        point_store.flush()


# ---------------------------------------------------------------------------
# Top-level worker functions for parallel generation (must be picklable)
# ---------------------------------------------------------------------------

def _probe_significant_range(
    l,
    m,
    typ,
    part,
    qr_max=4500.0,
    threshold=1e-10,
    n_probe=120,
    search_growth=2.0,
    search_max_multiplier=16.0,
    tail_fraction=0.2,
    tail_confirm_rounds=2,
    eval_fn=None,
):
    """Find effective qr range with progressive, channel-specific probing.

    The probe range is expanded until the high-qr tail is consistently below
    threshold * peak amplitude. This avoids relying on one fixed stop value for
    all channels while still keeping a safety cap.

    Returns (effective_qr_max, value_map) where:
    - effective_qr_max: first envelope crossing q where |I|_env <= threshold*peak,
      padded by a few probe points.
    - value_map: {qr: I_integral(qr,...)} for probe points up to effective_qr_max.
    """
    eval_q = eval_fn if eval_fn is not None else (lambda q: I_integral((q, l, m, typ, part)))

    base_qmax = max(float(qr_max), 1.0)
    current_qmax = max(base_qmax, 25.0)
    cap_qmax = max(current_qmax + 50.0, current_qmax * float(search_max_multiplier))

    n_probe = max(int(n_probe), 96)
    growth = max(float(search_growth), 1.25)
    tail_fraction = min(max(float(tail_fraction), 0.05), 0.5)
    rounds_needed = max(int(tail_confirm_rounds), 1)

    value_map = {}
    consecutive_tail_ok = 0
    last_cross_idx = None
    last_probe_pts = [0.0, current_qmax]
    max_rounds = 24

    for _ in range(max_rounds):
        log_start = max(current_qmax * 1e-5, 1e-3)
        log_pts = list(np.geomspace(log_start, current_qmax, num=n_probe))
        lin_pts = list(np.linspace(0.0, current_qmax, num=n_probe))
        probe_pts = sorted(set([0.0] + [float(v) for v in log_pts + lin_pts]))
        last_probe_pts = probe_pts

        for q in probe_pts:
            if q not in value_map:
                value_map[q] = eval_q(q)

        amps = [abs(complex(value_map[q])) for q in probe_pts]
        peak = max(amps) if amps else 0.0
        if peak <= 0.0:
            return current_qmax, {q: value_map[q] for q in probe_pts}

        peak_idx = int(np.argmax(np.array(amps)))
        cutoff = float(threshold) * peak

        env = [0.0] * len(amps)
        run = 0.0
        for i in range(len(amps) - 1, -1, -1):
            run = max(run, amps[i])
            env[i] = run

        cross_idx = None
        for i, e in enumerate(env):
            if e <= cutoff:
                cross_idx = i
                break

        tail_start = int((1.0 - tail_fraction) * (len(amps) - 1))
        tail_start = min(max(tail_start, peak_idx + 1), len(amps) - 1)
        tail_max = max(amps[tail_start:]) if tail_start < len(amps) else amps[-1]
        tail_ok = (tail_max <= cutoff) and (cross_idx is not None)

        if tail_ok:
            consecutive_tail_ok += 1
            last_cross_idx = cross_idx
        else:
            consecutive_tail_ok = 0

        reached_cap = current_qmax >= cap_qmax * (1.0 - 1e-12)
        if consecutive_tail_ok >= rounds_needed or reached_cap:
            if last_cross_idx is None:
                if cross_idx is not None:
                    last_cross_idx = cross_idx
                else:
                    last_cross_idx = len(probe_pts) - 1
            padded_idx = min(last_cross_idx + 3, len(probe_pts) - 1)
            effective_max = float(probe_pts[padded_idx])
            filtered = {q: v for q, v in value_map.items() if q <= effective_max * 1.001}
            return effective_max, filtered

        current_qmax = min(cap_qmax, current_qmax * growth)

    # Fallback (should be rare): return the explored range.
    effective_max = float(last_probe_pts[-1])
    filtered = {q: v for q, v in value_map.items() if q <= effective_max * 1.001}
    return effective_max, filtered


def _auto_detect_low_q_switch(value_map, effective_qr_max, fallback_switch=0.1, change_threshold=1e-4):
    """Detect a low-q patch handoff from probe values.

    The switch is the first q where |I(q)-I(0)| exceeds change_threshold times
    the peak amplitude over the effective range.
    """
    q_eff = max(float(effective_qr_max), 0.0)
    fallback = max(0.0, min(float(fallback_switch), q_eff))
    if q_eff <= 0.0 or 0.0 not in value_map:
        return fallback

    probe_q = sorted(float(q) for q in value_map if q >= 0.0 and q <= q_eff * 1.0000000001)
    if len(probe_q) < 4:
        return fallback

    base = complex(value_map[0.0])
    peak = max((abs(complex(value_map[q])) for q in probe_q), default=0.0)
    if peak <= 0.0:
        return fallback

    cutoff = max(float(change_threshold), 1e-12) * peak
    for q in probe_q[1:]:
        delta = abs(complex(value_map[q]) - base)
        if delta >= cutoff:
            return max(min(float(q), q_eff), 1e-4)

    return fallback


def _tail_envelope_stats(distance_cache, value_cache_real, value_cache_imag, thresholds=(1e-6, 1e-8, 1e-10, 1e-12)):
    """Compute robust tail cutoff recommendations from right-envelope amplitudes.

    The envelope at index i is max(|I(q_j)| for j>=i), which avoids false cutoffs
    caused by oscillatory zero crossings.
    """
    xs = [float(v) for v in distance_cache]
    vals = [complex(float(r), float(im)) for r, im in zip(value_cache_real, value_cache_imag)]
    amps = [abs(v) for v in vals]
    peak = max(amps) if amps else 0.0
    if not amps:
        return {
            "peak_abs": 0.0,
            "peak_qr": 0.0,
            "cutoffs": {str(t): 0.0 for t in thresholds},
        }

    env = [0.0] * len(amps)
    run = 0.0
    for i in range(len(amps) - 1, -1, -1):
        run = max(run, amps[i])
        env[i] = run

    peak_idx = int(np.argmax(np.array(amps)))
    cutoffs = {}
    if peak <= 0.0:
        for t in thresholds:
            cutoffs[str(t)] = float(xs[-1])
    else:
        for t in thresholds:
            target = float(t) * peak
            q_cut = float(xs[-1])
            for i, e in enumerate(env):
                if e <= target:
                    q_cut = float(xs[i])
                    break
            cutoffs[str(t)] = q_cut

    return {
        "peak_abs": float(peak),
        "peak_qr": float(xs[peak_idx]),
        "cutoffs": cutoffs,
    }


def _worker_linear(args):
    """Build one uniform-grid interpolator entry.  Returns (key, data_dict)."""
    n0, l, m, typ, part, qr_min, qr_max, num_qr = args
    qr_values = np.linspace(qr_min, qr_max, num=num_qr)
    interp = I_interpolator(n0, l, m, typ, part, qr_values=qr_values)
    key = _interpolator_key(n0, l, m, typ, part)
    data = {
        "distance_cache": [float(v) for v in interp.distance_cache_mpf],
        "value_cache_real": [float(v) for v in interp.value_cache_real_mpf],
        "value_cache_imag": [float(v) for v in interp.value_cache_imag_mpf],
    }
    return key, data


def _worker_adaptive(args):
    """Build one adaptive PCHIP interpolator entry.  Returns (key, data_dict)."""
    n0, l, m, typ, part, qr_min, qr_max, rel_tol, initial_points, max_points, auto_qr_max, auto_threshold, abs_tol, low_q_switch, low_q_anchor_points, auto_low_q_switch, low_q_change_threshold, phase_tol, use_q_scaling, phase_amp_rel_floor, min_interval_rel, min_interval_abs, max_new_points_per_iter, zero_component_threshold, use_tail_fit, tail_fit_start_q, fit_num_cycles, fit_scan_qr_max, fit_scan_points, use_point_cache, point_cache_dir = args
    interp = I_adaptive_interpolator(
        n0, l, m, typ, part,
        qr_min=qr_min, qr_max=qr_max,
        rel_tol=rel_tol, initial_points=initial_points, max_points=max_points,
        auto_qr_max=auto_qr_max, auto_threshold=auto_threshold, abs_tol=abs_tol,
        low_q_switch=low_q_switch, low_q_anchor_points=low_q_anchor_points,
        auto_low_q_switch=auto_low_q_switch, low_q_change_threshold=low_q_change_threshold,
        phase_tol=phase_tol, use_q_scaling=use_q_scaling,
        phase_amp_rel_floor=phase_amp_rel_floor,
        min_interval_rel=min_interval_rel, min_interval_abs=min_interval_abs,
        max_new_points_per_iter=max_new_points_per_iter,
        zero_component_threshold=zero_component_threshold,
        use_tail_fit=use_tail_fit,
        tail_fit_start_q=tail_fit_start_q,
        fit_num_cycles=fit_num_cycles,
        fit_scan_qr_max=fit_scan_qr_max,
        fit_scan_points=fit_scan_points,
        use_point_cache=use_point_cache,
        point_cache_dir=point_cache_dir,
    )
    key = _interpolator_key(n0, l, m, typ, part)
    data = {
        "distance_cache": [float(v) for v in interp.distance_cache_mpf],
        "value_cache_real": [float(v) for v in interp.value_cache_real_mpf],
        "value_cache_imag": [float(v) for v in interp.value_cache_imag_mpf],
        "qr_max_effective": interp.qr_max_effective,
        "low_q_switch_effective": interp.low_q_switch_effective,
        "use_q_scaling": bool(interp.use_q_scaling),
        "dominant_component": interp._dominant_component,
        "tail_fit": interp._tail_fit,
        "build_stats": dict(interp.build_stats),
    }
    data["tail_stats"] = _tail_envelope_stats(
        data["distance_cache"],
        data["value_cache_real"],
        data["value_cache_imag"],
    )
    return key, data


def _worker_cutoff_probe(args):
    n0, l, m, typ, part, qr_min, qr_max, n_probe, thresholds, integration_dps = args
    if integration_dps is not None:
        mpmath.mp.dps = int(integration_dps)
    log_start = max(qr_max * 1e-5, 1e-3)
    n_probe = max(int(n_probe), 80)
    log_pts = list(np.geomspace(log_start, qr_max, num=n_probe))
    lin_pts = list(np.linspace(qr_min, qr_max, num=n_probe))
    qs = sorted(set([float(qr_min)] + [float(q) for q in log_pts + lin_pts if q >= qr_min]))

    vals = [complex(I_integral((q, l, m, typ, part))) for q in qs]
    amps = [abs(v) for v in vals]
    peak = max(amps) if amps else 0.0
    peak_idx = int(np.argmax(np.array(amps))) if amps else 0

    env = [0.0] * len(amps)
    run = 0.0
    for i in range(len(amps) - 1, -1, -1):
        run = max(run, amps[i])
        env[i] = run

    cutoffs_strict = {}
    cutoffs_relaxed_5pct = {}
    cutoffs_relaxed_20pct = {}
    for thr in thresholds:
        if peak <= 0.0:
            key_thr = str(thr)
            cutoffs_strict[key_thr] = float(qr_min)
            cutoffs_relaxed_5pct[key_thr] = float(qr_min)
            cutoffs_relaxed_20pct[key_thr] = float(qr_min)
            continue
        target = float(thr) * peak

        # Strict: require all later points to stay below target.
        q_cut_strict = float(qr_max)
        for i, e in enumerate(env):
            if e <= target:
                q_cut_strict = float(qs[i])
                break

        # Robust: allow occasional high-qr spikes; use tail exceedance fraction.
        above = [1 if a > target else 0 for a in amps]
        tail_count = len(above)
        running_above = 0
        q_cut_5pct = float(qr_max)
        q_cut_20pct = float(qr_max)
        for i in range(len(above) - 1, -1, -1):
            running_above += above[i]
            tail_len = tail_count - i
            frac_above = running_above / float(tail_len)
            if frac_above <= 0.20:
                q_cut_20pct = float(qs[i])
            if frac_above <= 0.05:
                q_cut_5pct = float(qs[i])

        key_thr = str(thr)
        cutoffs_strict[key_thr] = q_cut_strict
        cutoffs_relaxed_5pct[key_thr] = q_cut_5pct
        cutoffs_relaxed_20pct[key_thr] = q_cut_20pct

    key = _interpolator_key(n0, l, m, typ, part)
    data = {
        "peak_abs": float(peak),
        "peak_qr": float(qs[peak_idx]) if qs else float(qr_min),
        "cutoffs": cutoffs_strict,
        "cutoffs_relaxed_5pct": cutoffs_relaxed_5pct,
        "cutoffs_relaxed_20pct": cutoffs_relaxed_20pct,
    }
    return key, data


# ---------------------------------------------------------------------------
# Builder helpers (used by interpolgrator)
# ---------------------------------------------------------------------------

def _build_interpolator(n0, l, m, typ, part, precomputed_data=None, qr_values=None):
    key = _interpolator_key(n0, l, m, typ, part)
    if precomputed_data is not None and key in precomputed_data:
        return I_interpolator(n0, l, m, typ, part, precomputed_cache=precomputed_data[key])
    return I_interpolator(n0, l, m, typ, part, qr_values=qr_values)


def _build_adaptive_interpolator(n0, l, m, typ, part, precomputed_data=None):
    key = _interpolator_key(n0, l, m, typ, part)
    if precomputed_data is not None and key in precomputed_data:
        return I_adaptive_interpolator(n0, l, m, typ, part, precomputed_cache=precomputed_data[key])
    return I_adaptive_interpolator(n0, l, m, typ, part)


# ---------------------------------------------------------------------------
# JSON generate / load / enable  —  linear (uniform) interpolator
# ---------------------------------------------------------------------------

def generate_interpolators_file(
    output_file="interpolators.json",
    n0_values=range(1, 6),
    l_max=10,
    qr_min=0.0,
    qr_max=4500.0,
    num_qr=4500,
    num_workers=None,
):
    """Create a JSON file with precomputed uniform-grid helper values.

    Parameters
    ----------
    num_workers : int or None
        Number of worker processes.  None uses all available CPU cores.
    """
    from multiprocessing import Pool, cpu_count

    tasks = [
        (n0, l, m, typ, part, qr_min, qr_max, num_qr)
        for n0 in n0_values
        for l in range(0, l_max + 1)
        for typ in ("s", "p", "d")
        if _n0_typ_valid(n0, typ)
        for m in get_m_values_for_type_l(typ, l)
        for part in get_parts_for_typ_m(typ, m)
    ]
    total = len(tasks)
    workers = num_workers if num_workers is not None else cpu_count()
    print("Building {} uniform-grid interpolators with {} workers ...".format(total, workers))

    results = {}
    with Pool(processes=workers) as pool:
        for done, (key, data) in enumerate(pool.imap_unordered(_worker_linear, tasks), start=1):
            results[key] = data
            print("  [{}/{}] {}".format(done, total, key))

    serializable = {
        "metadata": {
            "n0_values": list(n0_values),
            "l_values": list(range(0, l_max + 1)),
            "types": ["s", "p", "d"],
            "qr_min": float(qr_min),
            "qr_max": float(qr_max),
            "num_qr": int(num_qr),
        },
        "interpolators": results,
    }
    with open(output_file, "w") as out:
        json.dump(serializable, out)

    print("Saved {} interpolators to {}".format(len(results), output_file))
    return serializable


def load_interpolators_file(input_file="interpolators.json"):
    """Load precomputed uniform-grid helper values from JSON."""
    with open(input_file, "r") as inp:
        raw = json.load(inp)
    return raw.get("interpolators", {})


def enable_precomputed_interpolators(input_file="interpolators.json"):
    """Load and store uniform-grid interpolators in the module global."""
    global PRECOMPUTED_INTERPOLATORS
    PRECOMPUTED_INTERPOLATORS = load_interpolators_file(input_file)
    return len(PRECOMPUTED_INTERPOLATORS)


# ---------------------------------------------------------------------------
# JSON generate / load / enable  —  adaptive PCHIP interpolator
# ---------------------------------------------------------------------------

def generate_adaptive_interpolators_file(
    output_file="interpolator_adaptive.json",
    n0_values=range(1, 6),
    l_max=10,
    qr_min=0.0,
    qr_max=4500.0,
    rel_tol=1e-3,
    initial_points=24,
    max_points=8000,
    num_workers=None,
    auto_qr_max=True,
    auto_threshold=1e-10,
    abs_tol=1e-10,
    low_q_switch=0.1,
    low_q_anchor_points=32,
    auto_low_q_switch=True,
    low_q_change_threshold=1e-4,
    phase_tol=0.35,
    use_q_scaling=True,
    phase_amp_rel_floor=1e-3,
    min_interval_rel=1e-5,
    min_interval_abs=1e-6,
    max_new_points_per_iter=128,
    cutoff_file=DEFAULT_CUTOFF_FILE,
    auto_use_cutoffs=True,
    cutoff_field="cutoffs_relaxed_5pct",
    cutoff_threshold="0.2",
    cutoff_pad_factor=1.02,
    cutoff_min_qr_max=1.0,
    zero_component_threshold=1e-2,
    use_tail_fit=False,
    tail_fit_start_q=25.0,
    fit_num_cycles=0,
    fit_scan_qr_max=120.0,
    fit_scan_points=240,
    use_point_cache=True,
    point_cache_dir="integral_point_cache",
):
    """Create a JSON file with adaptive-grid real/imag helper values (PCHIP).

    Parameters
    ----------
    num_workers : int or None
        Number of worker processes.  None uses all available CPU cores.
    """
    from multiprocessing import Pool, cpu_count

    cutoff_data = _load_cutoff_recommendations_file(cutoff_file) if auto_use_cutoffs else {}
    tasks = []
    for n0 in n0_values:
        for l in range(0, l_max + 1):
            for typ in ("s", "p", "d"):
                if not _n0_typ_valid(n0, typ):
                    continue
                for m in get_m_values_for_type(typ):
                    if m > l:
                        continue
                    for part in get_parts_for_typ_m(typ, m):
                        key = _interpolator_key(n0, l, m, typ, part)
                        channel_qr_max = _resolve_channel_qr_max(
                            key,
                            qr_max,
                            cutoff_data=cutoff_data,
                            cutoff_field=cutoff_field,
                            cutoff_threshold=cutoff_threshold,
                            cutoff_pad_factor=cutoff_pad_factor,
                            min_qr_max=cutoff_min_qr_max,
                        )
                        tasks.append(
                            (
                                n0,
                                l,
                                m,
                                typ,
                                part,
                                qr_min,
                                channel_qr_max,
                                rel_tol,
                                initial_points,
                                max_points,
                                auto_qr_max,
                                auto_threshold,
                                abs_tol,
                                low_q_switch,
                                low_q_anchor_points,
                                auto_low_q_switch,
                                low_q_change_threshold,
                                phase_tol,
                                use_q_scaling,
                                phase_amp_rel_floor,
                                min_interval_rel,
                                min_interval_abs,
                                max_new_points_per_iter,
                                zero_component_threshold,
                                use_tail_fit,
                                tail_fit_start_q,
                                fit_num_cycles,
                                fit_scan_qr_max,
                                fit_scan_points,
                                use_point_cache,
                                point_cache_dir,
                            )
                        )
    total = len(tasks)
    workers = num_workers if num_workers is not None else cpu_count()
    print("Building {} adaptive interpolators with {} workers ...".format(total, workers))
    if auto_use_cutoffs:
        print(
            "Using cutoff file {} ({}, thr={}, pad={})".format(
                cutoff_file,
                cutoff_field,
                cutoff_threshold,
                cutoff_pad_factor,
            )
        )

    results = {}
    point_counts = []
    q_switches = []
    last_iter_rel_errs = []
    last_iter_phase_errs = []
    cutoff_1e10 = []
    cutoff_1e12 = []
    tolerance_met_count = 0
    capped_count = 0
    with Pool(processes=workers) as pool:
        for done, (key, data) in enumerate(pool.imap_unordered(_worker_adaptive, tasks), start=1):
            results[key] = data
            stats = data.get("build_stats", {})
            npts = int(len(data.get("distance_cache", [])))
            qsw = float(data.get("low_q_switch_effective", 0.0))
            rel_last = float(stats.get("last_iter_max_rel_err", 0.0))
            phase_last = float(stats.get("last_iter_max_phase_err", 0.0))
            tail = data.get("tail_stats", {})
            c10 = float(tail.get("cutoffs", {}).get("1e-10", data.get("qr_max_effective", 0.0)))
            c12 = float(tail.get("cutoffs", {}).get("1e-12", data.get("qr_max_effective", 0.0)))
            status = "OK" if stats.get("tolerance_met", False) else "CAPPED"

            point_counts.append(npts)
            q_switches.append(qsw)
            last_iter_rel_errs.append(rel_last)
            last_iter_phase_errs.append(phase_last)
            cutoff_1e10.append(c10)
            cutoff_1e12.append(c12)
            if stats.get("tolerance_met", False):
                tolerance_met_count += 1
            if stats.get("hit_point_cap", False):
                capped_count += 1

            print(
                "  [{}/{}] {} points={} qsw={:.6g} rel_last={:.3e} {}".format(
                    done, total, key, npts, qsw, rel_last, status
                )
            )
            print(
                "           tail q@1e-10={:.6g} q@1e-12={:.6g}".format(c10, c12)
            )

    if point_counts:
        point_arr = np.array(point_counts, dtype=float)
        qsw_arr = np.array(q_switches, dtype=float)
        rel_arr = np.array(last_iter_rel_errs, dtype=float)
        ph_arr = np.array(last_iter_phase_errs, dtype=float)
        c10_arr = np.array(cutoff_1e10, dtype=float)
        c12_arr = np.array(cutoff_1e12, dtype=float)
        print("Adaptive generation stats:")
        print(
            "  gridpoints min/median/max = {}/{:.1f}/{}".format(
                int(np.min(point_arr)), float(np.median(point_arr)), int(np.max(point_arr))
            )
        )
        print(
            "  low-q switch min/max = {:.6g}/{:.6g}".format(
                float(np.min(qsw_arr)), float(np.max(qsw_arr))
            )
        )
        print(
            "  achieved rel_last min/median/max = {:.3e}/{:.3e}/{:.3e}".format(
                float(np.min(rel_arr)), float(np.median(rel_arr)), float(np.max(rel_arr))
            )
        )
        print(
            "  achieved phase_last(rad) min/median/max = {:.3e}/{:.3e}/{:.3e}".format(
                float(np.min(ph_arr)), float(np.median(ph_arr)), float(np.max(ph_arr))
            )
        )
        print(
            "  tail cutoff q@1e-10 min/median/max = {:.6g}/{:.6g}/{:.6g}".format(
                float(np.min(c10_arr)), float(np.median(c10_arr)), float(np.max(c10_arr))
            )
        )
        print(
            "  tail cutoff q@1e-12 min/median/max = {:.6g}/{:.6g}/{:.6g}".format(
                float(np.min(c12_arr)), float(np.median(c12_arr)), float(np.max(c12_arr))
            )
        )
        print(
            "  tolerance_met = {}/{} ; point_cap_hit = {}".format(
                int(tolerance_met_count), int(len(point_counts)), int(capped_count)
            )
        )

    if bool(use_tail_fit):
        print("Tail fit parameters by key:")
        for key in sorted(results.keys()):
            fit = results[key].get("tail_fit", None)
            if isinstance(fit, dict):
                h = int(fit.get("harmonics", 1))
                if h >= 2:
                    print(
                        "  {} component={} q_start={:.6g} alpha={:.6g} omega={:.6g} c1_cos={:.6g} c1_sin={:.6g} c2_cos={:.6g} c2_sin={:.6g} rel_rmse={:.3e} cycles={} harmonics={}".format(
                            key,
                            str(fit.get("component", "both")),
                            float(fit.get("q_start", 0.0)),
                            float(fit.get("alpha", 0.0)),
                            float(fit.get("omega", 0.0)),
                            float(fit.get("c_cos", 0.0)),
                            float(fit.get("c_sin", 0.0)),
                            float(fit.get("c2_cos", 0.0)),
                            float(fit.get("c2_sin", 0.0)),
                            float(fit.get("rel_rmse", 0.0)),
                            int(fit.get("fit_num_cycles", 0)),
                            h,
                        )
                    )
                else:
                    print(
                        "  {} component={} q_start={:.6g} alpha={:.6g} omega={:.6g} c_cos={:.6g} c_sin={:.6g} rel_rmse={:.3e} cycles={} harmonics={}".format(
                            key,
                            str(fit.get("component", "both")),
                            float(fit.get("q_start", 0.0)),
                            float(fit.get("alpha", 0.0)),
                            float(fit.get("omega", 0.0)),
                            float(fit.get("c_cos", 0.0)),
                            float(fit.get("c_sin", 0.0)),
                            float(fit.get("rel_rmse", 0.0)),
                            int(fit.get("fit_num_cycles", 0)),
                            h,
                        )
                    )
            else:
                print("  {} NO_FIT".format(key))

    serializable = {
        "metadata": {
            "n0_values": list(n0_values),
            "l_values": list(range(0, l_max + 1)),
            "types": ["s", "p", "d"],
            "qr_min": float(qr_min),
            "qr_max": float(qr_max),
            "rel_tol": float(rel_tol),
            "initial_points": int(initial_points),
            "max_points": int(max_points),
            "auto_qr_max": bool(auto_qr_max),
            "auto_threshold": float(auto_threshold),
            "abs_tol": float(abs_tol),
            "low_q_switch": float(low_q_switch),
            "low_q_anchor_points": int(low_q_anchor_points),
            "auto_low_q_switch": bool(auto_low_q_switch),
            "low_q_change_threshold": float(low_q_change_threshold),
            "phase_tol": float(phase_tol),
            "phase_amp_rel_floor": float(phase_amp_rel_floor),
            "min_interval_rel": float(min_interval_rel),
            "min_interval_abs": float(min_interval_abs),
            "max_new_points_per_iter": int(max_new_points_per_iter) if max_new_points_per_iter is not None else None,
            "use_q_scaling": bool(use_q_scaling),
            "auto_use_cutoffs": bool(auto_use_cutoffs),
            "cutoff_file": cutoff_file,
            "cutoff_field": cutoff_field,
            "cutoff_threshold": cutoff_threshold,
            "cutoff_pad_factor": float(cutoff_pad_factor),
            "cutoff_min_qr_max": float(cutoff_min_qr_max),
            "use_tail_fit": bool(use_tail_fit),
            "tail_fit_start_q": float(tail_fit_start_q),
            "fit_num_cycles": int(fit_num_cycles),
            "fit_scan_qr_max": float(fit_scan_qr_max),
            "fit_scan_points": int(fit_scan_points),
            "use_point_cache": bool(use_point_cache),
            "point_cache_dir": str(point_cache_dir),
        },
        "interpolators": results,
    }
    with open(output_file, "w") as out:
        json.dump(serializable, out)

    print("Saved {} adaptive interpolators to {}".format(len(results), output_file))
    return serializable


def generate_cutoff_recommendations_file(
    output_file="interpolator_cutoffs.json",
    n0_values=range(1, 6),
    l_max=10,
    qr_min=0.0,
    qr_max=4500.0,
    n_probe=120,
    thresholds=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6),
    num_workers=None,
    integration_dps=35,
):
    """Build direct-probe per-channel qr cutoff recommendations.

    Cutoffs are computed from right-envelope amplitudes and represent the first
    qr where envelope <= threshold * peak.
    """
    from multiprocessing import Pool, cpu_count

    tasks = [
        (n0, l, m, typ, part, qr_min, qr_max, n_probe, tuple(float(t) for t in thresholds), integration_dps)
        for n0 in n0_values
        for l in range(0, l_max + 1)
        for typ in ("s", "p", "d")
        if _n0_typ_valid(n0, typ)
        for m in get_m_values_for_type_l(typ, l)
        for part in get_parts_for_typ_m(typ, m)
    ]

    total = len(tasks)
    workers = num_workers if num_workers is not None else cpu_count()
    print("Analyzing {} channels for cutoff recommendations with {} workers ...".format(total, workers))

    results = {}
    with Pool(processes=workers) as pool:
        for done, (key, data) in enumerate(pool.imap_unordered(_worker_cutoff_probe, tasks), start=1):
            results[key] = data
            t_first = str(float(thresholds[0]))
            t_last = str(float(thresholds[-1]))
            print("  [{}/{}] {} peak={:.3e} q@{}={:.6g} q@{}={:.6g}".format(
                done,
                total,
                key,
                float(data.get("peak_abs", 0.0)),
                t_first,
                float(data.get("cutoffs", {}).get(t_first, qr_max)),
                t_last,
                float(data.get("cutoffs", {}).get(t_last, qr_max)),
            ))
            print("           robust q@{} (5%/20%) = {:.6g}/{:.6g}".format(
                t_first,
                float(data.get("cutoffs_relaxed_5pct", {}).get(t_first, qr_max)),
                float(data.get("cutoffs_relaxed_20pct", {}).get(t_first, qr_max)),
            ))

    serializable = {
        "metadata": {
            "n0_values": list(n0_values),
            "l_values": list(range(0, l_max + 1)),
            "types": ["s", "p", "d"],
            "qr_min": float(qr_min),
            "qr_max": float(qr_max),
            "n_probe": int(n_probe),
            "thresholds": [float(t) for t in thresholds],
            "integration_dps": int(integration_dps) if integration_dps is not None else None,
        },
        "cutoffs": results,
    }
    with open(output_file, "w") as out:
        json.dump(serializable, out)

    print("Saved {} channel cutoffs to {}".format(len(results), output_file))
    return serializable


def load_adaptive_interpolators_file(input_file="interpolator_adaptive.json"):
    """Load precomputed adaptive-grid helper values from JSON."""
    with open(input_file, "r") as inp:
        raw = json.load(inp)
    return raw.get("interpolators", {})


def enable_precomputed_adaptive_interpolators(input_file="interpolator_adaptive.json"):
    """Load and store adaptive-grid interpolators in the module global."""
    global PRECOMPUTED_INTERPOLATORS_ADAPTIVE
    PRECOMPUTED_INTERPOLATORS_ADAPTIVE = load_adaptive_interpolators_file(input_file)
    return len(PRECOMPUTED_INTERPOLATORS_ADAPTIVE)


# ---------------------------------------------------------------------------
# Process-level initialisation (used by multiprocessing Pool)
# ---------------------------------------------------------------------------

def init_precomputed_interpolators(input_file=None, quiet=False):
    """Load interpolator cache for the current process and remember source file."""
    global ACTIVE_INTERPOLATORS_FILE
    resolved_file = resolve_interpolators_file(input_file)
    count = enable_precomputed_interpolators(resolved_file)
    ACTIVE_INTERPOLATORS_FILE = resolved_file
    if not quiet:
        print("Loaded {} interpolators from {}".format(count, resolved_file))
    return resolved_file, count


def _pool_worker_init(input_file=None):
    """Pool initializer: load precomputed interpolators once per worker process."""
    try:
        resolved_file, count = init_precomputed_interpolators(input_file=input_file, quiet=True)
        adaptive_file = "interpolator_adaptive.json"
        adaptive_count = 0
        if os.path.exists(adaptive_file):
            adaptive_count = enable_precomputed_adaptive_interpolators(adaptive_file)
        #print(
        #    "Worker loaded {} linear interpolators from {} and {} adaptive interpolators from {}".format(
        #        count,
        #        resolved_file,
        #        adaptive_count,
        #        adaptive_file,
        #    )
        #)
    except Exception as exc:
        print("Worker failed to load interpolator file: {}".format(exc))


# ---------------------------------------------------------------------------
# Accuracy testers
# ---------------------------------------------------------------------------

def _build_combinations(n0_values, l_max):
    """Return all valid (n0, l, m, typ, part) combinations."""
    combinations = []
    for n0 in n0_values:
        for l in range(0, l_max + 1):
            for typ in ("s", "p", "d"):
                if not _n0_typ_valid(n0, typ):
                    continue
                for m in get_m_values_for_type(typ):
                    if m > l:
                        continue
                    for part in get_parts_for_typ_m(typ, m):
                        combinations.append((n0, l, m, typ, part))
    return combinations


_DEFAULT_TEST_QR = [0.0, 1.93764e-4, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 37.5, 123.4, 777.0]


def _safe_plot_filename(key):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(key))


def _build_adaptive_test_qr_points(qr_min, qr_max, n_points):
    qr_min = float(qr_min)
    qr_max = float(qr_max)
    n_points = max(int(n_points), 24)
    if qr_max <= qr_min:
        return [qr_min]

    q_lo = max(qr_min, min(1e-4, qr_max * 1e-3))
    log_n = max(n_points // 2, 12)
    lin_n = max(n_points // 2, 12)
    log_pts = []
    if q_lo > 0.0:
        log_pts = [float(v) for v in np.geomspace(q_lo, qr_max, num=log_n)]
    lin_pts = [float(v) for v in np.linspace(qr_min, qr_max, num=lin_n)]
    pts = sorted(set([qr_min, qr_max] + log_pts + lin_pts))
    return pts


def _saved_qr_max_excluding_terminal_zero(entry):
    """Return max saved qr, ignoring a terminal zero-valued anchor if present."""
    xs = entry.get("distance_cache", []) if isinstance(entry, dict) else []
    rs = entry.get("value_cache_real", []) if isinstance(entry, dict) else []
    ims = entry.get("value_cache_imag", []) if isinstance(entry, dict) else []
    n = min(len(xs), len(rs), len(ims))
    if n <= 0:
        return 0.0

    idx = n - 1
    try:
        r_last = float(rs[idx])
        i_last = float(ims[idx])
        if idx > 0 and abs(r_last) <= 1e-30 and abs(i_last) <= 1e-30:
            idx -= 1
        return float(xs[idx])
    except Exception:
        return float(xs[-1]) if len(xs) > 0 else 0.0


def _write_adaptive_test_plot(
    key,
    qr_points,
    rel_diff_percent,
    manual_real,
    manual_imag,
    interp_real,
    interp_imag,
    threshold_percent,
    qr_max_effective,
    out_path,
):
    import matplotlib.pyplot as plt

    x = np.asarray(qr_points, dtype=float)
    rel = np.maximum(np.asarray(rel_diff_percent, dtype=float), 1e-16)
    man_re = np.asarray(manual_real, dtype=float)
    man_im = np.asarray(manual_imag, dtype=float)
    itp_re = np.asarray(interp_real, dtype=float)
    itp_im = np.asarray(interp_imag, dtype=float)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10.0, 6.8), sharex=True)
    ax0.semilogy(x, rel, color="#C21A1A", linewidth=1.3, label="rel diff [%]")
    ax0.axhline(float(threshold_percent), color="#1F5AA6", linestyle="--", linewidth=1.0, label="threshold")
    ax0.set_ylabel("Relative diff [%]")
    ax0.set_title("Adaptive Test Performance: {}".format(key))
    ax0.grid(True, which="both", alpha=0.28)
    ax0.legend(loc="upper right", fontsize=8)

    ax1.plot(x, man_re, color="#1E7A3B", linewidth=1.1, label="manual Re")
    ax1.plot(x, man_im, color="#1E7A3B", linewidth=1.1, linestyle="--", label="manual Im")
    ax1.plot(x, itp_re, color="#8A4A00", linewidth=1.1, label="pred Re")
    ax1.plot(x, itp_im, color="#8A4A00", linewidth=1.1, linestyle="--", label="pred Im")
    x_min = float(np.min(x)) if x.size else 0.0
    x_max = float(np.max(x)) if x.size else 0.0
    q_eff = float(qr_max_effective)
    q_line = min(max(q_eff, x_min), x_max)
    q_label = "qr_max_effective"
    if q_eff > x_max + max(1e-12, 1e-9 * max(1.0, x_max)):
        q_label = "qr_max_effective (off-plot)"
    elif q_eff < x_min - max(1e-12, 1e-9 * max(1.0, x_min)):
        q_label = "qr_max_effective (before plot)"
    ax1.axvline(q_line, color="#4C4C4C", linestyle=":", linewidth=1.0, label=q_label)
    ax1.set_xlabel("qr")
    ax1.set_ylabel("Value")
    ax1.grid(True, alpha=0.28)
    ax1.legend(loc="upper right", fontsize=8, ncol=2)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def test_precomputed_interpolators(
    threshold_percent=10.0,
    input_file="interpolators.json",
    qr_test_points=None,
    n0_values=range(1, 6),
    l_max=10,
    max_interpolators=25,
):
    """Compare I_integral directly against JSON-backed uniform-grid interpolator."""
    if qr_test_points is None:
        qr_test_points = _DEFAULT_TEST_QR

    precomputed_data = load_interpolators_file(input_file)
    threshold = float(threshold_percent) / 100.0
    combinations = _build_combinations(n0_values, l_max)
    if max_interpolators is not None:
        combinations = combinations[:max_interpolators]

    checked = 0
    flagged = []
    missing_cache = []

    for n0, l, m, typ, part in combinations:
        key = _interpolator_key(n0, l, m, typ, part)
        if key not in precomputed_data:
            missing_cache.append(key)
            continue

        interp_obj = I_interpolator(n0, l, m, typ, part, precomputed_cache=precomputed_data[key])
        for qr in qr_test_points:
            manual = I_integral((qr, l, m, typ, part))
            interpolated = interp_obj.interpolate(qr)
            manual_abs = abs(complex(manual))
            diff_abs = abs(complex(interpolated - manual))
            denom = max(manual_abs, 1e-14)
            rel_diff = diff_abs / denom
            checked += 1
            if rel_diff > threshold:
                flagged.append({
                    "key": key,
                    "qr": float(qr),
                    "relative_difference_percent": 100.0 * rel_diff,
                    "manual": complex(manual),
                    "interpolated": complex(interpolated),
                    "manual_abs": manual_abs,
                    "interpolated_abs": abs(complex(interpolated)),
                    "difference_abs": diff_abs,
                })

    print("Interpolator JSON tester summary:")
    print("  checked points:", checked)
    print("  flagged points (> {:.3f}%): {}".format(threshold_percent, len(flagged)))
    print("  missing cache entries:", len(missing_cache))
    for item in flagged[:20]:
        print(
            "FLAG", item["key"],
            "qr=", item["qr"],
            "reldiff%=", "{:.4f}".format(item["relative_difference_percent"]),
            "manual=", item["manual"],
            "interpolated=", item["interpolated"],
            "|manual|=", "{:.6e}".format(item["manual_abs"]),
            "|interp|=", "{:.6e}".format(item["interpolated_abs"]),
            "|diff|=", "{:.6e}".format(item["difference_abs"]),
        )
    if len(flagged) > 20:
        print("... {} more flagged points omitted".format(len(flagged) - 20))
    if missing_cache:
        print("First missing cache key:", missing_cache[0])
    return {
        "checked_points": checked,
        "flagged_count": len(flagged),
        "missing_cache_count": len(missing_cache),
        "flagged": flagged,
        "missing_cache": missing_cache,
    }


def _worker_test_adaptive_channel(args):
    """Evaluate one adaptive interpolator channel against direct integration.

    Returns (key, qr_points, rel_series, manual_real, manual_imag,
             interp_real, interp_imag, flagged_items).
    All values are plain Python types so the result is picklable.
    """
    n0, l, m, typ, part, entry, qr_points, threshold, use_point_cache, point_cache_dir = args
    key = _interpolator_key(n0, l, m, typ, part)
    interp_obj = I_adaptive_interpolator(n0, l, m, typ, part, precomputed_cache=entry)
    point_store = _ChannelIntegralPointStore(
        cache_dir=point_cache_dir,
        n0=n0,
        l=l,
        m=m,
        typ=typ,
        part=part,
        enabled=bool(use_point_cache),
    )
    rel_series = []
    manual_real_series = []
    manual_imag_series = []
    interp_real_series = []
    interp_imag_series = []
    flagged = []
    for qr in qr_points:
        manual_c = complex(point_store.evaluate(qr))
        interp_c = complex(interp_obj.interpolate(qr))
        manual_abs = abs(manual_c)
        interp_abs = abs(interp_c)
        diff_abs = abs(interp_c - manual_c)
        rel_diff = diff_abs / max(manual_abs, 1e-14)
        rel_pct = 100.0 * rel_diff
        rel_series.append(rel_pct)
        manual_real_series.append(manual_c.real)
        manual_imag_series.append(manual_c.imag)
        interp_real_series.append(interp_c.real)
        interp_imag_series.append(interp_c.imag)
        if rel_diff > threshold:
            flagged.append({
                "key": key,
                "qr": float(qr),
                "relative_difference_percent": rel_pct,
                "manual": manual_c,
                "interpolated": interp_c,
                "manual_abs": manual_abs,
                "interpolated_abs": interp_abs,
                "difference_abs": diff_abs,
            })
    point_store.flush()
    return key, qr_points, rel_series, manual_real_series, manual_imag_series, interp_real_series, interp_imag_series, flagged


def test_precomputed_adaptive_interpolators(
    threshold_percent=1.0,
    input_file="interpolator_adaptive.json",
    qr_test_points=None,
    n0_values=range(1, 6),
    l_max=10,
    max_interpolators=25,
    plot_output_dir=None,
    plot_points=120,
    plot_qr_max=None,
    num_workers=None,
    use_point_cache=True,
    point_cache_dir="integral_point_cache",
):
    """Compare I_integral directly against JSON-backed adaptive PCHIP interpolator."""
    from multiprocessing import Pool, cpu_count

    if qr_test_points is None and plot_output_dir is None:
        qr_test_points = _DEFAULT_TEST_QR

    precomputed_data = load_adaptive_interpolators_file(input_file)
    threshold = float(threshold_percent) / 100.0
    combinations = _build_combinations(n0_values, l_max)
    if max_interpolators is not None:
        combinations = combinations[:max_interpolators]

    if plot_output_dir:
        os.makedirs(plot_output_dir, exist_ok=True)

    missing_cache = []
    tasks = []
    task_meta = []  # (n0, l, m, typ, part, entry) parallel to tasks

    for n0, l, m, typ, part in combinations:
        key = _interpolator_key(n0, l, m, typ, part)
        if key not in precomputed_data:
            missing_cache.append(key)
            continue
        entry = precomputed_data[key]
        if plot_output_dir:
            qmin_local = float(entry.get("distance_cache", [0.0])[0]) if entry.get("distance_cache") else 0.0
            if plot_qr_max is not None:
                qmax_local = max(qmin_local, float(plot_qr_max))
            else:
                saved_qr_max = _saved_qr_max_excluding_terminal_zero(entry)
                qmax_local = max(qmin_local, 1.5 * float(saved_qr_max))
            channel_qr_points = _build_adaptive_test_qr_points(qmin_local, qmax_local, plot_points)
        else:
            channel_qr_points = list(qr_test_points)
        tasks.append((n0, l, m, typ, part, entry, channel_qr_points, threshold, bool(use_point_cache), str(point_cache_dir)))
        task_meta.append((n0, l, m, typ, part, entry))

    total = len(tasks)
    workers = num_workers if num_workers is not None else cpu_count()
    print("Testing {} adaptive interpolator channels with {} workers ...".format(total, workers))

    checked = 0
    flagged = []
    per_channel_summary = []
    # Results keyed by channel key for ordered output after parallel collection.
    channel_results = {}

    with Pool(processes=workers) as pool:
        for done, result in enumerate(pool.imap_unordered(_worker_test_adaptive_channel, tasks), start=1):
            key, qr_points, rel_series, manual_real, manual_imag, interp_real, interp_imag, ch_flagged = result
            channel_results[key] = (qr_points, rel_series, manual_real, manual_imag, interp_real, interp_imag, ch_flagged)
            checked += len(rel_series)
            flagged.extend(ch_flagged)
            status = "FLAGGED({})".format(len(ch_flagged)) if ch_flagged else "OK"
            print("  [{}/{}] {} pts={} max_reldiff%={:.3e} {}".format(
                done, total, key, len(rel_series),
                max(rel_series) if rel_series else 0.0, status))

    # Build per-channel summary and plots (plots must be in main process).
    for n0, l, m, typ, part, entry in task_meta:
        key = _interpolator_key(n0, l, m, typ, part)
        if key not in channel_results:
            continue
        qr_points, rel_series, manual_real, manual_imag, interp_real, interp_imag, ch_flagged = channel_results[key]
        if rel_series:
            per_channel_summary.append({
                "key": key,
                "points": int(len(rel_series)),
                "max_rel_diff_percent": float(max(rel_series)),
                "median_rel_diff_percent": float(np.median(np.array(rel_series))),
                "flagged_points": int(len(ch_flagged)),
                "qr_max_effective": float(entry.get("qr_max_effective", 0.0)),
            })
        if plot_output_dir and rel_series:
            plot_name = "adaptive_test_{}.png".format(_safe_plot_filename(key))
            plot_path = os.path.join(plot_output_dir, plot_name)
            _write_adaptive_test_plot(
                key=key,
                qr_points=qr_points,
                rel_diff_percent=rel_series,
                manual_real=manual_real,
                manual_imag=manual_imag,
                interp_real=interp_real,
                interp_imag=interp_imag,
                threshold_percent=threshold_percent,
                qr_max_effective=float(entry.get("qr_max_effective", 0.0)),
                out_path=plot_path,
            )

    print("Adaptive interpolator JSON tester summary:")
    print("  checked points:", checked)
    print("  flagged points (> {:.3f}%): {}".format(threshold_percent, len(flagged)))
    print("  missing cache entries:", len(missing_cache))
    for item in flagged[:20]:
        print(
            "FLAG", item["key"],
            "qr=", item["qr"],
            "reldiff%=", "{:.4f}".format(item["relative_difference_percent"]),
            "manual=", item["manual"],
            "interpolated=", item["interpolated"],
            "|manual|=", "{:.6e}".format(item["manual_abs"]),
            "|interp|=", "{:.6e}".format(item["interpolated_abs"]),
            "|diff|=", "{:.6e}".format(item["difference_abs"]),
        )
    if len(flagged) > 20:
        print("... {} more flagged points omitted".format(len(flagged) - 20))
    if missing_cache:
        print("First missing cache key:", missing_cache[0])

    if plot_output_dir and per_channel_summary:
        summary_path = os.path.join(plot_output_dir, "adaptive_test_summary.json")
        with open(summary_path, "w") as out:
            json.dump({
                "input_file": input_file,
                "threshold_percent": float(threshold_percent),
                "channels": per_channel_summary,
            }, out, indent=2)
        print("  wrote plots to:", plot_output_dir)
        print("  wrote plot summary:", summary_path)

    return {
        "checked_points": checked,
        "flagged_count": len(flagged),
        "missing_cache_count": len(missing_cache),
        "flagged": flagged,
        "missing_cache": missing_cache,
        "channel_summary": per_channel_summary,
        "plot_output_dir": plot_output_dir,
    }
