"""Context-aware weights for PLS spectral baseline fitting.

The Context method supplements point-wise asym2 weights with a low-cost,
dual-scale test for coherent positive residuals.  Window widths are expressed
in MHz, so this module is intended for fits along a frequency axis.
"""

from collections import namedtuple

import numpy as np
from scipy.ndimage import maximum_filter1d
from scipy.special import expit


_EPSILON = np.finfo(np.float64).eps


def validate_frequency_axis(freq):
    """Return a validated one-dimensional monotonic frequency array."""
    values = np.asarray(freq, dtype=np.float64)
    if values.ndim != 1 or values.size < 3:
        raise ValueError("PLS-context requires at least three frequency channels")
    if not np.all(np.isfinite(values)):
        raise ValueError("PLS-context requires finite frequency coordinates")
    differences = np.diff(values)
    if not (np.all(differences > 0) or np.all(differences < 0)):
        raise ValueError("PLS-context requires a strictly monotonic frequency axis")
    return values


def mhz_to_odd_window(freq, width_mhz):
    """Convert a physical full window width to an odd channel count."""
    if width_mhz <= 0:
        raise ValueError("Context window widths must be positive")
    values = validate_frequency_axis(freq)
    spacing = abs(float(np.median(np.diff(values))))
    window = max(1, int(round(float(width_mhz) / spacing)))
    return window if window % 2 else window + 1


def _window_geometry(size, window):
    half = window // 2
    index = np.arange(size)
    left = np.maximum(0, index - half)
    right = np.minimum(size, index + half + 1)
    return left, right


def _window_mean(cumulative_values, geometry, count):
    left, right = geometry
    mean = np.zeros(count.size, dtype=np.float64)
    total = cumulative_values[right] - cumulative_values[left]
    np.divide(total, count, out=mean, where=count > 0)
    return mean


def _robust_scale(values):
    if values.size == 0:
        return _EPSILON
    center = float(np.median(values))
    scale = 1.4826 * float(np.median(np.abs(values - center)))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.std(values))
    return max(scale, _EPSILON)


def asym2_weights_and_noise(residual, valid, offset=2.0):
    """Return exclusion-aware asym2 weights and a full-noise estimate."""
    residual = np.asarray(residual, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    negative = residual[valid & (residual < 0)]
    if negative.size >= 2:
        mean = float(np.mean(negative))
        scale = max(float(np.std(negative)), _EPSILON)
        noise_sigma = max(float(np.hypot(mean, scale)), _EPSILON)
    else:
        usable = residual[valid]
        scale = _robust_scale(usable)
        mean = float(np.mean(negative)) if negative.size else -scale
        noise_sigma = scale

    threshold = float(offset) * scale - mean
    weights = expit(-2.0 * (np.abs(residual) - threshold) / scale)
    weights[~valid] = 0.0
    return weights, noise_sigma


_ContextContrast = namedtuple(
    '_ContextContrast', ['contrasts', 'variance_scales', 'usable']
)


class ContextReweighter:
    """Compute the tested dual-scale Context update for one frequency grid."""

    def __init__(
        self,
        freq,
        *,
        small_mhz=(1.0, 2.0),
        large_mhz=(4.0, 8.0),
        thresholds=(2.5, 3.0),
        slope=2.0,
        expand_mhz=0.5,
        linefree_probability_max=0.9,
        offset=2.0,
    ):
        self.freq = validate_frequency_axis(freq)
        self.small_mhz = tuple(float(value) for value in small_mhz)
        self.large_mhz = tuple(float(value) for value in large_mhz)
        self.thresholds = tuple(float(value) for value in thresholds)
        if not (
            len(self.small_mhz)
            == len(self.large_mhz)
            == len(self.thresholds)
            and len(self.thresholds) > 0
        ):
            raise ValueError("Context scales and thresholds must have equal lengths")
        if slope <= 0:
            raise ValueError("Context logistic slope must be positive")
        if expand_mhz < 0:
            raise ValueError("Context expansion width cannot be negative")
        if not 0 <= linefree_probability_max < 1:
            raise ValueError("Context line-free probability limit must be in [0, 1)")

        self.slope = float(slope)
        self.expand_mhz = float(expand_mhz)
        self.linefree_probability_max = float(linefree_probability_max)
        self.offset = float(offset)
        self._geometry = []
        for small_width, large_width in zip(self.small_mhz, self.large_mhz):
            small_window = mhz_to_odd_window(self.freq, small_width)
            large_window = mhz_to_odd_window(self.freq, large_width)
            if large_window <= small_window:
                raise ValueError("Each large Context window must exceed its small window")
            self._geometry.append(
                (
                    _window_geometry(self.freq.size, small_window),
                    _window_geometry(self.freq.size, large_window),
                )
            )
        self._expand_window = (
            mhz_to_odd_window(self.freq, self.expand_mhz)
            if self.expand_mhz > 0
            else 1
        )
        self._cached_valid = None
        self._cached_mask_geometry = None

    def _prepare_mask_geometry(self, valid):
        if (
            self._cached_valid is not None
            and np.array_equal(self._cached_valid, valid)
        ):
            return self._cached_mask_geometry

        cumulative_count = np.empty(valid.size + 1, dtype=np.int64)
        cumulative_count[0] = 0
        np.cumsum(valid, dtype=np.int64, out=cumulative_count[1:])
        mask_geometry = []
        for small_geometry, large_geometry in self._geometry:
            small_left, small_right = small_geometry
            large_left, large_right = large_geometry
            small_count = (
                cumulative_count[small_right] - cumulative_count[small_left]
            )
            large_count = (
                cumulative_count[large_right] - cumulative_count[large_left]
            )
            usable = (small_count >= 2) & (large_count > small_count)
            variance = (
                1.0 / np.maximum(small_count, 1)
                - 1.0 / np.maximum(large_count, 1)
            )
            variance = np.sqrt(np.maximum(variance, _EPSILON))
            mask_geometry.append(
                (small_count, large_count, variance, usable)
            )
        self._cached_valid = np.copy(valid)
        self._cached_mask_geometry = tuple(mask_geometry)
        return self._cached_mask_geometry

    def _contrasts(self, residual, valid):
        contrasts = []
        variance_scales = []
        usable_scales = []
        cumulative_values = np.empty(residual.size + 1, dtype=np.float64)
        cumulative_values[0] = 0.0
        np.cumsum(
            np.where(valid, residual, 0.0), out=cumulative_values[1:]
        )
        mask_geometry = self._prepare_mask_geometry(valid)
        for (
            (small_geometry, large_geometry),
            (small_count, large_count, variance, usable),
        ) in zip(self._geometry, mask_geometry):
            small_mean = _window_mean(
                cumulative_values, small_geometry, small_count
            )
            large_mean = _window_mean(
                cumulative_values, large_geometry, large_count
            )
            contrasts.append(small_mean - large_mean)
            variance_scales.append(variance)
            usable_scales.append(usable)
        return _ContextContrast(
            tuple(contrasts), tuple(variance_scales), tuple(usable_scales)
        )

    def _probability(self, contrast, noise_sigma):
        probability = np.zeros(self.freq.size, dtype=np.float64)
        sigma = max(float(noise_sigma), _EPSILON)
        for values, variance, usable, threshold in zip(
            contrast.contrasts,
            contrast.variance_scales,
            contrast.usable,
            self.thresholds,
        ):
            score = values / (sigma * variance)
            scale_probability = expit(self.slope * (score - threshold))
            scale_probability[~usable] = 0.0
            probability = np.maximum(probability, scale_probability)
        if self._expand_window > 1:
            probability = maximum_filter1d(
                probability, size=self._expand_window, mode="nearest"
            )
        return np.clip(probability, 0.0, 1.0)

    def update(self, residual, exclude=None):
        """Return new weights, the baseline shift, and signal probability."""
        values = np.asarray(residual, dtype=np.float64)
        if values.ndim != 1 or values.size != self.freq.size:
            raise ValueError("Context residual must match the frequency axis")
        if exclude is None:
            valid = np.ones(values.size, dtype=bool)
        else:
            excluded = np.asarray(exclude, dtype=bool)
            if excluded.shape != values.shape:
                raise ValueError("Context exclusion mask must match the residual")
            valid = ~excluded
        if np.count_nonzero(valid) < 3:
            raise ValueError("PLS-context requires at least three usable channels")

        _, initial_sigma = asym2_weights_and_noise(
            values, valid, offset=self.offset
        )
        contrast = self._contrasts(values, valid)
        initial_probability = self._probability(contrast, initial_sigma)
        line_free = valid & (
            initial_probability <= self.linefree_probability_max
        )
        center_shift = (
            float(np.median(values[line_free]))
            if np.count_nonzero(line_free) >= 3
            else 0.0
        )

        centered = values - center_shift
        point_weights, centered_sigma = asym2_weights_and_noise(
            centered, valid, offset=self.offset
        )
        # A constant shift cancels in small_mean - large_mean.  Reusing the
        # contrasts avoids a second set of moving-window calculations while
        # retaining the shifted residual's noise normalization.
        probability = self._probability(contrast, centered_sigma)
        weights = point_weights * (1.0 - probability)
        weights[~valid] = 0.0
        return weights, center_shift, probability
