import unittest

import numpy as np

from hifast.bld import IO, create_parser
from hifast.core.baseline import BL_ContextPLS, get_baseline
from hifast.core.baseline_context import (
    ContextReweighter,
    asym2_weights_and_noise,
    mhz_to_odd_window,
    validate_frequency_axis,
)


class TestContextWeights(unittest.TestCase):
    def test_mhz_window_is_odd_for_both_frequency_directions(self):
        ascending = np.linspace(1400.0, 1410.0, 101)
        descending = ascending[::-1]
        self.assertEqual(mhz_to_odd_window(ascending, 1.0), 11)
        self.assertEqual(mhz_to_odd_window(descending, 1.0), 11)

    def test_frequency_must_be_strictly_monotonic(self):
        with self.assertRaisesRegex(ValueError, "strictly monotonic"):
            validate_frequency_axis([1400.0, 1400.1, 1400.05])

    def test_excluded_value_does_not_affect_context_update(self):
        freq = np.linspace(1400.0, 1410.0, 2001)
        rng = np.random.default_rng(3)
        residual = rng.normal(0.0, 0.02, freq.size)
        exclude = np.zeros(freq.size, dtype=bool)
        exclude[1000] = True
        changed = residual.copy()
        changed[1000] = 1e9

        reweighter = ContextReweighter(freq)
        weights_a, shift_a, probability_a = reweighter.update(
            residual, exclude=exclude
        )
        weights_b, shift_b, probability_b = reweighter.update(
            changed, exclude=exclude
        )

        np.testing.assert_allclose(weights_a, weights_b, atol=0, rtol=0)
        np.testing.assert_allclose(probability_a, probability_b, atol=0, rtol=0)
        self.assertEqual(shift_a, shift_b)
        self.assertEqual(weights_a[1000], 0.0)

    def test_expansion_only_increases_detection_probability(self):
        freq = np.linspace(1400.0, 1410.0, 2001)
        residual = np.zeros(freq.size)
        residual += 0.08 * np.exp(-0.5 * ((freq - 1405.0) / 0.35) ** 2)
        valid = np.ones(freq.size, dtype=bool)
        plain = ContextReweighter(freq, expand_mhz=0.0)
        expanded = ContextReweighter(freq, expand_mhz=0.5)
        _, sigma = asym2_weights_and_noise(residual, valid)
        plain_probability = plain._probability(
            plain._contrasts(residual, valid), sigma
        )
        expanded_probability = expanded._probability(
            expanded._contrasts(residual, valid), sigma
        )
        self.assertTrue(np.all(expanded_probability >= plain_probability))


class TestContextPLS(unittest.TestCase):
    @staticmethod
    def _reference_spectrum():
        freq = np.linspace(1400.0, 1420.0, 4097)
        rng = np.random.default_rng(12)
        spectrum = (
            0.04 * np.sin((freq - 1400.0) * np.pi / 10.0)
            + 0.025 * rng.normal(size=freq.size)
            + 0.06 * np.exp(-0.5 * ((freq - 1410.0) / 0.8) ** 2)
        )
        return freq, spectrum

    def test_matches_frozen_experimental_result(self):
        freq, spectrum = self._reference_spectrum()
        fitter = BL_ContextPLS()
        baseline = fitter.fit(x=freq, y=spectrum)
        indices = np.array([0, 128, 1024, 2048, 3072, 4096])
        expected = np.array([
            0.0030017260645885,
            0.009867951475815578,
            0.039191119534411276,
            0.002342640967814209,
            -0.039461928141893254,
            -0.009909831988152604,
        ])
        np.testing.assert_allclose(
            baseline[indices], expected, rtol=0, atol=1e-7
        )
        self.assertEqual(fitter.i + 1, 6)

    def test_descending_frequency_reverses_the_same_fit(self):
        freq, spectrum = self._reference_spectrum()
        ascending = BL_ContextPLS().fit(x=freq, y=spectrum)
        descending = BL_ContextPLS().fit(
            x=freq[::-1], y=spectrum[::-1]
        )
        np.testing.assert_allclose(
            ascending, descending[::-1], rtol=0, atol=2e-8
        )

    def test_frequency_averaging_uses_the_averaged_grid(self):
        freq, spectrum = self._reference_spectrum()
        baseline = get_baseline(
            freq,
            spectrum,
            average_every=4,
            method='PLS-context',
            bl_para={
                'lam': 1e9,
                'deg': 2,
                'offset': 2,
                'ratio': 0.01,
                'niter': 100,
                'rew': True,
            },
        )
        self.assertEqual(baseline.shape, spectrum.shape)
        self.assertTrue(np.all(np.isfinite(baseline)))

    def test_masked_nan_does_not_break_fit(self):
        freq, spectrum = self._reference_spectrum()
        spectrum = spectrum.copy()
        spectrum[2000:2010] = np.nan
        exclude = np.zeros(spectrum.size, dtype=bool)
        exclude[2000:2010] = True
        baseline = get_baseline(
            freq,
            spectrum,
            exclude=exclude,
            method='PLS-context',
            bl_para={
                'lam': 1e9,
                'deg': 2,
                'offset': 2,
                'ratio': 0.01,
                'niter': 100,
                'rew': True,
            },
        )
        self.assertTrue(np.all(np.isfinite(baseline)))

    def test_polarizations_are_fitted_independently(self):
        freq, spectrum = self._reference_spectrum()
        other = spectrum[::-1].copy()
        values = np.stack((spectrum, other), axis=-1)[None, :, :]
        baseline_a = get_baseline(
            freq,
            values,
            axis=1,
            method='PLS-context',
            bl_para={
                'lam': 1e9,
                'deg': 2,
                'offset': 2,
                'ratio': 0.01,
                'niter': 100,
                'rew': True,
            },
        )
        values[0, :, 1] += 10.0 * np.exp(
            -0.5 * ((freq - 1408.0) / 0.3) ** 2
        )
        baseline_b = get_baseline(
            freq,
            values,
            axis=1,
            method='PLS-context',
            bl_para={
                'lam': 1e9,
                'deg': 2,
                'offset': 2,
                'ratio': 0.01,
                'niter': 100,
                'rew': True,
            },
        )
        np.testing.assert_allclose(
            baseline_a[0, :, 0], baseline_b[0, :, 0], atol=0, rtol=0
        )

    def test_get_baseline_rejects_non_frequency_coordinates(self):
        values = np.zeros(5)
        with self.assertRaisesRegex(ValueError, "strictly monotonic"):
            get_baseline(
                np.array([0.0, 1.0, 1.0, 2.0, 3.0]),
                values,
                method='PLS-context',
                bl_para={
                    'lam': 1e9,
                    'deg': 2,
                    'offset': 2,
                    'ratio': 0.01,
                    'niter': 2,
                    'rew': True,
                },
            )

    def test_command_line_exposes_context_method(self):
        args = create_parser().parse_args([
            'input.hdf5', '--method', 'PLS-context'
        ])
        self.assertEqual(args.method, 'PLS-context')
        self.assertIsNone(args.lam)

    def test_bld_rejects_context_on_time_axis(self):
        args = create_parser().parse_args([
            'input.hdf5', '--method', 'PLS-context', '--trans', 'True'
        ])
        spectra = np.zeros((2, 9, 1))
        with self.assertRaisesRegex(ValueError, "cannot fit along the time axis"):
            IO.fit_baseline(
                spectra,
                np.linspace(1400.0, 1401.0, 9),
                np.array([60000.0, 60000.1]),
                args,
            )


if __name__ == '__main__':
    unittest.main()
