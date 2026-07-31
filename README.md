# HiFAST

**HiFAST** is a pipeline designed for the calibration and imaging of HI (neutral hydrogen) data from the Five-hundred-meter Aperture Spherical radio Telescope (FAST).

Documentation: https://hifast.readthedocs.io/

## Experimental quasi-periodic bump mask

The experimental `rfi_bump` command detects the approximately 1 MHz
quasi-periodic narrow bump RFI in pre-Doppler `flux` spectra. It reads the
input file without modifying it and writes accepted mask samples as `NaN` in
a new file:

```bash
python -m hifast.rfi_bump \
  example-flux.hdf5 \
  --outdir output
```

The output name is `example-flux-rfi_bump.hdf5`. XX and YY are located and
masked independently. Localization failures and weak detections are skipped;
pipeline errors stop the command. The method remains experimental because its
hard mask has not passed the science-signal protection test.

The command and file I/O are implemented in `hifast/rfi_bump.py`; detection
and mask functions are implemented in `hifast/core/rfi_bump.py`. The frozen
v5.3 helper programs and configurations are private package resources under
`hifast/rfi_bump_support/`.
