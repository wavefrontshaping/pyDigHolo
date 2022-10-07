# pyDigHolo
Python module to use [Joel Carpenter](https://github.com/joelacarpenter)'s [digHolo](https://github.com/joelacarpenter/digHolo)  C++ library.
[digHolo](https://github.com/joelacarpenter/digHolo) is a high-speed library for off-axis digital holography.

## Installation

First, install 
Then, install pyDigHolo, using `pip`

```bash
pip install git+https://github.com/wavefrontshaping/pyDigHolo.git
```

or clone the repository and install using in the module folder

```bash
python setup.py install
```
## Documentation 

This module more or less simply wraps the function from the `digHolo` C module, 
please refer to the offocial  [digHolo repository](https://github.com/joelacarpenter/digHolo) for documentation. 
## Basic usage

### Import the module

```python
from pyDigHolo import digHolo
```

### Instantiate the digHolo object

You need the provide the file path to fetch the `.dll` from `digHolo`.
If you copied the `.dll` into your working folder, simply use

```python
dh = digHolo('digHolo.dll')
```

### Set the experimental parameters

```python
dh.ConfigOffAxis(
    [frameWidth, frameHeight],
    [nx, ny], 
    resolutionMode,
    pixelSize,
    lambda0, 
    maxMG,
    polCount
)
```

### Configure auto-align procedure

Specify what you want the algorithm to automatically calibrate.
For instance:

```python
dh.ConfigSetAutoAlign(
    enable_align_beam_centre=True,
    enable_align_defocus=False,
    enable_align_tilt=True,
    enable_align_basis_waist=True,
    enable_align_fourier_win_radius=False,
)
```

### Give `pyDigHolo` a batch of interferograms to treat

```python
dh.SetBatch(frameCount, frames, dataType = 'Python')
```

### Perform the automatic calibration procedure

```python
dh.AutoAlign()
```

### Provide the intensity pattern of the reference

```python
dh.SetRefCalibrationIntensity(ref)
```

### Save the calibration

```python
dh.SaveConfig('config_pola1.npz')
```

### Load the calibration

Next time you want to use the off-axis procedure, 
if the system did not change, 
you can simply load a previous calibration file using:

```python
dh.Loadonfig('config_pola1.npz')
```

### Process a batch of data
It will treat the interferograms provided using `SetBatch()`.


```python
dh.ProcessBatch()
```

### Retrieve the complex field maps

```python
fields = dh.GetFields()
```




## Examples

See the two jupyter notebook:

* [Example 1: Simulated data](/examples/simulated_data.ipynb). 
  Generate test interferograms using digHolo and recover the complex fields.

* [Example 2: Experimental data ](/examples/experimental_data.ipynb). 
  Load experimentally measured with a non uniform reference intensity and recover the complex fields.