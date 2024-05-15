# Meteor camera auto calibration tool

## Requirements

For Ubuntu, install `libgirepository1.0-dev`, `astrometry.net` and `imagemagick`.

Create a Python virtual environment and install requirements:
```
pip install -r requirements.txt
```

Then install the RMS package separately:
```
pip install git+https://github.com/CroatianMeteorNetwork/RMS.git
```

Then, clone the amscams repository:

```
git clone https://github.com/pemo1065/amscams.git
```

## Calibration

Run the following command:

```
python initial_calibration.py calibrate -t "2024-04-27 22:04:26" -i /path/to/image.png -m /path/to/mask.png -c <calibration_method> -l <lim-magnitude>
```

`<calibration-mehod>` can be `nmn`, `rms` or `ams`.
