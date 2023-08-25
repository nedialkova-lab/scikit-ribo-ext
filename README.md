# scikit-ribo-ext
Extended and improved scikit-ribo (https://github.com/schatzlab/scikit-ribo) for accurate estimation and robust modelling of translation dynamics at codon resolution

[![Documentation Status](https://readthedocs.org/projects/scikit-ribo/badge/?version=latest)](http://scikit-ribo.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Please see original [paper](https://doi.org/10.1016/j.cels.2017.12.007), GitHub [repo](https://github.com/schatzlab/scikit-ribo/tree/master) and [documentation](https://scikit-ribo.readthedocs.io/en/latest/?badge=latest) for details on scikit-ribo method, functionality, and usage.

## New features and concepts


## Create environment from supplied yml file
```bash
conda env create -f scikit-ribo-ext.yml
```

## Activate environment
```bash
conda activate scikit-ribo-ext
```

## Running model fitting
**Critical:** Please run the script `scikit-ribo-run.py` that is located inside the `scikit-ribo-ext` directory of your cloned GitHub repository (for example, ```/path/to/repo/scikit-ribo-ext/scikit-ribo-ext/scikit-ribo-run.py```). This is crucial as inside your environment you also have an installed version of the original pipeline which can be run by simply calling `scikit-ribo-run.py` without specifying the path to this version. This will run the original pipeline and will be missing the new features and improvements listed above.