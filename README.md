# Quantum Optomechanics Papers

> A Collection of Solved Papers using [The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom).

## List of Solved Papers

* J. Opt. Soc. Am. B **33**, 1335 (2016)
* [New J. Phys. **22**, 013049](https://doi.org/10.1088/1367-2630/ab6522) (2020)
* [New J. Phys. **22**, 063041](https://doi.org/10.1088/1367-2630/ab90d2) (2020)
* [Opt. Lett. **41**, 2676](https://doi.org/10.1364/OL.41.002676) (2016)
* [Phys. Rev. A **101**, 053836](https://doi.org/10.1103/PhysRevA.101.053836) (2020)
* [Phys. Rev. Lett. **103**, 213603](https://doi.org/10.1103/PhysRevLett.103.213603) (2009)
* [Phys. Rev. Lett. **111**, 103605](https://doi.org/10.1103/PhysRevLett.111.103605) (2013)
* [Phys. Rev. Lett. **114**, 013601](https://doi.org/10.1103/PhysRevLett.114.013601) (2015)
* [Phys. Rev. Lett. **119**, 153901](https://doi.org/10.1103/PhysRevLett.119.153901) (2017)

## Structure of the Repository

```
ROOT_DIR/
|
├───gui_templates/
│   ├───foo/
│   │   ├───bar.py
│   │   └───...
│   └───...
|
├───notebooks/
│   ├───foo/
│   │   ├───bar.ipynb
│   │   └───...
│   └───...
|
│───scripts/
│   ├───foo/
│   │   ├───bar.py
│   │   └───...
│   └───...
|
├───systems/
│   ├───__init__.py
│   ├───Foo.py
│   └───...
│
├───.gitignore
├───CHANGELOG.md
├───GUI.py
└───README.md
```

## Execution

### Installing Dependencies

The project requires `Python 3.8+` installed via the [Anaconda distribution](https://www.anaconda.com/products/individual). 
An extensive guide to set up your python environment same can be found [here](https://sampreet.github.io/python-for-physicists/modules/m01-getting-started/m01t01-setting-up-python.html).

Once the installation is complete and `conda` is configured, it is preferable to create a new conda environment (say `qom`) and activate it using:

```bash
conda create -n qom python=3
conda activate qom
```

This project uses [The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom) via Python Package Index using `pip`:

```bash
pip install -i https://test.pypi.org/simple/ qom
```

Alternatively, [clone the repository](https://github.com/Sampreet/qom) or [download the sources](https://github.com/Sampreet/qom/archive/refs/heads/master.zip) as `.zip` and extract the contents.
Now, execute the following from *outside* the top-level directory, `ROOT_DIR`, inside which `setup.py` is located:

```bash
pip install -e ROOT_DIR
```

### Running the Scripts

To run the scripts, navigate *inside* the top-level directory, `ROOT_DIR`, and execute:

```bash
python scripts/foo_bar.py
```

Here, `foo_bar.py` is the name of the script.

To run in GUI mode using `PowerShell` or `bash`, navigate to `ROOT_DIR` and execute:

```bash
python -c 'from qom.ui.gui import run; run()'
```

Alternatively, run `GUI.py` from within the directory.