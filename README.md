# Quantum Optomechanics Papers

> A collection of solved papers using [The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom).

## Table of Solved Papers

paper link  | paper name    | qom system    | qom scripts
--------    | --------      | --------      | --------      
[New J. Phys. **22**, 013049](https://doi.org/10.1088/1367-2630/ab6522) | Nonlinear Dynamics of Weakly Dissipative Optomechanical Systems | [NewJPhys_22_013049](./systems/NewJPhys_22_013049.py) | [newjphys_22_013049](./scripts/newjphys_22_013049/)
[New J. Phys. **22**, 063041](https://doi.org/10.1088/1367-2630/ab90d2) | Stationary Quantum Entanglement between a Massive Mechanical Membrane and a Low Frequency LC Circuit | [NewJPhys_22_063041](./systems/NewJPhys_22_063041.py) | [newjphys_22_063041](./scripts/newjphys_22_063041/)
[Opt. Lett. **41**, 2676](https://doi.org/10.1364/OL.41.002676) | Solitons in Optomechanical Arrays | [OptLett_41_2676](./systems/OptLett_41_2676.py) | [optlett_41_2676](./scripts/optlett_41_2676/)
[Phys. Rev. A **101**, 053836](https://doi.org/10.1103/PhysRevA.101.053836) | Strong Mechanical Squeezing in a Standard Optomechanical System by Pump Modulation | [PhysRevA_101_053836](./systems/PhysRevA_101_053836.py) | [physreva_101_053836](./scripts/physreva_101_053836/)
[Phys. Rev. Lett. **103**, 213603](https://doi.org/10.1103/PhysRevLett.103.213603) | Gently Modulating Optomechanical Systems | [PhysRevLett_103_213603](./systems/PhysRevLett_103_213603.py) | [physrevlett_103_213603](./scripts/physrevlett_103_213603/)
[Phys. Rev. Lett. **111**, 103605](https://doi.org/10.1103/PhysRevLett.111.103605) | Measures of Quantum Synchronization in Continuous Variable Systems | [PhysRevLett_111_103605](./systems/PhysRevLett_111_103605.py) | [physrevlett_111_103605](./scripts/physrevlett_111_103605/)
[Phys. Rev. Lett. **114**, 013601](https://doi.org/10.1103/PhysRevLett.114.013601) | Route to Chaos in Optomechanics | [PhysRevLett_114_013601](./systems/PhysRevLett_114_013601.py) | [physrevlett_114_013601](./scripts/physrevlett_114_013601/)
[Phys. Rev. Lett. **119**, 153901](https://doi.org/10.1103/PhysRevLett.119.153901) | Kuznetsov-Ma Soliton Dynamics based on the Mechanical Effect of Light | [PhysRevLett_119_153901](./systems/PhysRevLett_119_153901.py) | [physrevlett_119_153901](./scripts/physrevlett_119_153901/)

## Structure of the Repository

```
ROOT_DIR/
|
├───gui_templates/
│   ├───foobar_123/
│   │   ├───baz.py
│   │   └───...
│   └───...
|
├───notebooks/
│   ├───foobar_123/
│   │   ├───baz.ipynb
│   │   └───...
│   └───...
|
│───scripts/
│   ├───foobar_123/
│   │   ├───baz.py
│   │   └───...
│   └───...
|
├───systems/
│   ├───__init__.py
│   ├───FooBar_123.py
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
python scripts/foobar_123/baz.py
```

Here, `foobar_123` is the name of the folder and `baz.py` is the name of the file.

To run in GUI mode using `PowerShell` or `bash`, navigate to `ROOT_DIR` and execute:

```bash
python -c 'from qom.ui.gui import run; run()'
```

Alternatively, run `GUI.py` from within the directory.