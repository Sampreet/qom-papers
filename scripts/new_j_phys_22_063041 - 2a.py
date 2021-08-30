# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import NewJPhys22_063041

# all parameters
params = {
    'looper': {
        'X': {
            'var': 'Delta_norm',
            'min': 0.4,
            'max': 1.6,
            'dim': 121
        },
        'Y': {
            'var': 'omega_LC_norm',
            'min': 0.4,
            'max': 1.6,
            'dim': 121
        }
    },
    'solver': {
        'measure_type': 'entan_ln',
        'idx_e': (1, 2)
    },
    'system': {},
    'plotter': {
        'type': 'pcolormesh',
        'palette': 'Greens',
        'bins': 11,
        'title': '$G = g = 3 \\kappa$',
        'x_label': '$\\Delta / \\omega_{m}$',
        'x_ticks': [0.4, 1.0, 1.6],
        'y_label': '$\\omega_{LC} / \\omega_{m}$',
        'y_ticks': [0.4, 1.0, 1.6],
        'cbar_title': '$E_{N}$',
        'cbar_ticks': [0, 0.05, 0.1]
    }
}

# wrapper
wrap_looper(SystemClass=NewJPhys22_063041, params=params, func='measure_stationary', looper='xy_looper', plot=True)