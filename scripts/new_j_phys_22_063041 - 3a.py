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
            'var': 'g_norm',
            'min': 0,
            'max': 8,
            'dim': 81
        },
        'Y': {
            'var': 'G_norm',
            'min': 0,
            'max': 8,
            'dim': 81
        }
    },
    'solver': {
        'measure_type': 'entan_ln',
        'idx_e': (1, 2)
    },
    'system': {
        'G_norm': 5, 
        'g_norm': 5,
        't_Delta': 'relative'
    },
    'plotter': {
        'type': 'pcolormesh',
        'palette': 'Greens',
        'title': '$\\Delta = \\omega_{LC}$',
        'x_label': '$g / \\kappa$',
        'x_ticks': [0, 3, 6],
        'y_label': '$G / \\kappa$',
        'y_ticks': [0, 3, 6],
        'cbar_title': '$E_{N}$',
        'cbar_ticks': [0, 0.1, 0.2, 0.3]
    }
}

# wrapper
wrap_looper(SystemClass=NewJPhys22_063041, params=params, func='measure_stationary', looper='xy_looper', plot=True)