# dependencies
import os 
import sys

# qom modules
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import NewJPhys_22_013049

# all parameters
params = {
    'looper': {
        'X': {
            'var': 'Delta_norm',
            'min': -1.0,
            'max': 1.0,
            'dim': 201
        },
        'Y': {
            'var': 'P',
            'min': -4,
            'max': 0,
            'dim': 201,
            'scale': 'log'
        }
    },
    'system': {
        'Delta_norm': 0.0, 
        'gamma_norm': 1e-4,
        'kappa_norm': 0.1,
        'P': 0.1
    },
    'plotter': {
        'type': 'pcolormesh',
        'palette': 'Blues',
        'x_label': '$\\Delta / \\Omega_{m}$',
        'x_ticks': [-1.0, -0.5, 0.0, 0.5, 1.0],
        'y_label': '$P$',
        'y_scale': 'log',
        'y_ticks': [0.0001, 0.001, 0.01, 0.1, 1.0],
        'y_tick_labels': ['$10^{' + str(i - 4) + '}$' for i in range(5)],
        'show_cbar': True,
        'cbar_title': 'Stability Zone of $\\left| \\alpha \\right|^{2}$',
        'cbar_position': 'top',
        'cbar_ticks': list(range(6)),
        'cbar_tick_labels': ['0S1U', '1S0U', '0S3U', '1S2U', '2S1U', '3S0U']
    }
}

# wrapper
looper = wrap_looper(SystemClass=NewJPhys_22_013049, params=params, func='osz', looper='XYLooper', file_path_prefix='data/newjphys_22_013049/2b_osz', plot=True)