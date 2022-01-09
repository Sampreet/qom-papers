# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevLett111_103605

# all parameters
params = {
    'looper': {
        'show_progress': True,
        'X': {
            'var': 'n_b',
            'min': 0,
            'max': 20,
            'dim': 11
        }
    },
    'solver': {
        'cache': True,
        'method': 'zvode',
        'idx_e': (1, 3),
        'range_min': 371,
        'range_max': 1001,
        't_min': 0.0,
        't_max': 100.0,
        't_dim': 1001
    },
    'system': {
        'E': 320,
        'g': 0.005,
        'gamma': 0.005,
        'kappa': 0.15,
        'mu': 0.02,
        'n_b': 0,
        'omega_1': 1.0,
        'omega_2': 1.0
    },
    'plotter': {
        'type': 'scatters',
        'x_label': '$n_{b}$',
        'x_ticks': [0, 5, 10, 15, 20],
        'y_colors': ['b', 'g'],
        'y_sizes': [25, 25],
        'y_styles': ['o', 's'],
        'v_label': '$\\bar{S}_{c}, \\bar{S}_{p}$',
        'v_ticks': [0, 0.1, 0.2, 0.3],
        'height': 4.0,
        'width': 6.0
    }
}

# get average complete synchronization
params['solver']['measure_type'] = 'sync_c'
looper = wrap_looper(SystemClass=PhysRevLett111_103605, params=params, func='mav', looper='x_looper')
Sync_c_avg = looper.results['V']

# get average phase synchronization
params['solver']['measure_type'] = 'sync_p'
looper = wrap_looper(SystemClass=PhysRevLett111_103605, params=params, func='mav', looper='x_looper')
Sync_p_avg = looper.results['V']

# plotter
X = looper.axes['X']['val']
plotter = MPLPlotter(axes={
    'X': X,
    'Y': ['$\\bar{S}_{c}$', '$\\bar{S}_{p}$']
}, params=params['plotter'])
plotter.update(xs=[X, X], vs=[Sync_c_avg, Sync_p_avg])
plotter.show(True)