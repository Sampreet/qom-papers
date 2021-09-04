# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui import init_log

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import OptLett41_2676

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        't_min': 0.0,
        't_max': 4.0,
        't_dim': 401
    },
    'system': {
        'n': 301,
        'x_0': 20.0,
        'order': 2
    },
    'plotter': {
        'type': 'surface_cz',
        'x_label': '$x$',
        'x_ticks': [0, 150, 300],
        'y_label': '$\\tau$',
        'y_ticks': [0, np.pi / 2, np.pi, 4],
        'y_tick_labels': ['0', '$\\pi / 2$', '$\\pi$', ''],
        'v_label': '$|\\alpha| (10^{3})$',
        'v_ticks': [0, 1000, 2000],
        'v_tick_labels': [0, 1, 2],
        'show_cbar': False,
        'width': 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = OptLett41_2676(params=params['system'])

# get mode amplitude dynamics
amps, T, X = system.get_mode_amplitude_dynamics(solver_params=params['solver'], plot=True, plotter_params=params['plotter'])