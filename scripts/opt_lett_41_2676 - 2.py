# dependencies
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
        't_max': 2.0,
        't_dim': 201
    },
    'system': {
        'n': 151,
        'x_0': 10.0
    },
    'plotter': {
        'type': 'surface_cz',
        'x_label': '$x$',
        'x_bound': 'both',
        'x_ticks': [0, 75, 150],
        'y_label': '$\\tau$',
        'y_bound': 'both',
        'y_ticks': [0, 1, 2],
        'v_label': '$|\\alpha| (10^{3})$',
        'v_bound': 'both',
        'v_ticks': [0, 1000, 2000],
        'v_tick_labels': [0, 1, 2],
        'show_cbar': False,
        'width': 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = OptLett41_2676(params['system'])

# get mode amplitude dynamics
amps, T, X = system.get_mode_amplitude_dynamics(solver_params=params['solver'], plot=True, plotter_params=params['plotter'])