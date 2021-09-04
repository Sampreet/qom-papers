# dependencies
import os
import sys

# qom modules
from qom.ui import init_log

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevLett119_153901

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        't_min': 0.0,
        't_max': 4.0,
        't_dim': 4001
    },
    'system': {
        't_alphas': 'Gaussian'
    },
    'plotter': {
        'type': 'pcolormesh',
        'x_label': '$x$',
        'x_ticks': [0, 200, 400],
        'x_tick_labels': [-1, 0, 1],
        'y_label': '$\\tau$',
        'y_ticks': [0, 2, 4],
        'show_cbar': True,
        'cbar_title': '$|\\alpha| (10^{3})$',
        'cbar_ticks': [0, 4000, 8000],
        'cbar_tick_labels': [0, 4, 8]
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett119_153901(params=params['system'])
# get mode amplitude dynamics
amps, T, X = system.get_mode_amplitude_dynamics(solver_params=params['solver'], plot=True, plotter_params=params['plotter'])