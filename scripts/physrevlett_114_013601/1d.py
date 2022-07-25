# dependencies
import numpy as np
import os 
import sys

# qom modules
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevLett_114_013601

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': False,
        'method': 'zvode',
        'measure_type': 'mode_amp',
        'idx_e': 1,
        't_min': 0.0,
        't_max': 200.0,
        't_dim': 2001
    },
    'system': {
        'Delta_norm': -0.75, 
        'Gamma_norm': 1e-3,
        'kappa_norm': 1.0,
        'P': 1.3
    },
    'plotter': {
        'type': 'line',
        'x_label': '$2 \\pi \\tau$',
        'x_ticks': [10 * np.pi * i for i in range(5)],
        'x_tick_labels': [5 * i for i in range(5)],
        'v_label': '$x$',
        'v_ticks': [-3, -2, -1, 0, 1, 2]
    }
}

# initilize log
init_log()

# initialize system
system = PhysRevLett_114_013601(params=params['system'])

# get dynamics
M, T = system.get_measure_dynamics(solver_params=params['solver'])
# extract positions
vs = [2 * np.real(m) for m in np.transpose(M)[0]]

# plotter
plotter = MPLPlotter(axes={
    'X': T
}, params=params['plotter'])
plotter.update(xs=T, vs=vs)
plotter.show(True)