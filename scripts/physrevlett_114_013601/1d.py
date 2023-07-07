# dependencies
import numpy as np
import os 
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_114_013601

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        'indices'       : [1],
        't_min'         : 0.0,
        't_max'         : 200.0,
        't_dim'         : 2001
    },
    'system'    : {
        'Delta_norm'    : -0.75, 
        'Gamma_norm'    : 1e-3,
        'kappa_norm'    : 1.0,
        'P'             : 1.3
    },
    'plotter'   : {
        'type'          : 'line',
        'x_label'       : '$2 \\pi \\tau$',
        'x_ticks'       : [10 * np.pi * i for i in range(5)],
        'x_tick_labels' : [5 * i for i in range(5)],
        'v_label'       : '$x$',
        'v_ticks'       : [-2, -1, 0, 1],
        'height'        : 4.0
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett_114_013601(
    params=params['system']
)

# initialize solver
solver = HLESolver(
    system=system,
    params=params['solver']
)
# get times and modes
xs  = solver.get_times()
vs  = 2.0 * np.real(solver.get_mode_indices()[:, 0]) / np.sqrt(2.0)

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    xs=xs,
    vs=vs
)
plotter.show()