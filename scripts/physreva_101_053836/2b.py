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
from systems import PhysRevA_101_053836

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : np.pi * 50.0,
        't_dim'         : 1001,
        'indices'       : [(3, 3)]
    },
    'system'    : {
        'Delta_a_norm'  : 1.0,
        'G_norms'       : [0.1, 0.01, 0.05],
        'gamma_m_norm'  : 1e-6,
        'kappa_norm'    : 0.1,
        'ns'            : [0.0, 10.0],
        'Omega_norm'    : 2.0,
        't_rwa'         : True
    },
    'plotter'   : {
        'type'          : 'lines',
        'colors'        : ['b', 'r'],
        'x_label'       : '$t / \\tau$',
        'x_ticks'       : [10 * np.pi * i for i in range(6)],
        'x_tick_labels' : [10 * i for i in range(6)],
        'v_label'       : '$\\langle \\delta Y_{b}^{2} \\rangle$',
        'v_ticks'       : [0.0, 1.0, 2.0, 3.0, 4.0],
        'width'         : 8.0,
        'height'        : 4.0
    }
}

# initialize logger
init_log()

# update parameters without RWA
params['system']['t_rwa'] = False
# initialize system
system = PhysRevA_101_053836(
    params=params['system']
)

# initialize solver
solver = HLESolver(
    system=system,
    params=params['solver']
)
# get times and variances
T   = solver.get_times()
M_0 = solver.get_corr_indices().transpose()[0]

# update parameters with RWA
params['system']['t_rwa'] = True
# initialize system
system = PhysRevA_101_053836(
    params=params['system']
)

# get variances
M_1 = HLESolver(
    system=system,
    params=params['solver']
).get_corr_indices().transpose()[0]

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    xs=T,
    vs=[M_0, M_1]
)
plotter.show()