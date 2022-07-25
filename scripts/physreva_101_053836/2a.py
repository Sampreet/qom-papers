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
from systems import PhysRevA_101_053836

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': False,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (2, 2),
        'range_min': 0,
        'range_max': 1001,
        't_min': 0.0,
        't_max': np.pi * 50.0,
        't_dim': 1001
    },
    'system': {
        'Delta_a_norm': 1.0,
        'G_norms': [0.1, 0.01, 0.05],
        'gamma_m_norm': 1e-6,
        'kappa_norm': 0.1,
        'ns': [0.0, 10.0],
        'Omega_norm': 2.0,
        't_rwa': True
    },
    'plotter': {
        'type': 'lines',
        'x_label': '$t / \\tau$',
        'x_ticks': [10 * np.pi * i for i in range(6)],
        'x_tick_labels': [10 * i for i in range(6)],
        'y_colors': ['b', 'r'],
        'v_label': '$\\langle \\delta X_{b}^{2} \\rangle$',
        'v_ticks': [0.0, 1.0, 2.0, 3.0, 4.0],
        'show_legend': True,
        'width': 8.0,
        'height': 4.0
    }
}

# initialize logger
init_log()

# initialize system without RWA
params['system']['t_rwa'] = False
system = PhysRevA_101_053836(params=params['system'])
# get measure dynamics
M_0, T = system.get_measure_dynamics(solver_params=params['solver'])

# initialize system with RWA
params['system']['t_rwa'] = True
system = PhysRevA_101_053836(params=params['system'])
# get measure dynamics
M_1, T = system.get_measure_dynamics(solver_params=params['solver'])

# plotter
plotter = MPLPlotter(axes={
    'X': T,
    'Y': ['without RWA', 'with RWA']
}, params=params['plotter'])
plotter.update(xs=T, vs=[np.transpose(M_0)[0], np.transpose(M_1)[0]])
plotter.show(True)