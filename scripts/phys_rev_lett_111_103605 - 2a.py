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
from systems import PhysRevLett111_103605

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        'idx_e': (1, 3),
        'range_min': 0,
        'range_max': 1001,
        't_min': 0.0,
        't_max': 100.0,
        't_dim': 1001
    },
    'system': {
        'mu': 0.02
    },
    'plotter': {
        'type': 'lines',
        'x_label': '$t / \\tau$',
        'x_ticks': [0, 20, 40, 60, 80, 100],
        'v_label': '$S_{c}, S_{p}$',
        'y_colors': ['b', 'g'],
        'v_ticks': [0.0, 0.1, 0.2, 0.3],
        'height': 4.0,
        'width': 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett111_103605(params=params['system'])

# get complete synchronization
params['solver']['measure_type'] = 'sync_c'
M_0, T = system.get_measure_dynamics(solver_params=params['solver'])
sync_c_avg = np.mean(M_0[371:])
# get phase synchronization
params['solver']['measure_type'] = 'sync_p'
M_1, T = system.get_measure_dynamics(solver_params=params['solver'])
sync_p_avg = np.mean(M_1[371:])

# plotter
plotter = MPLPlotter(axes={
    'X': T,
    'Y': {
        'var': 'M',
        'val': [0, 2]
    }
}, params=params['plotter'])
axis = plotter.get_current_axis()
axis.plot([sync_c_avg for i in range(len(T))], linestyle='--', color='b')
axis.plot([sync_p_avg for i in range(len(T))], linestyle='--', color='g')
_xs = [T, T]
_vs = [M_0, M_1]
plotter.update(xs=_xs, vs=_vs)
plotter.show(True)