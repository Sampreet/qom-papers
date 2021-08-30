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
from systems import PhysRevLett103_213603

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        'measure_type': 'mode_amp',
        'idx_e': 1,
        'range_min': 0,
        'range_max': 5001,
        't_min': 0,
        't_max': 50,
        't_dim': 5001
    },
    'system': {
        'P_1': 2e-3
    },
    'plotter': {
        'type': 'line',
        'x_label': '$\\langle q \\rangle (10^{3})$',
        'x_ticks': [8e3, 14e3, 20e3],
        'x_tick_labels': [8, 14, 20],
        'v_label': '$\\langle p \\rangle (10^{3})$',
        'v_ticks': [-10e3, -1e3, 8e3],
        'v_tick_labels': [-10, -1, 8],
        'width': 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett103_213603(params['system'])

# get measure dynamics
M, T = system.get_measure_dynamics(params['solver'])
qs = np.real(np.transpose(M)[0]).tolist()
ps = np.imag(np.transpose(M)[0]).tolist()

# plotter
plotter = MPLPlotter({
    'X': qs
}, params['plotter'])
plotter.update(xs=qs, vs=ps)
plotter.show(True)