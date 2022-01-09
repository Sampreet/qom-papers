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
        't_min': 0.00,
        't_max': 50.00,
        't_dim': 5001
    },
    'system': {
        'F': 1.4e4,
        'L': 25e-3,
        'lambda_l': 1064e-9,
        'm': 150e-12,
        'omega_m': 2e6 * np.pi,
        'P_0': 10e-3,
        'P_1': 2e-3,
        'Q': 1e6,
        'T': 0.1
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
system = PhysRevLett103_213603(params=params['system'])

# get measure dynamics
M, T = system.get_measure_dynamics(solver_params=params['solver'])
qs = np.real(np.transpose(M)[0]).tolist()
ps = np.imag(np.transpose(M)[0]).tolist()

# plotter
plotter = MPLPlotter(axes={
    'X': qs
}, params=params['plotter'])
plotter.update(xs=qs, vs=ps)
plotter.show(True)