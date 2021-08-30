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
        'measure_type': 'entan_ln',
        'idx_e': (0, 1),
        'range_min': 3000,
        'range_max': 3201,
        't_min': 0,
        't_max': 50,
        't_dim': 5001
    },
    'system': {},
    'plotter': {
        'type': 'lines',
        'show_legend': True,
        'x_label': '$t / \\tau$',
        'y_name': '$P_{\\pm 1}$',
        'y_unit': '$\\mathrm{mW}$',
        'y_colors': ['b', 'g'],
        'v_label': '$E_{N}$',
        'v_ticks': [0.0, 0.1, 0.2, 0.3, 0.4],
        'width': 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett103_213603(params['system'])

# get measure dynamics without modulation
system.params['P_1'] = 0.0
M_0, T = system.get_measure_dynamics(params['solver'])

# get measure dynamics with modulation
system.params['P_1'] = 2e-3
M_1, T = system.get_measure_dynamics(params['solver'])

# plotter
plotter = MPLPlotter({
    'X': T,
    'Y': [0, 2]
}, params['plotter'])
plotter.update(xs=[T, T], vs=[M_0, M_1])
plotter.show(True)