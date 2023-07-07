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
from systems import PhysRevLett_103_213603

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.00,
        't_max'         : 50.00,
        't_dim'         : 5001,
        't_index_min'   : 0,
        't_index_max'   : 5001,
        'indices'       : [1]
    },
    'system'    : {
        'F'         : 1.4e4,
        'L'         : 25e-3,
        'lambda_l'  : 1064e-9,
        'm'         : 150e-12,
        'omega_m'   : 2.0 * np.pi * 1e6,
        'Ps'        : [10e-3, 2e-3],
        'Q'         : 1e6,
        'T'         : 0.1
    },
    'plotter'   : {
        'type'          : 'line',
        'x_label'       : '$\\langle q \\rangle / 10^{3}$',
        'x_limits'      : [10e3, 30e3],
        'x_tick_labels' : [10, 20, 30],
        'x_ticks'       : [10e3, 20e3, 30e3],
        'x_ticks_minor' : [i * 2e3 + 10e3 for i in range(10)],
        'y_colors'      : ['k', 'g'],
        'v_label'       : '$\\langle p \\rangle / 10^{3}$',
        'v_limits'      : [-12e3, 10e3],
        'v_tick_labels' : [-10, 0, 10],
        'v_ticks'       : [-10e3, 0, 10e3],
        'v_ticks_minor' : [i * 2e3 - 12e3 for i in range(12)],
        'width'         : 4.2,
        'height'        : 4.0,
        'annotations'   : [{
            'text'  : '(a)',
            'xy'    : [0.29, 0.84]
        }]
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett_103_213603(
    params=params['system']
)

# get modes
Modes = HLESolver(
    system=system,
    params=params['solver']
).get_mode_indices().transpose()[0]
qs = np.real(Modes) * np.sqrt(2.0)
ps = np.imag(Modes) * np.sqrt(2.0)

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    xs=[qs, qs[-628:]],
    vs=[ps, ps[-628:]]
)
plotter.show()