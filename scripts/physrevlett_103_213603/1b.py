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
        'indices'       : [0]
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
        'type'          : 'lines',
        'x_label'       : 'Im $\\langle a \\rangle / 10^{4}$',
        'x_tick_labels' : [-4, 2, 8],
        'x_ticks'       : [-40e3, 20e3, 80e3],
        'x_ticks_minor' : [i * 2e4 - 4e4 for i in range(7)],
        'y_colors'      : ['b', 'r'],
        'v_label'       : 'Re $\\langle a \\rangle / 10^{4}$',
        'v_tick_labels' : [-10, -6, -2],
        'v_ticks'       : [-100e3, -60e3, -20e3],
        'v_ticks_minor' : [i * 1e4 - 1e5 for i in range(10)],
        'width'         : 4.2,
        'height'        : 4.0,
        'annotations'   : [{
            'text'  : '(b)',
            'xy'    : [0.31, 0.84]
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
qs = np.real(Modes)
ps = np.imag(Modes)

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