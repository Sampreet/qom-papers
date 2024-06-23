# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import wrap_looper
from qom.utils.solvers import get_func_quantum_correlation_measures

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_063846_01

# all parameters
params = {
    'looper'    : {
        'show_progress' : True,
        'X'             : {
            'var': 'J_norm',
            'val': [0.5 / np.sqrt(2), 1.0 / np.sqrt(2), 1.0]
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 0.4,
        't_dim'         : 4001
    },
    'system'    : {
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    },
    'plotter'   : {
        'type'          : 'lines',
        'colors'        : ['m', 'r', 'b'],
        'x_label'       : '$\\Gamma t$',
        'x_ticks'       : [0.0, 0.2, 0.4],
        'x_ticks_minor' : [0.04 * i for i in range(11)],
        'v_label'       : '$S$',
        'v_ticks'       : [0, 1, 2],
        'v_ticks_minor' : [0.2 * i for i in range(11)],
        'show_legend'   : True,
        'legend_labels' : [
            '$J = \\Gamma / 2 \\sqrt{2}$',
            '$J = \\Gamma / \\sqrt{2}$',
            '$J = \\Gamma$'
        ],
        'width'         : 4.8,
        'height'        : 4.8
    }
}

# function to obtain the correlation matrix
def func(system_params):
    # initialize system
    system = PhysRevA_100_063846_01(
        params=system_params
    )
    # initialize solver
    solver = HLESolver(
        system=system,
        params=params['solver']
    )
    # get correlation matrices
    corrs = solver.get_corrs()

    return corrs

# looper
looper = wrap_looper(
    looper_name='XLooper',
    func=func,
    params=params['looper'],
    params_system=params['system']
)
vs = [[v[0][0] + v[2][2] / 2.0 + v[4][4] / 2.0 - np.sqrt(2) * v[0][2] - np.sqrt(2) * v[0][4] + v[2][4] + v[1][1] + v[3][3] / 2.0 + v[5][5] / 2.0 + np.sqrt(2) * v[1][3] + np.sqrt(2) * v[1][5] + v[3][5] for v in row] for row in looper.results['V']]
xs = np.linspace(params['solver']['t_min'], params['solver']['t_max'], params['solver']['t_dim'])
ys = looper.results['X']

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    vs=vs,
    xs=xs,
    ys=ys
)
plotter.show()