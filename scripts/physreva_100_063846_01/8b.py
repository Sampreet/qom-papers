# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_063846_01

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physreva_100_063846_01/8b',
        'X'                 : {
            'var': 'J_norm',
            'val': np.linspace(0.5, 1.0, 501)
        },
        'Y'                 : {
            'var': 'n_th',
            'val': np.linspace(0.0, 300.0, 301)
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 0.1,
        't_dim'         : 101
    },
    'system'    : {
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    },
    'plotter'   : {
        'type'          : 'pcolormesh',
        'bins'          : 41,
        'x_label'       : '$J / \\Gamma$',
        'x_ticks'       : [0.5, 0.75, 1.0],
        'y_label'       : '$n_{th}$',
        'y_ticks'       : [0, 150, 300],
        'show_cbar'     : True,
        'cbar_title'    : '$S$',
        'cbar_ticks'    : [0.9, 1.0, 1.1],
        'width'         : 5.5,
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
    # get correlation matrix at final time
    corrs = solver.get_corrs()[-1]

    return corrs

if __name__ == '__main__':
    # looper
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system']
    )
    # tripartite entanglement measure
    vs = [[v[0][0] + v[2][2] / 2.0 + v[4][4] / 2.0 - np.sqrt(2) * v[0][2] - np.sqrt(2) * v[0][4] + v[2][4] + v[1][1] + v[3][3] / 2.0 + v[5][5] / 2.0 + np.sqrt(2) * v[1][3] + np.sqrt(2) * v[1][5] + v[3][5] for v in row] for row in looper.results['V']]
    # plotter
    plotter = MPLPlotter(
        axes={
            'X' : looper.axes['X']['val'],
            'Y' : looper.axes['Y']['val']
        },
        params=params['plotter']
    )
    plotter.update(
        vs=vs
    )
    plotter.show()
    