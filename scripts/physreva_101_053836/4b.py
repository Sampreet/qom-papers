# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevA_101_053836

# all parameters
params = {
    'looper': {
        'show_progress_x': True,
        'X': {
            'var': 'G_p1_norm',
            'min': 0.0,
            'max': 1.0,
            'dim': 1001
        },
        'Y': {
            'var': 'ns',
            'idx': 1,
            'val': [10.0, 100.0]
        }
    },
    'solver': {
        'cache': False,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': [(2, 2), (3, 3)],
        'range_min': 990,
        'range_max': 1001,
        't_min': 0.0,
        't_max': np.pi * 100.0,
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
        'x_label': '$G_{1} / G_{0}$',
        'x_ticks': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        'y_colors': ['b', 'r'],
        'y_legend': ['$n_{m} = 10$', '$n_{m} = 100$'],
        'v_label': '$\\langle \\beta^{\\dagger} \\beta \\rangle$',
        'v_ticks': [1e-4, 1, 1e4],
        'v_tick_labels': ['$10^{-4}$', '$10^{0}$', '$10^{4}$'],
        'v_scale': 'log',
        'show_legend': True,
        'width': 8.0,
        'height': 4.0
    }
}

# function to calculate the phonon number in the Bogoluibov mode
def func_n_beta(system_params, val, logger, results):
    # update parameter
    system_params['G_norms'][2] = system_params['G_p1_norm'] * system_params['G_norms'][0]
    # initialize system
    system = PhysRevA_101_053836(params=system_params)
    # get dynamics of correlation elements
    M, _ = system.get_measure_dynamics(solver_params=params['solver'])
    # calculate hyperbolic angles
    r = np.arctanh(system_params['G_p1_norm'])
    chr = np.cosh(r)
    shr = np.sinh(r)
    # get phonon number in the Bogoluibov mode
    n_beta = np.mean([(chr**2 + shr**2) * (v[0] + v[1] - 1) / 2 + shr**2 + chr * shr * (v[0] - v[1]) for v in M])
    # update results
    results.append((val, n_beta))

# wrapper
looper = wrap_looper(SystemClass=None, params=params, func=func_n_beta, looper='XYLooper', file_path_prefix='data/physreva_101_053836/4b_n_beta', plot=True)