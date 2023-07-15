# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_101_053836

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physreva_101_053836/4a',
        'X'                 : {
            'var'   : 'G_p1_norm',
            'min'   : 0.0,
            'max'   : 1.0,
            'dim'   : 1001
        },
        'Y'             : {
            'var'   : 'ns',
            'idx'   : 1,
            'val'   : [10.0, 100.0]
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 2 * np.pi * 100.0,
        't_dim'         : 10001,
        't_index_min'   : 9000,
        't_index_max'   : 10001,
        'indices'       : [(2, 2)]
    },
    'system'    : {
        'Delta_a_norm'  : 1.0,
        'G_norms'       : [0.1, 0.01, 0.05],
        'gamma_m_norm'  : 1e-6,
        'kappa_norm'    : 0.1,
        'ns'            : [0.0, 10.0],
        'Omega_norm'    : 2.0,
        't_rwa'         : True
    },
    'plotter'   : {
        'type'          : 'lines',
        'colors'        : ['b', 'r'],
        'x_label'       : '$G_{1} / G_{0}$',
        'x_ticks'       : [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        'v_label'       : '$\\langle \\delta \\tilde{X}_{b}^{2} \\rangle$',
        'v_ticks'       : [0.0, 5.0, 10.0, 15.0, 20.0],
        'show_legend'   : True,
        'legend_labels' : ['$n_{m} = 10$', '$n_{m} = 100$'],
        'width'         : 8.0,
        'height'        : 4.0
    }
}

# function to calculate the average variance in decibels
def func(system_params):
    # update parameter
    system_params['G_norms'][2] = system_params['G_p1_norm'] * system_params['G_norms'][0]
    # initialize system
    system = PhysRevA_101_053836(
        params=system_params
    )

    # get average correlation element
    corrs = np.mean(HLESolver(
        system=system,
        params=params['solver']
    ).get_corr_indices(), axis=0)

    # convert to decibels and return
    return - 10 * np.log10(corrs[0] / 0.5)

if __name__ == '__main__':
    # wrap looper and plot
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=True,
        params_plotter=params['plotter']
    )