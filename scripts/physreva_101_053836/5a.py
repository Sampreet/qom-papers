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
from systems import PhysRevA_101_053836

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physreva_101_053836/5',
        'X'             : {
            'var'   : 'G_p1_norm',
            'min'   : 0.0,
            'max'   : 1.0,
            'dim'   : 1001
        },
        'Y'             : {
            'var'   : 'ns',
            'idx'   : 1,
            'val'   : [10.0, 100.0]
        },
        'Z'             : {
            'var'   : 'G_0_norm',
            'min'   : 0.0,
            'max'   : 2.0,
            'dim'   : 101
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
        'x_label'       : '$G_{0} / \\kappa$',
        'x_ticks'       : [0.2 + 0.2 * i for i in range(10)],
        'x_tick_labels' : ['{:0.1f}'.format(0.2 + 0.2 * i) for i in range(10)],
        'y_colors'      : ['b', 'r'],
        'y_legend'      : ['$n_{m} = 10$', '$n_{m} = 100$'],
        'v_label'       : '$\\langle \\delta \\tilde{X}_{b}^{2} \\rangle_{max}$ (dB)',
        'v_ticks'       : [0, 5, 10, 15, 20, 25],
        'show_legend'   : True,
        'width'         : 8.0,
        'height'        : 4.0
    }
}

# function to calculate the average variance in decibels
def func(system_params):
    # update parameter
    system_params['G_norms'][0] = system_params['G_0_norm'] * system_params['kappa_norm']
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
    # wrap looper
    looper = run_loopers_in_parallel(
        looper_name='XYZLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )
    # extract values
    xs = looper.axes['Z']['val']
    vs = np.max(looper.results['V'], axis=2).transpose()

    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=xs,
        vs=vs
    )
    plotter.show()