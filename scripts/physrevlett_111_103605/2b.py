# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper
from qom.utils.solvers import get_func_quantum_correlation_measures

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_111_103605

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physrevlett_111_103605/2b',
        'X'                 : {
            'var'   : 'mu_norm',
            'min'   : 0.00,
            'max'   : 0.04,
            'dim'   : 11
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        'measure_codes' : ['sync_c', 'sync_p', 'discord_G'],
        'indices'       : (1, 3),
        't_min'         : 0.0,
        't_max'         : 200.0,
        't_dim'         : 2001,
        't_index_min'   : 1371,
        't_index_max'   : 2001
    },
    'system'    : {
        'E_norm'        : 320.0,
        'g_norm'        : 0.005,
        'gamma_norm'    : 0.005,
        'kappa_norm'    : 0.15,
        'mu_norm'       : 0.02,
        'n_b'           : 0.0,
        'omega_2_norm'  : 1.0
    },
    'plotter'   : {
        'type'          : 'scatters',
        'x_label'       : '$\\mu$',
        'x_ticks'       : [0.00, 0.01, 0.02, 0.03, 0.04],
        'x_ticks_minor' : [i * 0.002 for i in range(21)],
        'y_colors'      : ['b', 'g', 'y'],
        'y_legend'      : ['$\\bar{S}_{c}$', '$\\bar{S}_{p}$', '$\\bar{D}_{G}$'],
        'y_sizes'       : [50, 50, 50],
        'y_styles'      : ['o', 's', 'D'],
        'v_label'       : '$\\bar{S}_{c}, \\bar{S}_{p}, \\bar{D}_{G}$',
        'v_ticks'       : [0, 0.1, 0.2, 0.3],
        'v_ticks_minor' : [i * 0.02 for i in range(16)],
        'show_legend'   : True,
        'width'         : 6.0,
        'height'        : 4.0
    }
}

# function to obtain the average quantum correlation measures
def func(system_params):
    Measures = get_func_quantum_correlation_measures(
        SystemClass=PhysRevLett_111_103605,
        params=params['solver'],
        steady_state=False
    )(system_params)
    return np.mean(Measures, axis=0)

if __name__ == '__main__':
    # wrap looper
    looper = run_loopers_in_parallel(
        looper_name='XLooper',
        func=func,
        params=params['looper'],
        params_system=params['system']
    )
    # plotter
    plotter = MPLPlotter(
        axes={},
        params=params['plotter']
    )
    plotter.update(
        xs=looper.axes['X']['val'],
        vs=looper.results['V'].transpose()
    )
    plotter.show()