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
        'file_path_prefix'  : 'data/physrevlett_111_103605/2c',
        'X'                 : {
            'var'   : 'n_b',
            'min'   : 0.0,
            'max'   : 20.0,
            'dim'   : 11
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        'measure_codes' : ['sync_c', 'sync_p'],
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
        'colors'        : ['b', 'g'],
        'sizes'         : [50, 50],
        'styles'        : ['o', 's'],
        'x_label'       : '$n_{b}$',
        'x_ticks'       : [0, 5, 10, 15, 20],
        'x_ticks_minor' : [i for i in range(21)],
        'v_label'       : '$\\bar{S}_{c}, \\bar{S}_{p}$',
        'v_ticks'       : [0, 0.1, 0.2, 0.3],
        'v_ticks_minor' : [i * 0.02 for i in range(16)],
        'show_legend'   : True,
        'legend_labels' : ['$\\bar{S}_{c}$', '$\\bar{S}_{p}$'],
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
        vs=np.transpose(looper.results['V'])
    )
    plotter.show()