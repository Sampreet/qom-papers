# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper
from qom.utils.solvers import get_func_quantum_correlation_measures

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_063846_00

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physreva_100_063846_00/8a',
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
        't_max'         : 0.5,
        't_dim'         : 101,
        'measure_codes' : ['entan_ln_2'],
        'indices'       : (0, 1)
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
        'cbar_title'    : '$E_{N}$',
        'cbar_ticks'    : [0.0, 0.2, 0.4],
        'width'         : 5.5,
        'height'        : 4.8
    }
}

# function to obtain the entanglement
def func(system_params):
    return get_func_quantum_correlation_measures(
        SystemClass=PhysRevA_100_063846_00,
        params=params['solver'],
        steady_state=False
    )(system_params)[-1, 0]

if __name__ == '__main__':
    # looper
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=True,
        params_plotter=params['plotter']
    )