# dependencies
import numpy as np
import os 
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper
from qom.utils.solvers import get_func_Lyapunov_exponents

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_114_013601

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physrevlett_114_013601/2b',
        'X'                 : {
            'var'   : 'Delta_norm',
            'min'   : -1.2,
            'max'   : -0.4,
            'dim'   : 801
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : True,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 10000.0,
        't_dim'         : 100001,
        'num_steps'     : 100000,
        'step_size'     : 0.1,
        'use_svd'       : True
    },
    'system'    : {
        'Delta_norm'    : 0.0, 
        'Gamma_norm'    : 1e-3,
        'kappa_norm'    : 1.0,
        'P'             : 1.4
    },
    'plotter'   : {
        'type'          : 'line',
        'x_label'       : '$\\Delta / \\Omega$',
        'x_ticks'       : [-1.2, -0.8, -0.4],
        'x_ticks_minor' : [i * 0.05 - 1.2 for i in range(17)],
        'v_label'       : '$\\lambda_{max}$',
        'v_limits'      : [-0.25, 0.05],
        'v_ticks'       : [-0.2, -0.1, 0.0],
        'v_ticks_minor' : [i * 0.05 -0.2 for i in range(9)],
        'width'         : 6.0,
        'height'        : 3.0
    }
}

# wrapper function for parallel processes
def func(system_params):
    return get_func_Lyapunov_exponents(
        SystemClass=PhysRevLett_114_013601,
        params=params['solver'],
        steady_state=False
    )(system_params)

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
        vs=np.max(looper.results['V'], axis=1)
    )
    plotter.show()