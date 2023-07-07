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
from systems import NewJPhys_22_013049

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/newjphys_22_013049/5',
        'X'                 : {
            'var'   : 'Delta_norm',
            'min'   : -1.0,
            'max'   : 1.0,
            'dim'   : 201
        },
        'Y'                 : {
            'var'   : 'P',
            'min'   : 0.0,
            'max'   : 0.5,
            'dim'   : 501
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : True,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 1000.0,
        't_dim'         : 10001,
        'num_steps'     : 10000,
        'step_size'     : 0.1,
        'use_svd'       : True
    },
    'system'    : {
        'Delta_norm'    : 0.0, 
        'gamma_norm'    : 1e-4,
        'kappa_norm'    : 0.1,
        'P'             : 0.1
    },
    'plotter'   : {
        'type'          : 'pcolormesh',
        'x_label'       : '$\\Delta / \\Omega_{m}$',
        'x_ticks'       : [-1.0, 0.0, 1.0],
        'y_label'       : '$P$',
        'y_ticks'       : [0.00, 0.25, 0.50],
        'show_cbar'     : True,
        'cbar_title'    : '$\\lambda_{max}$',
        'cbar_ticks'    : [-0.005, 0, 0.005],
        'cbar_position' : 'top'
    }
}

# wrapper function for parallel processes
def func(system_params):
    return get_func_Lyapunov_exponents(
        SystemClass=NewJPhys_22_013049,
        params=params['solver'],
        steady_state=False
    )(system_params)

if __name__ == '__main__':
    # wrap looper
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=False
    )

    # plotter
    plotter = MPLPlotter(
        axes={
            'X' : looper.axes['X']['val'],
            'Y' : looper.axes['Y']['val']
        },
        params=params['plotter']
    )
    plotter.update(
        vs=np.max(looper.results['V'], axis=2)
    )
    plotter.show()