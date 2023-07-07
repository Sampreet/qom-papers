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
        'file_path_prefix'  : 'data/physrevlett_114_013601/4',
        'X'                 : {
            'var'   : 'Delta_norm',
            'min'   : -1.5,
            'max'   : 0.0,
            'dim'   : 151
        },
        'Y': {
            'var'   : 'P',
            'min'   : 1.0,
            'max'   : 1.6,
            'dim'   : 151
        }
    },
    'solver': {
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
    'system': {
        'Delta_norm'    : 0.0, 
        'Gamma_norm'    : 1e-3, 
        'kappa_norm'    : 1.0,
        'P'             : 1.0
    },
    'plotter': {
        'type'          : 'pcolormesh',
        'x_label'       : '$\\Delta / \\Omega$',
        'x_ticks'       : [-1.5, -1.0, -0.5, 0.0],
        'y_label'       : '$P$',
        'y_ticks'       : [1.0, 1.2, 1.4, 1.6],
        'cbar_title'    : '$\\lambda_{max}$',
        'cbar_ticks'    : [-0.005, 0.00, 0.005],
        'show_cbar'     : True,
        'cbar_position' : 'top'
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
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system']
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