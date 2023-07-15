# dependencies
import os 
import sys

# qom modules
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper
from qom.utils.solvers import get_func_stability_zone

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import NewJPhys_22_013049

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/newjphys_22_013049/2b',
        'X'                 : {
            'var'   : 'Delta_norm',
            'min'   : -1.0,
            'max'   : 1.0,
            'dim'   : 1001
        },
        'Y'                 : {
            'var'   : 'P',
            'min'   : 1e-4,
            'max'   : 1e0,
            'dim'   : 1001,
            'scale' : 'log'
        }
    },
    'solver'    : {
        'show_progress' : False,
        'method'        : 'eig'
    },
    'system'    : {
        'Delta_norm'    : 0.0, 
        'gamma_norm'    : 1e-4,
        'kappa_norm'    : 0.1,
        'P'             : 0.1
    },
    'plotter'   : {
        'type'              : 'pcolormesh',
        'palette'           : 'Blues',
        'x_label'           : '$\\Delta / \\Omega_{m}$',
        'x_ticks'           : [-1.0, -0.5, 0.0, 0.5, 1.0],
        'y_label'           : '$P$',
        'y_scale'           : 'log',
        'y_ticks'           : [0.0001, 0.001, 0.01, 0.1, 1.0],
        'y_tick_labels'     : ['$10^{' + str(i - 4) + '}$' for i in range(5)],
        'show_cbar'         : True,
        'cbar_position'     : 'top',
        'cbar_title'        : 'Stability Zone of $\\left| \\alpha \\right|^{2}$',
        'cbar_tick_labels'  : ['0S1U', '1S0U', '0S3U', '1S2U', '2S1U', '3S0U'],
        'cbar_ticks'        : list(range(6))
    }
}

# wrapper function for parallel processes
def func(system_params):
    return get_func_stability_zone(
        SystemClass=NewJPhys_22_013049,
        params=params['solver'],
        steady_state=True,
        use_rhc=False
    )(system_params)

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