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
from systems import NewJPhys_22_063041

# all parameters
params = {
    'looper'    : {
        'show_progress' : True,
        'X'             : {
            'var'   : 'Delta_norm',
            'min'   : 0.4,
            'max'   : 1.6,
            'dim'   : 121
        },
        'Y'             : {
            'var'   : 'omega_LC_norm',
            'min'   : 0.4,
            'max'   : 1.6,
            'dim'   : 121
        }
    },
    'solver'    : {
        'measure_codes' : ['entan_ln'],
        'indices'       : (1, 2)
    },
    'system'    : {
        'Delta_norm'    : 1.0,
        'G_norm'        : 3.0, 
        'g_norm'        : 3.0,
        'gamma_LC_norm' : 1e-5,
        'gamma_m_norm'  : 1e-6,
        'kappa_norm'    : 0.1,
        'omega_LC_norm' : 1.0,
        'omega_m'       : 2.0 * np.pi * 1e6,
        'T_LC'          : 1e-2,
        'T_m'           : 1e-2,
        't_Delta'       : 'absolute'
    },
    'plotter'   : {
        'type'          : 'pcolormesh',
        'palette'       : 'Greens',
        'title'         : '$G = g = 3 \\kappa$',
        'x_label'       : '$\\Delta / \\omega_{m}$',
        'x_ticks'       : [0.4, 1.0, 1.6],
        'y_label'       : '$\\omega_{LC} / \\omega_{m}$',
        'y_ticks'       : [0.4, 1.0, 1.6],
        'show_cbar'     : True,
        'cbar_title'    : '$E_{N}$',
        'cbar_ticks'    : [0, 0.05, 0.1]
    }
}

# function to obtain the steady state entanglement
def func(system_params):
    return get_func_quantum_correlation_measures(
        SystemClass=NewJPhys_22_063041,
        params=params['solver'],
        steady_state=True
    )(system_params)[0, 0]

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