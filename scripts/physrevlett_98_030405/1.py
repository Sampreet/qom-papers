# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import wrap_looper
from qom.utils.solvers import get_func_quantum_correlation_measures

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_98_030405

# all parameters
params = {
    'looper'    : {
        'show_progress' : True,
        'X'             : {
            'var'   : 'Delta_norm',
            'min'   : 0.0,
            'max'   : 5.0,
            'dim'   : 1001
        }
    },
    'solver'    : {
        'use_system_method' : True,
        'measure_codes'     : ['entan_ln_2'],
        'indices'           : (0, 1)
    },
    'system'    : {
        'Delta_norm'    : 1.0,
        'F'             : 1.07e4,
        'gamma_m'       : 2.0 * np.pi * 1e2,
        'L'             : 1e-3,
        'lamb'          : 810e-9,
        'm'             : 5e-12,
        'omega_m'       : 2.0 * np.pi * 10e6,
        'P'             : 50e-3,
        'T'             : 400e-3
    },
    'plotter'   : {
        'type'          : 'lines',
        'colors'        : ['b', 'r'],
        'styles'        : ['-', '--'],
        'x_label'       : '$\\Delta / w_{m}$',
        'x_ticks'       : [i for i in range(6)],
        'x_ticks_minor' : [i * 0.25 for i in range(21)],
        'v_label'       : '$E_{N}$',
        'v_limits'      : [0.0, 0.32],
        'v_ticks'       : [0.0, 0.1, 0.2, 0.3],
        'v_ticks_minor' : [i * 0.02 for i in range(32)],
        'width'         : 6.0
    }
}

# update parameters for low mass
params['system']['F']   = 1.07e4
params['system']['m']   = 5e-12
# wrap looper
looper = wrap_looper(
    looper_name='XLooper',
    func=get_func_quantum_correlation_measures(
        SystemClass=PhysRevLett_98_030405,
        params=params['solver'],
        steady_state=True
    ),
    params=params['looper'],
    params_system=params['system']
)
# extract values
xs  = looper.results['X']
M_0 = looper.results['V']

# update parameters for high mass
params['system']['F']   = 3.4e4
params['system']['m']   = 50e-12
# wrap looper
looper = wrap_looper(
    looper_name='XLooper',
    func=get_func_quantum_correlation_measures(
        SystemClass=PhysRevLett_98_030405,
        params=params['solver'],
        steady_state=True
    ),
    params=params['looper'],
    params_system=params['system']
)
# extract values
M_1 = looper.results['V']

# plotter
plotter = MPLPlotter(
    axes={},
    params=params['plotter']
)
plotter.update(
    xs=xs,
    vs=[M_0, M_1]
)
plotter.show()