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
from systems import PhysRevA_100_063846_00

# all parameters
params = {
    'looper'    : {
        'show_progress' : True,
        'X'             : {
            'var': 'J_norm',
            'val': [0.5, 0.75, 1.0]
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 1.0,
        't_dim'         : 1001,
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
        'type'          : 'lines',
        'colors'        : ['m', 'r', 'b'],
        'x_label'       : '$\\Gamma t$',
        'x_ticks'       : [0.0, 0.5, 1.0],
        'x_ticks_minor' : [0.1 * i for i in range(11)],
        'y_name'        : '$J$',
        'y_unit'        : '$\\Gamma$',
        'v_label'       : '$E_{N}$',
        'v_ticks'       : [0, 1, 2],
        'v_ticks_minor' : [0.2 * i for i in range(11)],
        'show_legend'   : True,
        'width'         : 4.8,
        'height'        : 4.8
    }
}

# looper
looper = wrap_looper(
    looper_name='XLooper',
    func=get_func_quantum_correlation_measures(
        SystemClass=PhysRevA_100_063846_00,
        params=params['solver'],
        steady_state=False
    ),
    params=params['looper'],
    params_system=params['system']
)
vs = [[m[0] for m in row] for row in looper.results['V']]
xs = np.linspace(params['solver']['t_min'], params['solver']['t_max'], params['solver']['t_dim'])
ys = looper.results['X']

# plotter
plotter = MPLPlotter(
    axes={
        'X' : xs,
        'Y' : ys
    },
    params=params['plotter']
)
plotter.update(
    vs=vs,
    xs=xs,
    ys=ys
)
plotter.show()