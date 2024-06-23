# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.solvers.measure import get_Wigner_distributions_two_mode
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_063846_00

# frequently used variables
xs = np.linspace(-5, 5, 101)
ys = np.linspace(-5, 5, 101)

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'X'                 : {
            'var': 'J_norm',
            'val': [1.0, 0.5]
        }
    },
    'solver'    : {
        'show_progress' : False,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 1.0,
        't_dim'         : 1001,
        't_index_min'   : 750,
        't_index_max'   : 750,
        'indices'       : [(0, 0), (1, 0)],
        'wigner_xs'     : xs,
        'wigner_ys'     : ys
    },
    'system'    : {
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    },
    'plotter'   : {
        'type'              : 'pcolormesh',
        'x_label'           : '$q_{1}$',
        'x_tick_position'   : 'both-out',
        'x_ticks'           : [-5, 0, 5],
        'x_ticks_minor'     : [i - 5 for i in range(11)],
        'y_label'           : '$q_{2}$',
        'y_tick_position'   : 'both-out',
        'y_ticks'           : [-5, 0, 5],
        'y_ticks_minor'     : [i - 5 for i in range(11)],
        'show_cbar'         : True,
        'cbar_label'        : '$W$',
        'cbar_ticks'        : [0, 0.6e-2, 1.2e-2],
        'width'             : 5.5,
        'height'            : 4.8
    }
}

# function to obtain the wigners
def func(system_params):
    # initialize system
    system = PhysRevA_100_063846_00(
        params=system_params
    )
    # initialize solver
    solver = HLESolver(
        system=system,
        params=params['solver']
    )
    # get correlations
    Corrs = solver.get_corrs()
    # get wigners
    W = get_Wigner_distributions_two_mode(
        Corrs=Corrs,
        params=params['solver']
    )[0]
    return W

# looper
looper = wrap_looper(
    looper_name='XLooper',
    func=func,
    params=params['looper'],
    params_system=params['system']
)
# plotter
for i in range(2):
    params['plotter']['title'] = '$J = {} \\Gamma$'.format(params['looper']['X']['val'][i])
    plotter = MPLPlotter(
        axes={
            'X' : xs,
            'Y' : ys
        },
        params=params['plotter']
    )
    plotter.update(
        vs=looper.results['V'][i]
    )
    plotter.show()

