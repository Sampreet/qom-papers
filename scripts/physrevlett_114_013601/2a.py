# dependencies
import numpy as np
import os 
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.solvers.measure import get_bifurcation_amplitudes
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_114_013601

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/physrevlett_114_013601/2a',
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
        't_index_min'   : 99371,
        't_index_max'   : 100001,
        'indices'       : [1]
    },
    'system'    : {
        'Delta_norm'    : 0.0, 
        'Gamma_norm'    : 1e-3,
        'kappa_norm'    : 1.0,
        'P'             : 1.4
    },
    'plotter'   : {
        'type'          : 'scatter',
        'colors'        : ['k'] * 20,
        'sizes'         : [0.1] * 20,
        'styles'        : ['.'] * 20,
        'x_label'       : '$\\Delta / \\Omega$',
        'x_ticks'       : [-1.2, -1.0, -0.8, -0.6, -0.4],
        'x_ticks_minor' : [i * 0.05 - 1.2 for i in range(17)],
        'v_label'       : 'amplitudes',
        'v_limits'      : [0.25, 2.0],
        'v_ticks'       : [0.5, 1.0, 1.5, 2.0],
        'v_ticks_minor' : [i * 0.25 + 0.25 for i in range(8)],
        'width'         : 6.0,
        'height'        : 4.0
    }
}

# wrapper function for parallel processes
def func(system_params):
    # initialize system
    system = PhysRevLett_114_013601(
        params=system_params
    )

    # initialize solver
    solver = HLESolver(
        system=system,
        params=params['solver']
    )
    # get modes
    Modes = solver.get_mode_indices()

    # return bifurcation amplitudes of the mechanical position
    return get_bifurcation_amplitudes(
        Modes=Modes
    )[0]

if __name__ == '__main__':
    # wrap looper
    looper = wrap_looper(
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