# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.loopers import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_063846_00

# function to obtain the normalized frequencies
def func(system_params):
    # initialize system
    system = PhysRevA_100_063846_00(
        params=system_params
    )
    return np.imag(system.get_omega_norms(
        c=None
    ))

# looper
looper = wrap_looper(
    looper_name='XLooper',
    func=func,
    params={
        'show_progress' : True,
        'X'             : {
            'var'   : 'J_norm',
            'val'   : np.linspace(0.0, 1.0, 1001)
        }
    },
    params_system={
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    }
)
# plotter
plotter = MPLPlotter(
    axes={
        'X' : looper.results['X']
    },
    params={
        'type'          : 'lines',
        'colors'        : ['b', 'r'],
        'x_label'       : '$J / \\Gamma$',
        'x_ticks'       : [0, 0.5, 1.0],
        'x_ticks_minor' : [0.125 * i for i in range(9)],
        'v_label'       : 'Im$[ \\omega_{\\pm} ]$',
        'v_ticks'       : [-1, 0, 1],
        'v_ticks_minor' : [0.25 * i - 1.0 for i in range(9)],
        'width'         : 4.8,
        'height'        : 4.8
    }
)
plotter.update(
    vs=np.transpose(looper.results['V']),
    xs=looper.results['X']
)
plotter.show()