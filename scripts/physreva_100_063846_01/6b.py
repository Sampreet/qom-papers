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
from systems import PhysRevA_100_063846_01

# function to obtain the normalized frequencies
def func(system_params):
    # initialize system
    system = PhysRevA_100_063846_01(
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
            'val'   : np.linspace(0.3, 0.4, 1001)
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
        'colors'        : ['r', 'm', 'b'],
        'x_label'       : '$J / \\Gamma$',
        'x_ticks'       : [0.3, 0.35, 0.4],
        'x_ticks_minor' : [0.0125 * i + 0.3 for i in range(9)],
        'v_label'       : 'Im$[ \\omega_{n} ]$',
        'v_ticks'       : [-0.4, 0.0, 0.4],
        'v_ticks_minor' : [0.1 * i - 0.4 for i in range(9)],
        'width'         : 4.8,
        'height'        : 4.8
    }
)
plotter.update(
    vs=np.transpose(looper.results['V']),
    xs=looper.results['X']
)
plotter.show()