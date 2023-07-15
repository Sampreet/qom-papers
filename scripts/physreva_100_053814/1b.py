# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import LLESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevA_100_053814

# frequently used variables
N       = 2**12
tau_max = 100.0

# all parameters
params = {
    'solver': {
        'show_progress' : True,
        'cache'         : True,
        't_min'         : 0.0,
        't_max'         : 20.0,
        't_dim'         : 201,
        'indices'       : [i * 2 for i in range(N)]
    },
    'system': {
        'N'         : N,
        'Delta'     : 3.0,
        'S_0'       : np.sqrt(3.5),
        'tau_max'   : tau_max,
        't_alphas'  : 'sech'
    },
    'plotter': {
        'type'              : 'surface_cz',
        'palette'           : 'RdBu_r',
        'x_label'           : '$\\tau$',
        'x_tick_pad'        : 2,
        'x_ticks'           : [-4, 0, 4],
        'x_ticks_minor'     : [-2, 2],
        'y_label'           : '$t$',
        'y_tick_pad'        : 2,
        'y_ticks'           : [0, 10, 20],
        'y_ticks_minor'     : [5, 15],
        'v_label'           : '$| u |^{2}$',
        'v_tick_pad'        : 2,
        'v_ticks'           : [0, 4, 8, 12],
        'v_ticks_minor'     : [2, 6, 10, 14],
        'show_cbar'         : False,
        'view_aspect'       : [1.5, 1.0, 1.0],
        'view_elevation'    : 32.0,
        'view_rotation'     : -45.0,
        'show_cbar'         : False,
        'cbar_position'     : 'left',
        'cbar_ticks'        : [0, 4, 8],
        'width'             : 5.0,
        'height'            : 4.4,
        'annotations'       : [{
            'text'  : '(b)',
            'xy'    : (0.19, 0.78)
        }],
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevA_100_053814(
    params=params['system']
)

# initialize solver
solver = LLESolver(
    system=system,
    params=params['solver']
)
# get times and intensities
T           = solver.T
Intensities = solver.get_mode_intensities()
# extract required region
X       = np.linspace(- tau_max, tau_max, N)
idxs    = np.argwhere((X > -4) & (X < 4)).transpose()[0]
xs      = X[idxs[0]:idxs[-1] + 1]
vs      = Intensities[:, idxs]

# plotter
plotter = MPLPlotter(
    axes={
        'X': xs,
        'Y': T
    },
    params=params['plotter']
)
plotter.update(
    vs=vs
)
plotter.show()