# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver, NLSESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import PhysRevLett_119_153901

# frequently used variables
N = 400

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'update_betas'  : True,
        'use_sources'   : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 4.0,
        't_dim'         : 4001,
        'indices'       : list(range(0, 2 * N, 2))
    },
    'system'    : {
        'N'         : N,
        'Gamma_m'   : 190 / 9.5e9,
        'g_0'       : 292e3 / 9.5e9,
        'J'         : - 95e9 / 9.5e9,
        'kappa'     : 3.8e6 / 9.5e9,
        'Omega_m'   : 1.0,
        'x_d'       : 40.0,
        't_alphas'  : 'square'
    },
    'plotter'   : {
        'type'              : 'pcolormesh',
        'x_label'           : '$x$',
        'x_ticks'           : [0, 200, 400],
        'x_tick_labels'     : [-1, 0, 1],
        'y_label'           : '$\\tau$',
        'y_ticks'           : [0, 2, 4],
        'show_cbar'         : True,
        'cbar_title'        : '$| \\alpha | / 10^{3}$',
        'cbar_tick_labels'  : [0, 4, 8],
        'cbar_ticks'        : [0, 4000, 8000]
    }
}

# initialize logger
init_log()

# initialize system
system = PhysRevLett_119_153901(
    params=params['system']
)

# # initialize solver
# solver = HLESolver(
#     system=system,
#     params=params['solver']
# )
# initialize solver
solver = NLSESolver(
    system=system,
    params=params['solver']
)
# get times, mode indices and intensities
ys  = solver.T
xs  = list(range(N))
vs  = np.sqrt(solver.get_mode_intensities())

# plotter
plotter = MPLPlotter(
    axes={
        'X' : xs,
        'Y' : ys
    },
    params=params['plotter']
)
plotter.update(
    vs=vs
)
plotter.show(
    hold=True
)