# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems import OptLett_41_2676

# frequently used variables
n = 300

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.0,
        't_max'         : 100.0,
        't_dim'         : 1001,
        'indices'       : list(range(0, 2 * n, 2))
    },
    'system'    : {
        'n'             : n,
        'Gamma_m'       : 0.0,
        'g_0'           : 1e-4,
        'gamma'         : 0.0,
        'J'             : 2.0,
        'Omega'         : 1.0,
        'x_0'           : 10.0,
        'n_solitons'    : 2,
        'dist_norm'     : 7.0,
        'phi'           : np.pi / 2.0,
        'order'         : 1
    },
    'plotter'   : {
        'type'          : 'surface_cz',
        'x_label'       : '$x$',
        'x_tick_pad'    : 2,
        'x_ticks'       : [0, 150, 300],
        'y_label'       : '$\\tau$',
        'y_tick_pad'    : 2,
        'y_ticks'       : [0, 50, 100],
        'v_label'       : '$|\\alpha|$\n$(10^{3})$',
        'v_tick_pad'    : 2,
        'v_ticks'       : [0, 1000, 2000],
        'v_tick_labels' : [0, 1, 2],
        'show_cbar'     : False,
        'width'         : 6.0
    }
}

# initialize logger
init_log()

# initialize system
system = OptLett_41_2676(
    params=params['system']
)

# initialize solver
solver = HLESolver(
    system=system,
    params=params['solver']
)
# get times, mode indices and intensities
ys  = solver.T
xs  = list(range(n))
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
plotter.show()