# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.solvers.deterministic import HLESolver
from qom.solvers.measure import QCMSolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevLett_103_213603

# all parameters
params = {
    'solver'    : {
        'show_progress' : True,
        'cache'         : False,
        'ode_method'    : 'vode',
        't_min'         : 0.00,
        't_max'         : 50.00,
        't_dim'         : 5001,
        't_index_min'   : 3000,
        't_index_max'   : 3201,
        'indices'       : [(2, 2)]
    },
    'system'    : {
        'F'         : 1.4e4,
        'L'         : 25e-3,
        'lambda_l'  : 1064e-9,
        'm'         : 150e-12,
        'omega_m'   : 2.0 * np.pi * 1e6,
        'Ps'        : [10e-3, 2e-3],
        'Q'         : 1e6,
        'T'         : 0.1
    },
    'plotter'   : {
        'type'          : 'lines',
        'show_legend'   : True,
        'x_label'       : '$t / \\tau$',
        'x_ticks'       : [30.0, 30.5, 31.0, 31.5, 32.0],
        'y_name'        : '$P_{\\pm 1}$',
        'y_unit'        : '$\\mathrm{mW}$',
        'y_colors'      : ['b', 'g'],
        'v_label'       : '$\\langle \\delta q^{2} \\rangle$',
        'v_ticks'       : [0.0, 0.4, 0.8, 1.2, 1.6],
        'width'         : 6.0
    }
}

# initialize logger
init_log()

# update parameters without modulation
params['system']['Ps'][1] = 0.0
# initialize system
system = PhysRevLett_103_213603(
    params=params['system']
)

# initialize solver
solver = HLESolver(
    system=system,
    params=params['solver']
)
# get times
xs  = solver.get_times()
# get variances
vs_0 = solver.get_corr_indices()[:, 0]

# update parameters without modulation
params['system']['Ps'][1] = 2e-3
# initialize system
system = PhysRevLett_103_213603(
    params=params['system']
)

# initialize solver
solver = HLESolver(
    system=system,
    params=params['solver']
)
# get variances
vs_1 = solver.get_corr_indices()[:, 0]

# plotter
plotter = MPLPlotter(
    axes={
        'Y' : [0, 2]
    },
    params=params['plotter']
)
plotter.update(
    xs=xs,
    vs=[vs_0, vs_1]
)
plotter.show()