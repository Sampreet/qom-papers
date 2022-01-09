# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevA101_053836_RWA

# all parameters
params = {
    'looper': {
        'show_progress': False,
        'X': {
            'var': 'G_p1_norm',
            'min': 0.0,
            'max': 1.0,
            'dim': 1001
        },
        'Y': {
            'var': 'G_0_norm',
            'min': 0.0,
            'max': 2.0,
            'dim': 101
        }
    },
    'solver': {
        'cache': False,
        'method': 'zvode',
        'measure_type': 'corr_ele',
        'idx_e': (2, 2),
        'range_min': 900,
        'range_max': 1001,
        't_min': 0.0,
        't_max': np.pi * 100.0,
        't_dim': 1001
    },
    'system': {
        'Delta_a': 1.0,
        'G_0': 0.1,
        'G_p1': 0.05,
        'G_m1': 0.01,
        'gamma_m': 1e-6,
        'kappa': 0.1,
        'omega_m': 1.0,
        'Omega': 2.0,
        'n_a': 0,
        'n_m': 10
    },
    'plotter': {
        'type': 'lines',
        'x_label': '$G_{0} / \\kappa$',
        'x_ticks': [0.2 + 0.2 * i for i in range(10)],
        'x_tick_labels': ['{:0.1f}'.format(0.2 + 0.2 * i) for i in range(10)],
        'y_colors': ['b', 'r'],
        'v_label': 'Optimal $G_{1} / G_{0}$',
        'v_ticks': [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'show_legend': True,
        'width': 8.0,
        'height': 4.0
    }
}

# function to calculate average measure in decibels
def func_mav_db(system_params, val, logger, results):
    # update parameters
    system_params['G_0'] = system_params['G_0_norm'] * system_params['kappa']
    system_params['G_p1'] = system_params['G_p1_norm'] * system_params['G_0']
    # initialize system
    system = PhysRevA101_053836_RWA(params=system_params)
    # get average measure
    value = system.get_measure_average(solver_params=params['solver'])
    # convert to dB scale
    value = 10 * np.log10(0.5 / value)
    # update results
    results.append((val, value))

# wrapper
looper = wrap_looper(SystemClass=None, params=params, func=func_mav_db, looper='xy_looper', file_path='data/phys_rev_a_101_053836_RWA/mav_db')
xs = looper.axes['Y']['val']
vs_10 = [looper.axes['X']['val'][np.argmax(v)] for v in looper.results['V']]

# wrapper
params['system']['n_m'] = 100
looper = wrap_looper(SystemClass=None, params=params, func=func_mav_db, looper='xy_looper', file_path='data/phys_rev_a_101_053836_RWA/mav_db')
vs_100 = [looper.axes['X']['val'][np.argmax(v)] for v in looper.results['V']]

# plotter
plotter = MPLPlotter(axes={
    'X': xs,
    'Y': ['$n_{m} = 10$', '$n_{m} = 100$']
}, params=params['plotter'])
plotter.update(xs=[xs, xs], vs=[vs_10, vs_100])
plotter.show(True)