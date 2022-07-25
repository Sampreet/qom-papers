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
from systems import PhysRevA_101_053836

# all parameters
params = {
    'looper': {
        'show_progres_yz': True,
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
        },
        'Z': {
            'var': 'ns',
            'idx': 1,
            'val': [10.0, 100.0]
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
        'Delta_a_norm': 1.0,
        'G_norms': [0.1, 0.01, 0.05],
        'gamma_m_norm': 1e-6,
        'kappa_norm': 0.1,
        'ns': [0.0, 10.0],
        'Omega_norm': 2.0,
        't_rwa': True
    },
    'plotter': {
        'type': 'lines',
        'x_label': '$G_{0} / \\kappa$',
        'x_ticks': [0.2 + 0.2 * i for i in range(10)],
        'x_tick_labels': ['{:0.1f}'.format(0.2 + 0.2 * i) for i in range(10)],
        'y_colors': ['b', 'r'],
        'v_label': '$\\langle \\delta \\tilde{X}_{b}^{2} \\rangle_{max}$ (dB)',
        'v_ticks': [0, 5, 10, 15, 20, 25],
        'show_legend': True,
        'width': 8.0,
        'height': 4.0
    }
}

# function to calculate average measure in decibels
def func_mav_db(system_params, val, logger, results):
    # update parameters
    system_params['G_norms'][0] = system_params['G_0_norm'] * system_params['kappa_norm']
    system_params['G_norms'][2] = system_params['G_p1_norm'] * system_params['G_norms'][0]
    # initialize system
    system = PhysRevA_101_053836(params=system_params)
    # get average measure
    value = system.get_measure_average(solver_params=params['solver'])
    # convert to dB scale
    value = 10 * np.log10(0.5 / value)
    # update results
    results.append((val, value))

# wrapper
looper = wrap_looper(SystemClass=None, params=params, func=func_mav_db, looper='XYZLooper', file_path_prefix='data/physreva_101_053836/5_mav_db')
xs = looper.axes['Y']['val']
vs_10 = [np.max(v) for v in looper.results['V'][0]]
vs_100 = [np.max(v) for v in looper.results['V'][1]]

# plotter
plotter = MPLPlotter(axes={
    'X': xs,
    'Y': ['$n_{m} = 10$', '$n_{m} = 100$']
}, params=params['plotter'])
plotter.update(xs=[xs, xs], vs=[vs_10, vs_100])
plotter.show(True)