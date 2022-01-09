# dependencies
import numpy as np

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        'measure_type': 'entan_ln',
        'idx_e': (0, 1),
        'range_min': 3000,
        'range_max': 3201,
        't_min': 0.00,
        't_max': 50.00,
        't_dim': 5001
    },
    'system': {
        'F': 1.4e4,
        'L': 25e-3,
        'lambda_l': 1064e-9,
        'm': 150e-12,
        'omega_m': 2e6 * np.pi,
        'P_0': 10e-3,
        'P_1': 2e-3,
        'Q': 1e6,
        'T': 0.1
    },
    'plotter': {
        'type': 'line',
        'x_label': '$t / \\tau$',
        'v_label': '$E_{N}$',
        'v_ticks': [0.0, 0.1, 0.2, 0.3, 0.4],
        'width': 6.0
    }
}