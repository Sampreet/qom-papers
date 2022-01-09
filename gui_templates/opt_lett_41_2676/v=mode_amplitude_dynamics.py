# dependencies
import numpy as np

# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        'mode_type': 'optical',
        'range_min': 0,
        'range_max': 401,
        't_min': 0.00,
        't_max': 4.00,
        't_dim': 401
    },
    'system': {
        'n': 301,
        'Gamma_m': 0.0,
        'g_0': 1e-4,
        'gamma': 0.0,
        'J': 2.0,
        'Omega': 1.0,
        'x_0': 20.0,
        'n_solitons': 1,
        'dist_norm': 0.0,
        'phi': 0.0,
        'order': 2
    },
    'plotter': {
        'type': 'surface_cz',
        'x_label': '$x$',
        'x_ticks': [0, 150, 300],
        'y_label': '$\\tau$',
        'y_ticks': [0, np.pi / 2, np.pi, 4],
        'y_tick_labels': ['0', '$\\pi / 2$', '$\\pi$', ''],
        'v_label': '$|\\alpha| (10^{3})$',
        'v_ticks': [0, 1000, 2000],
        'v_tick_labels': [0, 1, 2],
        'show_cbar': False,
        'width': 6.0
    }
}