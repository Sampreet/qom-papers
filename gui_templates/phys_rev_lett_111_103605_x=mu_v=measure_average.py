# all parameters
params = {
    'looper': {
        'show_progress': True,
        'X': {
            'var': 'mu',
            'min': 0.00,
            'max': 0.04,
            'dim': 11
        }
    },
    'solver': {
        'cache': True,
        'method': 'zvode',
        'measure_type': 'sync_p',
        'idx_e': (1, 3),
        'range_min': 371,
        'range_max': 1001,
        't_min': 0,
        't_max': 100,
        't_dim': 1001
    },
    'system': {
        'omega_2': 1
    },
    'plotter': {
        'type': 'scatter',
        'palette': 'RdBu_r',
        'x_label': '$\\mu$',
        'x_ticks': [0.00, 0.01, 0.02, 0.03, 0.04],
        'v_label': '$\\bar{S}_{p}$',
        'v_ticks': [0, 0.1, 0.2, 0.3]
    }
}