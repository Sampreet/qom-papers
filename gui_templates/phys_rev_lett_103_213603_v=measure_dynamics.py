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
        't_min': 0,
        't_max': 50,
        't_dim': 5001
    },
    'system': {
        'P_1': 2e-3
    },
    'plotter': {
        'type': 'line',
        'x_label': '$t / \\tau$',
        'v_label': '$E_{N}$',
        'v_ticks': [0.0, 0.1, 0.2, 0.3, 0.4],
        'width': 6.0
    }
}