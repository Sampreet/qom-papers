# all parameters
params = {
    'solver': {
        'show_progress': True,
        'cache': True,
        'method': 'zvode',
        't_min': 0.0,
        't_max': 4.0,
        't_dim': 4001
    },
    'system': {
        't_alphas': 'Gaussian'
    },
    'plotter': {
        'type': 'pcolormesh',
        'palette': 'RdBu_r',
        'x_label': '$x$',
        'x_ticks': [0, 200, 400],
        'x_tick_labels': [-1, 0, 1],
        'y_label': '$\\tau$',
        'y_ticks': [0, 2, 4],
        'show_cbar': True,
        'cbar_title': '$|\\alpha| (10^{3})$',
        'cbar_ticks': [0, 4000, 8000],
        'cbar_tick_labels': [0, 4, 8]
    }
}