# dependencies
import os 
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import NewJPhys_22_013049

# all parameters
params = {
    'looper': {
        'X': {
            'var': 'Delta_norm',
            'min': -1.0,
            'max': 1.0,
            'dim': 101
        },
        'Y': {
            'var': 'P',
            'min': 0.0,
            'max': 0.5,
            'dim': 51
        }
    },
    'solver': {
        'cache': True,
        'method': 'zvode',
        'method_le': 'svd',
        'num_iters': 10000,
        't_min': 0.0,
        't_max': 1000.0,
        't_dim': 10001
    },
    'system': {
        'Delta_norm': 0.0, 
        'gamma_norm': 1e-4,
        'kappa_norm': 0.1,
        'P': 0.1
    },
    'plotter': {
        'type': 'pcolormesh',
        'x_label': '$\\Delta / \\Omega_{m}$',
        'x_ticks': [-1.0, 0.0, 1.0],
        'y_label': '$P$',
        'y_ticks': [0.00, 0.25, 0.50],
        'show_cbar': True,
        'cbar_title': '$\\lambda_{max}$',
        'cbar_ticks': [-0.005, 0.00, 0.005],
        'cbar_position': 'top'
    }
}

# wrapper
looper = wrap_looper(SystemClass=NewJPhys_22_013049, params=params, func='les', looper='XYLooper', file_path_prefix='data/newjphys_22_013049/5c_les')

# get maximum lyapunov exponent zones
vs = [[max(e) for e in r] for r in looper.results['V']]

# plotter
plotter = MPLPlotter(axes={
    'X': looper.axes['X'],
    'Y': looper.axes['Y']
}, params=params['plotter'])
plotter.update(vs=vs)
plotter.show(True)