# dependencies
import os 
import sys

# qom modules
from qom.ui.plotters import MPLPlotter
from qom.utils.looper import wrap_looper

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-papers')))
# import system
from systems import PhysRevLett114_013601

# all parameters
params = {
    'looper': {
        'X': {
            'var': 'Delta_norm',
            'min': -1.5,
            'max': 0.0,
            'dim': 76
        },
        'Y': {
            'var': 'P',
            'min': 1.0,
            'max': 1.6,
            'dim': 61
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
        'Gamma_norm': 1e-3, 
        'kappa_norm': 1.0,
        'P': 1.0
    },
    'plotter': {
        'type': 'pcolormesh',
        'x_label': '$\\Delta / \\Omega_{m}$',
        'x_ticks': [-1.5, -1.0, -0.5, 0.0],
        'y_label': '$P$',
        'y_ticks': [1.0, 1.2, 1.4, 1.6],
        'cbar_title': '$\\lambda_{max}$',
        'cbar_ticks': [-0.005, 0.00, 0.005],
        'cbar_position': 'top'
    }
}

# wrapper
looper = wrap_looper(SystemClass=PhysRevLett114_013601, params=params, func='les', looper='xy_looper', file_path='data/phys_rev_lett_114_013601/lyapunov_exponents')

# get maximum lyapunov exponent zones
vs = [[max(e) for e in r] for r in looper.results['V']]

# plotter
plotter = MPLPlotter(axes={
    'X': looper.axes['X'],
    'Y': looper.axes['Y']
}, params=params['plotter'])
plotter.update(vs=vs)
plotter.show(True)