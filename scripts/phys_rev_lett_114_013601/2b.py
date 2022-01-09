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
        'show_progress': True,
        'X': {
            'var': 'Delta_norm',
            'min': -1.2,
            'max': -0.4,
            'dim': 161
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
        'gamma_norm': 1e-3,
        'kappa_norm': 1.0,
        'P': 1.4
    },
    'plotter': {
        'type': 'line',
        'x_label': '$\\Delta / \\Omega_{m}$',
        'x_bound': 'both',
        'x_ticks': [-1.2, -0.8, -0.4],
        'v_label': '$\\lambda_{max}$',
        'v_bound': 'both',
        'v_ticks': [-0.2, -0.1, 0.0, 0.1]
    }
}

# wrapper
looper = wrap_looper(SystemClass=PhysRevLett114_013601, params=params, func='les', looper='x_looper', file_path='data/phys_rev_lett_114_013601/lyapunov_exponents')

# get maximum lyapunov exponent zones
vs = [max(e) for e in looper.results['V']]

# plotter
plotter = MPLPlotter(axes={
    'X': looper.axes['X']
}, params=params['plotter'])
plotter.update(xs=looper.axes['X']['val'], vs=vs)
plotter.show(True)