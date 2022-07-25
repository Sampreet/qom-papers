#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the N-cell OM array system in Phys. Rev. Lett. **119**, 153901 (2017)."""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-08-09'
__updated__ = '2022-07-25'
__version__ = '0.8.5'

# dependencies
import numpy as np

# qom modules
from qom.systems import SOMASystem

class PhysRevLett_119_153901(SOMASystem):
    r"""Class to simulate the N-cell OM array system in Phys. Rev. Lett. **119**, 153901 (2017).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        N               (int) number of cells :math:`N`. Default is :math:`401`.
        Gamma_m         (float) normalized normalized mechanical dampling rate :math:`\Gamma_{m}` in Hertz. Default is :math:`190.0 / 9.5 \times 10^{9}`.
        g_0             (float) normalized optomechanical coupling strength :math:`g_{0}` in Hertz. Default is :math:`292 \times 10^{3} / 9.5 \times 10^{9}`.
        J               (float) normalized photon hopping strength :math:`J` in Hertz. Default is :math:`- 95 \times 10^{9} / 9.5 \times 10^{9}`.
        kappa           (float) normalized optical decay rate :math:`\kappa`. Default is :math:`3.8 \times 10^{6} / 9.5 \times 10^{9}`.
        Omega_m         (float) normalized mechanical freqency :math:`\Omega_{m}`. Default is :math:`1.0`.
        x_d             (float) envelope width :math:`x_{d}`. Default is :math:`40.0`.
        t_alphas        (str) type of initial occupancies ``t_alphas``. Default is "Gaussian".
        ============    ========================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """
    
    # default system parameters
    system_defaults = {
        'N'         : 401,
        'Gamma_m'   : 190 / 9.5e9,
        'g_0'       : 292e3 / 9.5e9,
        'J'         : - 95e9 / 9.5e9,
        'kappa'     : 3.8e6 / 9.5e9,
        'Omega_m'   : 1.0,
        'x_d'       : 40.0,
        't_alphas'  : 'Gaussian'
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_119_153901."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'physrevlett_119_153901'
        self.name = 'Array System in Phys. Rev. Lett. 119, 153901'

        # update parameters
        self.params = dict()
        for key in self.system_defaults:
            self.params[key] = params.get(key, self.system_defaults[key])

        # update properties
        self.num_modes = 2 * params.get('N', self.system_defaults['N'])

    def get_beta_rates(self, t, betas, c):
        """Method to obtain the mechanical mode rates.
        
        Parameters
        ----------
        t : float
            Time at which the rates are calculated.
        betas : list
            Mechanical modes of the system.
        c : list
            Optical modes and parameters of the system.

        Returns
        -------
        beta_rates : list
            Mechanical mode rates.
        """

        # extract frequently used variables
        N       = int(self.num_modes / 2)
        alphas  = c[:N]
        Gamma_m, g_0, J, _, Omega_m, x_d = c[N:]
        divisor = np.abs(J) / x_d**2
        
        # calculate mechanical mode rates
        beta_rates = list()
        for i in range(N):
            beta_rates.append((1j * g_0 * np.conjugate(alphas[i]) * alphas[i] - (Gamma_m + 1j * Omega_m) * betas[i]) / divisor)

        return beta_rates

    def get_ivc(self):
        r"""Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        m : list
            Initial values of modes.
            The :math:`2n` elements contain the optical and mechanical modes of each cavity.

        params : list
            Constants of the IVP.
            The elements contain the system parameters ``params`` in the following order:
            ========    =============================================
            index       parameter
            ========    =============================================
            0           normalized mechanical damping rate :math:`\Gamma_{m}`.
            1           normalized optomechanical coupling strength :math:`g_{0}`.
            2           normalized photon hopping strength :math:`J`.
            3           normalized optical decay rate :math:`\kappa`.
            4           normalized mechanical frequency :math:`\Omega_{m}`.
            5           envelope width :math:`x_{d}`.
            ========    =============================================
        """

        # extract frequently used variables
        _N          = self.params['N']
        Gamma_m     = self.params['Gamma_m']
        g_0         = self.params['g_0']
        J           = self.params['J']
        kappa       = self.params['kappa']
        Omega_m     = self.params['Omega_m']
        x_d         = self.params['x_d']
        t_alphas    = self.params['t_alphas']
        gamma       = 2 * g_0**2 * Omega_m / (Gamma_m**2 + Omega_m**2)
        R           = 1 / x_d / np.sqrt(gamma / np.abs(J))
        c_1         = 1
        c_2         = 2.5
        c_3         = 5
 
        # set initial optical amplitudes
        def gauss(X, mu, sigma):
            return np.exp(- np.power(X - mu, 2) / (2 * np.power(sigma, 2)))
        xs = np.linspace(-1, 1, _N)
        ys = xs / x_d / (xs[1] - xs[0])
        # default initial condition
        default = 1e3 - 2.5e3 * gauss(xs, 0, 0.1) + 7.5e3 * gauss(xs, 0, 0.05)

        # types of initial conditions
        alphas_0 = {
            'Gaussian' : R * (c_1 + c_2 * np.exp(- c_3 * ys**2)),
            'Lorentzian': R * c_2 / (c_1 + c_3 * ys**2)**2,
            'square': np.array([c_2 * R if np.abs(x) <= 0.15 else R for x in xs]),
            'soliton': R * (- 6 / (1 - 2 * np.cosh(ys)) - 1),
            'default': default
        }.get(t_alphas, default).tolist()

        # initial mode values as 1D list
        modes_0 = np.zeros(2 * _N, dtype=np.complex_).tolist()
        # add alphas
        for i in range(_N):
            modes_0[2 * i] = alphas_0[i]
        
        # constant parameters
        params = [Gamma_m, g_0, J, kappa, Omega_m, x_d]

        return modes_0, params

    def get_mode_rates(self, modes, params, t):
        """Function to obtain the rates of the optical and mechanical modes.

        Parameters
        ----------
        modes : list
            Values of the modes.
        params : list
            Constants parameters.
        t : float
            Time at which the rates are calculated.
        
        Returns
        -------
        mode_rates : list
            Rate for each mode.
        """

        # extract frequently used variables
        N       = int(self.num_modes / 2)
        Gamma_m, g_0, J, kappa, Omega_m, x_d = params
        alphas  = [modes[2 * i] for i in range(N)]
        betas   = [modes[2 * i + 1] for i in range(N)]
        divisor = np.abs(J) / x_d**2

        # mode rates
        mode_rates = list()
        for i in range(N):
            mode_rates.append((- 1j * J / 2 * (alphas[i - 1] if i > 0 else 0) + (- kappa / 2 + 2j * g_0 * np.real(betas[i]) + 1j * J) * alphas[i] - 1j * J / 2 * (alphas[i + 1] if i < N - 1 else 0)) / divisor)
            mode_rates.append((1j * g_0 * np.conjugate(alphas[i]) * alphas[i] - (Gamma_m + 1j * Omega_m) * betas[i]) / divisor)

        return mode_rates
        
    def get_op_d(self, params, ps, x_ss):
        """Method to get the dispersion operator.
        
        Parameters
        ----------
        params : list
            Parameters for the system.
        ps : numpy.ndarray
            Frequencies of the cells.
        x_ss : float
            Step-size of the cells.

        Returns
        -------
        op_D : numpy.ndarray
            Dispersion operator.
        """
        
        # extract frequently used variables
        _, _, J, kappa, _, x_d = params
        divisor = np.abs(J) / x_d**2
        
        # calculate dispersion operator
        op_D = (- kappa / 2 + 1j * J / 2 * x_ss**2 * ps**2) / divisor

        return op_D

    def get_op_n(self, params, betas):
        """Method to get the nonlinear operator.
        
        Parameters
        ----------
        params : list
            Parameters for the system.
        betas : list
            Mechanical modes of the system.

        Returns
        -------
        op_N : numpy.ndarray
            Nonlinear operator.
        """

        # extract frequently used variables
        _, g_0, J, _, _, x_d = params
        divisor = np.abs(J) / x_d**2

        # calculate nonlinear operator
        op_N = 2j * g_0 * np.real(betas) / divisor

        return op_N