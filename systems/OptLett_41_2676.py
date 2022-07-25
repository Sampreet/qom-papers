#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the n-cell OM array system in Opt. Lett. **41**, 2676 (2016)."""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-08-15'
__updated__ = '2022-07-24'
__version__ = '0.8.5'

# dependencies
import numpy as np

# qom modules
from qom.systems import SOMASystem

class OptLett_41_2676(SOMASystem):
    r"""Class to simulate the n-cell OM array system in Opt. Lett. **41**, 2676 (2016).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        n               (int) number of cells :math:`n`. Default is :math:`151`.
        Gamma_m         (float) normalized mechanical dampling rate :math:`\Gamma_{m}`. Default is :math:`0.0`.
        g_0             (float) normalized optomechanical coupling strength :math:`g_{0}`. Default is :math:`10^{-4}`.
        gamma           (float) normalized optical decay rate :math:`\gamma`. Default is :math:`0.0`.
        J               (float) normalized photon hopping strength :math:`J`. Default is :math:`2.0`.
        Omega           (float) normalized mechanical freqency :math:`\Omega`. Default is :math:`1.0`.
        x_0             (float) envelope width :math:`x_{0}`. Default is :math:`10.0`.
        n_solitons      (int) number of solitions ``n_solitons``. Default is :math:`1`.
        dist_norm       (float) normalized initial distance between the solitons :math:`d / x_0`. Default is :math:`0.0`.
        phi             (float) phase difference between the solitions :math:`\phi`. Default is :math:`0.0`.
        order           (int) order of the soliton :math:`N`. Default is :math:`1`.
        ============    ========================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default system parameters
    system_defaults = {
        'n'         : 151,
        'Gamma_m'   : 0.0,
        'g_0'       : 1e-4,
        'gamma'     : 0.0,
        'J'         : 2.0,
        'Omega'     : 1.0,
        'x_0'       : 10.0,
        'n_solitons': 1,
        'dist_norm' : 0.0,
        'phi'       : 0.0,
        'order'     : 1
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for OptLett_41_2676."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)
        
        # set attributes
        self.code = 'optlett_41_2676'
        self.name = 'Array System in Opt. Lett. 41, 2676'

        # update parameters
        self.params = dict()
        for key in self.system_defaults:
            self.params[key] = params.get(key, self.system_defaults[key])

        # update properties
        self.num_modes = 2 * params.get('n', self.system_defaults['n'])

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
        n      = int(self.num_modes / 2)
        alphas = c[:n]
        Gamma_m, g_0, _, J, Omega, x_0 = c[n:]
        divisor = J / x_0**2

        # calculate mechanical mode rates
        beta_rates = list()
        for i in range(n):
            beta_rates.append((1j * g_0 * np.conjugate(alphas[i]) * alphas[i] - (Gamma_m + 1j * Omega) * betas[i]) / divisor)

        return beta_rates

    def get_ivc(self):
        r"""Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        modes_0 : list
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
            2           normalized optical decay rate :math:`\gamma`.
            3           normalized photon hopping strength :math:`J`.
            4           normalized mechanical frequency :math:`\Omega`.
            5           envelope width :math:`x_{0}`.
            ========    =============================================
        """

        # extract frequently used variables
        n           = self.params['n']
        Gamma_m     = self.params['Gamma_m']
        g_0         = self.params['g_0']
        gamma       = self.params['gamma']
        J           = self.params['J']
        Omega       = self.params['Omega']
        x_0         = self.params['x_0']
        n_solitons  = self.params['n_solitons']
        dist_norm   = self.params['dist_norm']
        phi         = self.params['phi']
        order       = self.params['order']
 
        # set default optical amplitudes
        temp = order * np.sqrt(Omega * J / 2 / g_0**2 / x_0**2) / np.cosh(np.linspace(- (n - 1) / 2, (n - 1) / 2, n) / x_0)

        # double solitons
        if int(n_solitons) == 2:
            offset = int(dist_norm * x_0 / 2)
            alpha_0s = np.roll(temp, - offset) + np.roll(temp, offset) * np.exp(1j * phi)
        # single soliton
        else:
            alpha_0s = temp

        # initial mode values as 1D list
        modes_0 = np.zeros(self.num_modes, dtype=np.complex_).tolist()
        # add alphas
        for i in range(n):
            modes_0[2 * i] = alpha_0s[i]
        
        # constant parameters
        params = [Gamma_m, g_0, gamma, J, Omega, x_0]

        return modes_0, params

    def get_mode_rates(self, modes, params, t=None):
        """Function to obtain the rates of the optical and mechanical modes.

        Parameters
        ----------
        modes : list
            Values of the modes.
        params : list
            Constants parameters.
        t : float, optional
            Time at which the rates are calculated.
        
        Returns
        -------
        mode_rates : list
            Rate for each mode.
        """

        # extract frequently used variables
        n       = int(self.num_modes / 2)
        Gamma_m, g_0, gamma, J, Omega, x_0 = params
        alphas  = [modes[2 * i] for i in range(n)]
        betas   = [modes[2 * i + 1] for i in range(n)]
        divisor = J / x_0**2
        
        # mode rates
        mode_rates = list()
        for i in range(0, n):
            # optical mode
            mode_rates.append((1j * J / 2 * (alphas[i - 1] if i > 0 else 0) + (- gamma + 2j * g_0 * np.real(betas[i]) - 1j * J) * alphas[i] + 1j * J / 2 * (alphas[i + 1] if i < n - 1 else 0)) / divisor)
            # mechanical mode
            mode_rates.append((1j * g_0 * np.conjugate(alphas[i]) * alphas[i] - (Gamma_m + 1j * Omega) * betas[i]) / divisor)

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
        _, _, gamma, J, _, x_0 = params
        divisor = J / x_0**2

        # calculate dispersion operator
        op_D = (- gamma - 1j * J / 2 * x_ss**2 * ps**2) / divisor

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
        _, g_0, _, J, _, x_0 = params
        divisor = J / x_0**2

        # calculate nonlinear operator
        op_N = 2j * g_0 * np.real(betas) / divisor

        return op_N