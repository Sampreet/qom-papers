#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the n-cell OM array system in Opt. Lett. **41**, 2676 (2016)."""

__authors__ = ["Sampreet Kalita"]
__toolbox__ = "qom-v1.1.0"
__created__ = "2021-08-15"
__updated__ = "2025-03-11"
__all__     = ['OptLett_41_2676']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class OptLett_41_2676(BaseSystem):
    r"""Class to simulate the n-cell OM array system in Opt. Lett. **41**, 2676 (2016).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        n               (*int*) number of cells :math:`n`. Default is :math:`150`.
        Gamma_m         (*float*) normalized mechanical dampling rate :math:`\Gamma_{m}`. Default is :math:`0.0`.
        g_0             (*float*) normalized optomechanical coupling strength :math:`g_{0}`. Default is :math:`10^{-4}`.
        gamma           (*float*) normalized optical decay rate :math:`\gamma`. Default is :math:`0.0`.
        J               (*float*) normalized photon hopping strength :math:`J`. Default is :math:`2.0`.
        Omega           (*float*) normalized mechanical freqency :math:`\Omega`. Default is :math:`1.0`.
        x_0             (*float*) envelope width :math:`x_{0}`. Default is :math:`10.0`.
        n_solitons      (*int) number of solitions `n_solitons`. Default is :math:`1`.
        dist_norm       (*float*) normalized initial distance between the solitons :math:`d / x_0`. Default is :math:`0.0`.
        phi             (*float*) phase difference between the solitions :math:`\phi`. Default is :math:`0.0`.
        order           (*int*) order of the soliton :math:`N`. Default is :math:`1`.
        ============    ========================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'n'         : 150,
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
        super().__init__(
            params=params,
            name='OptLett_41_2676',
            desc="Array System in Opt. Lett. 41, 2676",
            num_modes=2 * params.get('n', self.system_defaults['n']),
            cb_update=cb_update
        )

    def get_beta_rates(self, modes, c, t):
        """Method to obtain the rates of change of the mechanical modes.

        Parameters
        ----------
        modes : *numpy.ndarray*
            Classical modes.
        c : *numpy.ndarray*
            Derived constants and controls.
        t : *float*
            Time at which the values are calculated.

        Returns
        -------
        beta_rates : *numpy.ndarray*
            Rates of change of the mechanical modes.
        """

        # extract frequently used variables
        divisor = self.params['J'] / self.params['x_0']**2

        # return rates
        return (1.0j * self.params['g_0'] * np.conjugate(modes[::2]) * modes[::2] - (self.params['Gamma_m'] + 1.0j * self.params['Omega']) * modes[1::2]) / divisor
    
    def get_ivc(self):
        """Method to obtain the initial values of the modes, correlations and derived constants and controls.
        
        Returns
        -------
        iv_modes : *numpy.ndarray*
            Initial values of the classical modes.
        iv_corrs : *numpy.ndarray*
            Initial values of the quantum correlations.
        c : *numpy.ndarray*
            Derived constants and controls.
        """

        # extract frequently used variables
        n   = self.params['n']
        x_0 = self.params['x_0']
 
        # set default optical amplitudes
        temp = self.params['order'] * np.sqrt(self.params['Omega'] * self.params['J'] / 2 / self.params['g_0']**2 / x_0**2) / np.cosh(np.linspace(- (n - 1.0) / 2.0, (n - 1.0) / 2.0, n) / x_0)

        # initial values of the modes
        iv_modes = np.zeros(self.num_modes, dtype=np.complex128)
        # double solitons
        if int(self.params['n_solitons']) == 2:
            offset          = int(self.params['dist_norm'] * x_0 / 2.0)
            iv_modes[::2]   = np.roll(temp, - offset) + np.roll(temp, offset) * np.exp(1.0j * self.params['phi'])
        # single soliton
        else:
            iv_modes[::2] = temp

        return iv_modes, None, None

    def get_mode_rates(self, modes, c, t):
        """Method to obtain the rates of change of the modes.

        Parameters
        ----------
        modes : *numpy.ndarray*
            Classical modes.
        c : *numpy.ndarray*
            Derived constants and controls.
        t : *float*
            Time at which the values are calculated.
        
        Returns
        -------
        mode_rates : *numpy.ndarray*
            Rates of change of the modes.
        """

        # extract frequently used variables
        g_0     = self.params['g_0']
        J       = self.params['J']
        alphas  = modes[::2]
        betas   = modes[1::2]
        divisor = J / self.params['x_0']**2

        # initialize mode rates
        mode_rates = np.zeros_like(modes, dtype=np.complex128)
        
        # update rates for optical mode
        for i in range(len(alphas)):
            mode_rates[2 * i] = (1.0j * J / 2.0 * (alphas[i - 1] if i > 0 else 0.0) + (- self.params['gamma'] + 2.0j * g_0 * np.real(betas[i]) - 1.0j * J) * alphas[i] + 1.0j * J / 2.0 * (alphas[i + 1] if i < len(alphas) - 1 else 0.0)) / divisor
        # update rates for mechanical modes
        mode_rates[1::2] = (1.0j * g_0 * np.conjugate(alphas) * alphas - (self.params['Gamma_m'] + 1.0j * self.params['Omega']) * betas) / divisor

        return mode_rates
        
    def get_coeffs_dispersion(self, modes, c, t):
        r"""Method to get the coefficients of :math:``( i \omega )^{j}`` in the dispersions.
        
        Parameters
        ----------
        modes : *numpy.ndarray*
            Classical modes.
        c : *numpy.ndarray*
            Derived constants and controls.
        t : *float*
            Time at which the values are calculated.

        Returns
        -------
        coeffs : *numpy.ndarray*
            Coefficients in the dispersion operator.
        """

        # extract frequently used variables
        J       = self.params['J']
        divisor = J / self.params['x_0']**2

        # return coefficients
        return np.array([0.0, 0.0, - 1.0j * J / 2 * 1.0**2 / divisor], dtype=np.complex128)

    def get_nonlinearities(self, modes, c, t):
        """Method to get the nonlinearities.
        
        Parameters
        ----------
        modes : *numpy.ndarray*
            Classical modes.
        c : *numpy.ndarray*
            Derived constants and controls.
        t : *float*
            Time at which the values are calculated.

        Returns
        -------
        nonlinearities : *numpy.ndarray*
            Nonlinearities.
        """

        # extract frequently used variables
        divisor = self.params['J'] / self.params['x_0']**2

        # return nonlinearities
        return - self.params['gamma'] + 2.0j * self.params['g_0'] * np.real(modes[1::2]) / divisor