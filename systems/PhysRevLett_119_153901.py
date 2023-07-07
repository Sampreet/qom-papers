#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the N-cell OM array system in Phys. Rev. Lett. **119**, 153901 (2017)."""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.0'
__created__ = '2021-08-09'
__updated__ = '2023-07-07'

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevLett_119_153901(BaseSystem):
    r"""Class to simulate the N-cell OM array system in Phys. Rev. Lett. **119**, 153901 (2017).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        N               (*int*) number of cells :math:`N`. Default is :math:`400`.
        Gamma_m         (*float*) normalized normalized mechanical dampling rate :math:`\Gamma_{m}` in Hertz. Default is :math:`190.0 / 9.5 \times 10^{9}`.
        g_0             (*float*) normalized optomechanical coupling strength :math:`g_{0}` in Hertz. Default is :math:`292 \times 10^{3} / 9.5 \times 10^{9}`.
        J               (*float*) normalized photon hopping strength :math:`J` in Hertz. Default is :math:`- 95 \times 10^{9} / 9.5 \times 10^{9}`.
        kappa           (*float*) normalized optical decay rate :math:`\kappa`. Default is :math:`3.8 \times 10^{6} / 9.5 \times 10^{9}`.
        Omega_m         (*float*) normalized mechanical freqency :math:`\Omega_{m}`. Default is :math:`1.0`.
        x_d             (*float*) envelope width :math:`x_{d}`. Default is :math:`40.0`.
        t_alphas        (*str*) type of initial occupancies `t_alphas`. Default is "Gaussian".
        ============    ========================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """
    
    # default parameters of the system
    system_defaults = {
        'N'         : 400,
        'Gamma_m'   : 190 / 9.5e9,
        'g_0'       : 292e3 / 9.5e9,
        'J'         : - 95e9 / 9.5e9,
        'kappa'     : 3.8e6 / 9.5e9,
        'Omega_m'   : 1.0,
        'x_d'       : 40,
        't_alphas'  : 'Gaussian'
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_119_153901."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevLett_119_153901',
            desc='Array System in Phys. Rev. Lett. 119, 153901',
            num_modes=2 * params.get('N', self.system_defaults['N']),
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
        divisor = np.abs(self.params['J']) / self.params['x_d']**2

        # return rates
        return (1.0j * self.params['g_0'] * np.conjugate(modes[::2]) * modes[::2] - (self.params['Gamma_m'] + 1.0j * self.params['Omega_m']) * modes[1::2]) / divisor

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
        N       = self.params['N']
        Omega_m = self.params['Omega_m']
        x_d     = self.params['x_d']
        gamma   = 2 * self.params['g_0']**2 * Omega_m / (self.params['Gamma_m']**2 + Omega_m**2)
        R       = 1 / x_d / np.sqrt(gamma / np.abs(self.params['J']))
        c_1     = 1
        c_2     = 2.5
        c_3     = 5
 
        # set initial optical amplitudes
        def gauss(X, mu, sigma):
            return np.exp(- np.power(X - mu, 2) / (2 * np.power(sigma, 2)))
        xs  = np.linspace(-1, 1, N)
        ys  = xs / x_d / (xs[1] - xs[0])
        # default initial values of optical modes
        default = 1e3 - 2.5e3 * gauss(xs, 0, 0.1) + 7.5e3 * gauss(xs, 0, 0.05)

        # initial values of the modes
        iv_modes        = np.zeros(2 * N, dtype=np.complex_)
        iv_modes[::2]   = {
            'Gaussian'  : R * (c_1 + c_2 * np.exp(- c_3 * ys**2)),
            'Lorentzian': R * c_2 / (c_1 + c_3 * ys**2)**2,
            'square'    : np.array([c_2 * R if np.abs(x) <= 0.15 else R for x in xs]),
            'soliton'   : R * (- 6.0 / (1.0 - 2.0 * np.cosh(ys)) - 1.0),
            'default'   : default
        }.get(self.params['t_alphas'], default)

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
        divisor = np.abs(J) / self.params['x_d']**2

        # initialize mode rates
        mode_rates = np.zeros_like(modes, dtype=np.complex_)

        # update rates for optical modes
        for i in range(len(alphas)):
            mode_rates[2 * i] = (- 1.0j * J / 2.0 * (alphas[i - 1] if i > 0 else 0.0) + (- self.params['kappa'] / 2.0 + 2.0j * g_0 * np.real(betas[i]) + 1.0j * J) * alphas[i] - 1.0j * J / 2.0 * (alphas[i + 1] if i < len(alphas) - 1 else 0.0)) / divisor
        # update rates for mechanical modes
        mode_rates[1::2] = (1.0j * g_0 * np.conjugate(alphas) * alphas - (self.params['Gamma_m'] + 1.0j * self.params['Omega_m']) * betas) / divisor

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
        divisor = np.abs(J) / self.params['x_d']**2

        # return coefficients
        return np.array([0.0, 0.0, - 1.0j * J / 2.0 * 1.0**2 / divisor], dtype=np.complex_)

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
        divisor = np.abs(self.params['J']) / self.params['x_d']**2

        # return nonlinearities
        return - self.params['kappa'] / 2.0 + 2.0j * self.params['g_0'] * np.real(modes[1::2]) / divisor