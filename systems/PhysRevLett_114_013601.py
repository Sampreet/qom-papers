#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the chaotic OM system in Phys. Rev. Lett. **114**, 013601 (2015)."""

__authors__ = ["Sampreet Kalita"]
__toolbox__ = "qom-v1.1.0"
__created__ = "2021-07-27"
__updated__ = "2025-03-11"
__all__     = ['PhysRevLett_114_013601']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevLett_114_013601(BaseSystem):
    r"""Class to simulate the chaotic OM system in Phys. Rev. Lett. **114**, 013601 (2015).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        Delta_norm  (*float*) normalized laser detuning :math:`\Delta / \Omega`. Default is :math:`0.0`.
        Gamma_norm  (*float*) normalized mechanical damping rate :math:`\Gamma / \Omega`. Default is :math:`10^{-3}`.
        kappa_norm  (*float*) normalized optical decay rate :math:`\kappa / \Omega`. Default is :math:`1.0`.
        P           (*float*) pump parameter :math:`P`. Default is `1.4`.
        ========    ============================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is a float and ``reset`` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'Delta_norm'    : 0.0,
        'Gamma_norm'    : 1e-3,
        'kappa_norm'    : 1.0,
        'P'             : 1.4
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_114_013601."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevLett_114_013601',
            desc="Chaotic System in Phys. Rev. Lett. 114, 013601",
            num_modes=2,
            cb_update=cb_update
        )

    def get_A(self, modes, c, t):
        """Method to obtain the drift matrix.

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
        A : *numpy.ndarray*
            Drift matrix.
        """

        # extract frequently used variables
        alpha, beta = modes

        # optical mode
        self.A[0][0]    = - self.params['kappa_norm'] / 2.0
        self.A[0][1]    = - self.params['Delta_norm'] + 2.0 * np.real(beta)
        self.A[0][2]    = 2.0 * np.imag(alpha)
        self.A[1][0]    = self.params['Delta_norm'] - 2.0 * np.real(beta)
        self.A[1][1]    = - self.params['kappa_norm'] / 2.0
        self.A[1][2]    = - 2.0 * np.real(alpha)
        # mechanical mode
        self.A[2][2]    = - self.params['Gamma_norm'] / 2.0
        self.A[2][3]    = 1.0
        self.A[3][0]    = - self.params['P'] * np.real(alpha)
        self.A[3][1]    = - self.params['P'] * np.imag(alpha)
        self.A[3][2]    = - 1.0
        self.A[3][3]    = - self.params['Gamma_norm'] / 2.0

        return self.A
    
    def get_coeffs_N_o(self, c):
        """Method to obtain coefficients of the polynomial in mean optical occupancy.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns 
        -------
        coeffs : *numpy.ndarray*
            Coefficients of the polynomial in mean optical occupancy.
        """
        
        # drive amplitude
        A_l_norm = - 0.5j
        # Coefficient of the mean optical occupancies
        C = 4.0 * self.params['P'] / (self.params['Gamma_norm']**2 + 4.0)
        
        # get coefficients
        coeffs      = np.zeros(2 * self.num_modes, dtype=np.float64)
        coeffs[0]   = 4.0 * C**2
        coeffs[1]   = 8.0 * C * self.params['Delta_norm']
        coeffs[2]   = 4.0 * self.params['Delta_norm']**2 + self.params['kappa_norm']**2
        coeffs[3]   = - 4.0 * np.real(np.conjugate(A_l_norm) * A_l_norm)

        return coeffs
    
    def get_D(self, modes, corrs, c, t):
        """Method to obtain the noise matrix.
        
        Parameters
        ----------
        modes : *numpy.ndarray*
            Classical modes.
        corrs : *numpy.ndarray*
            Quantum correlations.
        c : *numpy.ndarray*
            Derived constants and controls.
        t : *float*
            Time at which the values are calculated.
        
        Returns
        -------
        D : *numpy.ndarray*
            Noise matrix.
        """

        # update noise matrix
        self.D[0][0]    = self.params['kappa_norm']
        self.D[1][1]    = self.params['kappa_norm']
        self.D[2][2]    = self.params['Gamma_norm']
        self.D[3][3]    = self.params['Gamma_norm']

        return self.D

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

        # initial values of the modes
        iv_modes = np.zeros(self.num_modes, dtype=np.complex128)

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float64)
        iv_corrs[0][0]  = 0.5
        iv_corrs[1][1]  = 0.5
        iv_corrs[2][2]  = 0.5
        iv_corrs[3][3]  = 0.5

        return iv_modes, iv_corrs, None

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
        alpha, beta = modes

        # calculate mode rates
        dalpha_dt = (1.0j * self.params['Delta_norm'] - self.params['kappa_norm'] / 2.0) * alpha - 2.0j * alpha * np.real(beta) - 1.0j / 2.0
        dbeta_dt = (- 1.0j - self.params['Gamma_norm'] / 2.0) * beta - 1.0j * self.params['P'] / 2.0 * np.conjugate(alpha) * alpha

        return np.array([dalpha_dt, dbeta_dt], dtype=np.complex128)

    def get_modes_steady_state(self, c):
        """Method to obtain the steady state modes.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns 
        -------
        Modes : *numpy.ndarray*
            Steady state modes.
        """

        # get mean optical occupancy
        N_o = self.get_mean_optical_occupancies()[0]

        # calculate mechanical mode position
        beta_real = - 2.0 * self.params['P'] / (self.params['Gamma_norm']**2 + 4.0) * N_o

        # calculate steady state modes
        alpha = - 1.0j / (self.params['kappa_norm'] - 2.0j * (self.params['Delta_norm'] - 2.0 * beta_real))
        beta = - self.params['P'] * N_o * (2.0 + 1.0j * self.params['Gamma_norm']) / (self.params['Gamma_norm']**2 + 4.0)
        
        return np.array([[alpha, beta]], dtype=np.complex128)