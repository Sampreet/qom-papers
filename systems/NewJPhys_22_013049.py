#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the weakly dissipative OM system in New J. Phys. **22**, 013049 (2020)."""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.2'
__created__ = '2021-07-27'
__updated__ = '2024-06-23'
__all__     = ['NewJPhys_22_013049']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class NewJPhys_22_013049(BaseSystem):
    r"""Class to simulate the weakly dissipative OM system in New J. Phys. **22**, 013049 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        Delta_norm  (*float*) normalized laser detuning :math:`\Delta / \Omega_{m}`. Default is :math:`0.0`.
        gamma_norm  (*float*) normalized mechanical damping rate :math:`\gamma / \Omega_{m}`. Default is :math:`10^{-4}`.
        kappa_norm  (*float*) normalized optical decay rate :math:`\kappa / \Omega_{m}`. Default is :math:`0.1`.
        P           (*float*) normalized coupling strength :math:`P`. Default is :math:`0.1`.
        ========    ============================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is an integer and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'Delta_norm': 0.0,
        'gamma_norm': 1e-4,
        'kappa_norm': 0.1,
        'P'         : 0.1
    }

    def __init__(self, params={}, cb_update=None):
        """Class constructor for NewJPhys_22_013049."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='NewJPhys_22_013049',
            desc='Weakly Dissipative System in New J. Phys. 22, 013049',
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
        alpha, beta = [modes[0], modes[1]]
        
        # update drift matrix
        # optical mode
        self.A[0][0] = - self.params['kappa_norm'] / 2
        self.A[0][1] = - self.params['Delta_norm'] - 2 * np.real(beta)
        self.A[0][2] = - 2.0 * np.imag(alpha)
        self.A[1][0] = self.params['Delta_norm'] + 2 * np.real(beta)
        self.A[1][1] = - self.params['kappa_norm'] / 2
        self.A[1][2] = 2.0 * np.real(alpha)
        # mechanical mode
        self.A[2][2] = - self.params['gamma_norm'] / 2
        self.A[2][3] = 1.0
        self.A[3][0] = self.params['P'] * np.real(alpha)
        self.A[3][1] = self.params['P'] * np.imag(alpha)
        self.A[3][2] = - 1.0
        self.A[3][3] = - self.params['gamma_norm'] / 2

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

        # frequently used variables
        A_l_norm, Delta_norm, kappa_norm, C = self.get_params_steady_state(
            c=c
        )
        
        # get coefficients
        coeffs      = np.zeros(2 * self.num_modes, dtype=np.float_)
        coeffs[0]   = 4.0 * C**2
        coeffs[1]   = 8.0 * C * Delta_norm
        coeffs[2]   = 4.0 * Delta_norm**2 + kappa_norm**2
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
        self.D[0][0] = self.params['kappa_norm'] / 2.0
        self.D[1][1] = self.params['kappa_norm'] / 2.0
        self.D[2][2] = self.params['gamma_norm'] / 2.0
        self.D[3][3] = self.params['gamma_norm'] / 2.0
        
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
        iv_modes = np.zeros(self.num_modes, dtype=np.complex_)

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
            Rates of change of each mode.
        """
        
        # extract frequently used variables
        alpha, beta = [modes[0], modes[1]]

        # initialize lists
        dalpha_dt = (1.0j * self.params['Delta_norm'] - self.params['kappa_norm'] / 2.0) * alpha + 1.0j * alpha * (beta + np.conjugate(beta)) + 0.5
        dbeta_dt = (- 1.0j - self.params['gamma_norm'] / 2.0) * beta + 1.0j * self.params['P'] / 2.0 * np.conjugate(alpha) * alpha

        # rearrange per system
        return np.array([dalpha_dt, dbeta_dt], dtype=np.complex_)

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

        # frequently used variables
        _, Delta_norm, kappa_norm, _ = self.get_params_steady_state(
            c=c
        )

        # initialize lists
        Modes = list()

        # get mean optical occupancies
        N_os = self.get_mean_optical_occupancies()

        # for each mean optical occupancy
        for N_o in N_os:
            # mechanical mode position
            beta_real = 2.0 * self.params['P'] / (self.params['gamma_norm']**2 + 4.0) * N_o

            # calculate mode amplitudes
            alpha = 1.0 / (kappa_norm - 2.0j * (Delta_norm + 2.0 * beta_real))
            beta = self.params['P'] * N_o * (2.0 + 1.0j * self.params['gamma_norm']) / (self.params['gamma_norm']**2 + 4.0)

            # append to list
            Modes.append([alpha, beta])

        return np.array(Modes, dtype=np.complex_)

    def get_params_steady_state(self, c):
        r"""Method to obtain the parameters required to calculate the optical steady states.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns
        -------
        A_l_norm : *float*
            Normalized amplitude of the laser.
        Delta_0_norm : *float*
            Normalized detuning of the laser.
        kappa_norm : *float*
            Normalized optical decay rate.
        C : *float*
            Coefficient of :math:`| \alpha |^{2}`.
        """

        # Coefficient of the mean optical occupancies
        C = 4.0 * self.params['P'] / (self.params['gamma_norm']**2 + 4.0)

        return 0.5, self.params['Delta_norm'], self.params['kappa_norm'], C