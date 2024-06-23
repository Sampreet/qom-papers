#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the coupled QOM system in Phys. Rev. Lett. **111**, 103605 (2017).

# Measures of Quantum Synchronization in Continuous Variable Systems

==============  =====================================================
Author          Affiliation
==============  =====================================================
A. Mari         NEST, Scuola Normale Superiore and Istituto Nanoscienze-CNR, I-56127 Pisa, Italy
A. Farace       NEST, Scuola Normale Superiore and Istituto Nanoscienze-CNR, I-56127 Pisa, Italy
N. Didier   Universite de Sherbrooke, Sherbrooke, Quebec JIK 2R1, Canada
                McGrill University, Montreal, Quebec H3A 2T8, Canada
V. Giovannetti  NEST, Scuola Normale Superiore and Istituto Nanoscienze-CNR, I-56127 Pisa, Italy
R. Fazio        NEST, Scuola Normale Superiore and Istituto Nanoscienze-CNR, I-56127 Pisa, Italy
==============  =====================================================

## Objectives
* Introduce quantum synchronization measures for coupled continuous variable systems.
* Study of dynamics of coupled optomechanical systems.

## Novelty
* Extrapolation of notions of complete and phase synchronization in classical models.
* Bounds of quantum mechanics on the achievable level of synchronization.

## Assumptions
* Equal driving and coupling strengths.
* Variable mechanical frequencies.
* Linearized description.

## Results
* Values of complete and phase synchronization are consistent with the fundamental limit.
* No functional relationship between discord and synchronization observed.
"""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.2'
__created__ = '2021-05-18'
__updated__ = '2024-06-23'
__all__     = ['PhysRevLett_111_103605']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevLett_111_103605(BaseSystem):
    r"""Class to simulate the coupled QOM system in Phys. Rev. Lett. **111**, 103605 (2017).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ============================================================
        key             meaning
        ============    ============================================================
        E_norm          (*float*) normalized amplitude of the laser :math:`E / \kappa`. Default is :math:`320.0`.
        g_norm          (*float*) normalized optomechanical coupling strength :math:`g / \omega_{1}`. Default is :math:`0.005`.
        gamma_norm      (*float*) normalized mechanical damping rate :math:`\gamma / \omega_{1}`. Default is :math:`0.005`.
        kappa_norm      (*float*) normalized optical decay rate :math:`\kappa / \omega_{1}`. Default is :math:`0.15`.
        mu_norm         (*float*) normalized coupling strength :math:`\mu / \omega_{1}`. Default is :math:`0.2`.
        n_b             (*float*) quanta of thermal phonons :math:`n_{b}`. Default is :math:`0.0`.
        omega_2_norm    (*float*) normalized frequency of the second mechanical oscillator :math:`\omega_{2} / \omega_{1}`. Default is :math:`1.005`.
        ============    ============================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'E_norm'        : 320.0,
        'g_norm'        : 0.005,
        'gamma_norm'    : 0.005,
        'kappa_norm'    : 0.15,
        'mu_norm'       : 0.02,
        'n_b'           : 0.0,
        'omega_2_norm'  : 1.005
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_111_103605."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevLett_111_103605',
            desc='Coupled System in Phys. Rev. Lett. 111, 103605',
            num_modes=4,
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
        g, gamma, kappa, mu = c[1:5]
        omegas              = c[5:7]
        alphas              = [modes[0], modes[2]]
        betas               = [modes[1], modes[3]]

        # effective values
        Deltas  = [omegas[i] + 2.0 * g * np.real(betas[i]) for i in range(2)]
        Gs      = [g * alphas[i] for i in range(2)]

        # update drift matrix
        for i in range(2):
            # optical position quadrature
            self.A[4 * i + 0][4 * i + 0]        = - kappa
            self.A[4 * i + 0][4 * i + 1]        = - Deltas[i]
            self.A[4 * i + 0][4 * i + 2]        = - 2.0 * np.imag(Gs[i])
            # optical momentum quadrature
            self.A[4 * i + 1][4 * i + 0]        = Deltas[i]
            self.A[4 * i + 1][4 * i + 1]        = - kappa
            self.A[4 * i + 1][4 * i + 2]        = 2.0 * np.real(Gs[i])
            # mechanical position quadrature
            self.A[4 * i + 2][4 * i + 2]        = - gamma
            self.A[4 * i + 2][4 * i + 3]        = omegas[i]
            self.A[4 * i + 2][4 * (1 - i) + 3]  = - mu
            # mechanical momentum quadrature
            self.A[4 * i + 3][4 * i + 0]        = 2.0 * np.real(Gs[i])
            self.A[4 * i + 3][4 * i + 1]        = 2.0 * np.imag(Gs[i])
            self.A[4 * i + 3][4 * i + 2]        = - omegas[i]
            self.A[4 * i + 3][4 * i + 3]        = - gamma
            self.A[4 * i + 3][4 * (1 - i) + 2]  = mu
        # normalize
        self.A *= 2.0 * np.pi / omegas[0]

        return self.A
    
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

        # extract frequently used variables
        gamma, kappa    = c[2:4]
        omega_1         = c[5]

        # update noise matrix
        for i in range(2):
            self.D[4* i + 0][4* i + 0]  = kappa
            self.D[4* i + 1][4* i + 1]  = kappa
            self.D[4* i + 2][4* i + 2]  = gamma * (2.0 * self.params['n_b'] + 1.0)
            self.D[4* i + 3][4* i + 3]  = gamma * (2.0 * self.params['n_b'] + 1.0)
        # normalize
        self.D *= 2.0 * np.pi / omega_1
        
        return self.D

    def get_ivc(self):
        r"""Method to obtain the initial values of the modes, correlations and derived constants and controls.
        
        Returns
        -------
        iv_modes : *numpy.ndarray*
            Initial values of the classical modes.
        iv_corrs : *numpy.ndarray*
            Initial values of the quantum correlations.
        c : *numpy.ndarray*
            Derived constants and controls in the following order:
            ========    =============================================
            index       parameter
            ========    =============================================
            0           laser amplitude :math:`E`.
            1           optomechanical coupling strength :math:`g`.
            2           mechanical damping rate :math:`\gamma`.
            3           optical decay rate :math:`\kappa`.
            4           coupling strength :math:`\mu`.
            5           first mechanical frequency :math:`\omega_{1}`.
            6           second mechanical frequency :math:`\omega_{2}`.
            ========    =============================================
        """

        # extract frequently used variables
        omega_1 = 1.0
        E       = self.params['E_norm'] * self.params['kappa_norm'] * omega_1
        g       = self.params['g_norm'] * omega_1
        gamma   = self.params['gamma_norm'] * omega_1
        kappa   = self.params['kappa_norm'] * omega_1
        mu      = self.params['mu_norm'] * omega_1
        n_b     = self.params['n_b']
        omega_2 = self.params['omega_2_norm'] * omega_1
        dim     = (2 * self.num_modes, 2 * self.num_modes)
 
        # initial values of the modes
        iv_modes = np.zeros(self.num_modes, dtype=np.complex_)

        # initial quadrature correlations
        iv_corrs = np.zeros(dim, dtype=np.float_)
        for i in range(2):
            iv_corrs[4* i + 0][4* i + 0]    = 0.5
            iv_corrs[4* i + 1][4* i + 1]    = 0.5
            iv_corrs[4* i + 2][4* i + 2]    = n_b + 0.5
            iv_corrs[4* i + 3][4* i + 3]    = n_b + 0.5
        
        # derived constants
        c = np.array([E, g, gamma, kappa, mu, omega_1, omega_2], dtype=np.float_)

        return iv_modes, iv_corrs, c

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
        E, g, gamma, kappa, mu  = c[0:5]
        omegas                  = c[5:7]
        alphas                  = [modes[0], modes[2]]
        betas                   = [modes[1], modes[3]]

        # effective values
        Deltas  = [omegas[i] + 2.0 * g * np.real(betas[i]) for i in range(2)]
        Gs      = [g * alphas[i] for i in range(2)]

        # calculate rates
        dalpha_dts  = [(- kappa + 1.0j * Deltas[i]) * alphas[i] + E for i in range(2)]
        dbeta_dts   = [1.0j * Gs[i] * np.conjugate(alphas[i]) + (- gamma - 1.0j * omegas[i]) * betas[i] + 1.0j * mu * betas[1 - i] for i in range(2)]

        # rearrange rates, normalize and return
        return np.array([dalpha_dts[0], dbeta_dts[0], dalpha_dts[1], dbeta_dts[1]], dtype=np.complex_) * 2.0 * np.pi / omegas[0]