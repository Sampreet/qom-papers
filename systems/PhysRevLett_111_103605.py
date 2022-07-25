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
__created__ = '2021-05-18'
__updated__ = '2022-07-25'
__version__ = '0.8.5'

# dependencies
import numpy as np

# qom modules
from qom.systems import DODMSystem

class PhysRevLett_111_103605(DODMSystem):
    r"""Class to simulate the coupled QOM system in Phys. Rev. Lett. **111**, 103605 (2017).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        E_norm      (float) normalized amplitude of the laser :math:`E / \kappa`. Default is :math:`320.0`.
        g_norm      (float) normalized optomechanical coupling strength :math:`g / \omega_{1}`. Default is :math:`0.005`.
        gamma_norm  (float) normalized mechanical damping rate :math:`\gamma / \omega_{1}`. Default is :math:`0.005`.
        kappa_norm  (float) normalized optical decay rate :math:`\kappa / \omega_{1}`. Default is :math:`0.15`.
        mu_norm     (float) normalized coupling strength :math:`\mu / \omega_{1}`. Default is :math:`0.2`.
        n_b         (float) quanta of thermal phonons :math:`n_{b}`. Default is :math:`0.0`.
        omega_2_norm(float) normalized frequency of the second mechanical oscillator :math:`\omega_{2} / \omega_{1}`. Default is :math:`1.005`.
        ========    ============================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default system parameters
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
        super().__init__(params, cb_update=cb_update)

        # set attributes
        self.code = 'physrevlett_111_103605'
        self.name = 'Coupled System in Phys. Rev. Lett. 111, 103605'

        # update parameters
        self.params = dict()
        for key in self.system_defaults:
            self.params[key] = params.get(key, self.system_defaults[key])

        # matrices
        self.A = None

    def get_A(self, modes, params, t):
        """Function to obtain the drift matrix.

        Parameters
        ----------
        modes : list
            Values of the modes.
        params : list
            Constants parameters.
        t : float
            Time at which the drift matrix is calculated.
        
        Returns
        -------
        A : list
            Drift matrix.
        """

        # extract frequently used variables
        E, g, gamma, kappa, mu, omega_1, omega_2 = params
        omegas = [omega_1, omega_2]
        alphas = [modes[0], modes[2]]
        betas = [modes[1], modes[3]]
        tau = 2 * np.pi / omega_1

        # effective values
        Deltas = list()
        Gs = list()
        for i in range(2):
            Deltas.append(omegas[i] + 2 * g * np.real(betas[i]))
            Gs.append(g * alphas[i])

        # initialize drift matrix
        if self.A is None or np.shape(self.A) != (8, 8):
            self.A = np.zeros([8, 8], dtype=np.float_)
        for i in range(2):
            # optical position quadrature
            self.A[4 * i + 0][4 * i + 0] = - kappa
            self.A[4 * i + 0][4 * i + 1] = - Deltas[i]
            self.A[4 * i + 0][4 * i + 2] = - 2 * np.imag(Gs[i])
            # optical momentum quadrature
            self.A[4 * i + 1][4 * i + 0] = Deltas[i]
            self.A[4 * i + 1][4 * i + 1] = - kappa
            self.A[4 * i + 1][4 * i + 2] = 2 * np.real(Gs[i])
            # mechanical position quadrature
            self.A[4 * i + 2][4 * i + 2] = - gamma
            self.A[4 * i + 2][4 * i + 3] = omegas[i]
            self.A[4 * i + 2][4 * (1 - i) + 3] = - mu
            # mechanical momentum quadrature
            self.A[4 * i + 3][4 * i + 0] = 2 * np.real(Gs[i])
            self.A[4 * i + 3][4 * i + 1] = 2 * np.imag(Gs[i])
            self.A[4 * i + 3][4 * i + 2] = - omegas[i]
            self.A[4 * i + 3][4 * i + 3] = - gamma
            self.A[4 * i + 3][4 * (1 - i) + 2] = mu

        return self.A * tau

    def get_ivc(self):
        r"""Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.
            First element contains the optical mode of first cavity.
            Next element contains the mechanical mode of first cavity.
            First element contains the optical mode of second cavity.
            Next element contains the mechanical mode of second cavity.
            Next :math:`4 \times 4^{2}` elements contain the correlations.

        c : list
            Constants of the IVP.
            First :math:`4 \times 4^{2}` elements contain the noise matrix.
            The elements contain the system parameters ``params`` in the following order:
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
        tau     = 2 * np.pi / omega_1
 
        # initial mode values as 1D list
        modes_0 = np.zeros(self.num_modes, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros(dim, dtype=np.float_)
        for i in range(2):
            corrs_0[4* i + 0][4* i + 0] = 0.5
            corrs_0[4* i + 1][4* i + 1] = 0.5
            corrs_0[4* i + 2][4* i + 2] = (n_b + 0.5)
            corrs_0[4* i + 3][4* i + 3] = (n_b + 0.5)

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros(dim, dtype=np.float_)
        for i in range(2):
            D[4* i + 0][4* i + 0] = kappa
            D[4* i + 1][4* i + 1] = kappa
            D[4* i + 2][4* i + 2] = gamma * (2 * n_b + 1)
            D[4* i + 3][4* i + 3] = gamma * (2 * n_b + 1)
        D = D * tau
        
        # constant parameters
        params = [E, g, gamma, kappa, mu, omega_1, omega_2]

        # all constants
        c = D.flatten().tolist() + params

        return iv, c

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
        E, g, gamma, kappa, mu, omega_1, omega_2 = params
        omegas = [omega_1, omega_2]
        alphas = [modes[0], modes[2]]
        betas = [modes[1], modes[3]]
        tau = 2 * np.pi / omega_1

        # effective values
        Deltas = list()
        Gs = list()
        for i in range(2):
            Deltas.append(omegas[i] + 2 * g * np.real(betas[i]))
            Gs.append(g * alphas[i])

        # calculate rates
        dalpha_dts = list()
        dbeta_dts = list()
        for i in range(2):
            dalpha_dts.append(((- kappa + 1j * Deltas[i]) * alphas[i] + E))
            dbeta_dts.append((1j * Gs[i] * np.conjugate(alphas[i]) + (- gamma - 1j * omegas[i]) * betas[i] + 1j * mu * betas[1 - i]))

        # rearrange rates
        mode_rates = [dalpha_dts[0] * tau, dbeta_dts[0] * tau, dalpha_dts[1] * tau, dbeta_dts[1] * tau]

        return mode_rates