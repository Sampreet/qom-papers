#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the systems described in Phys. Rev. Lett. **111**, 103605 (2017).

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
__updated__ = '2021-08-30'
__version__ = '0.8.0'

# dependencies
import numpy as np

# qom modules
from qom.systems import DODMSystem

class PhysRevLett111_103605(DODMSystem):
    """Class to simulate the coupled QOM system in Phys. Rev. Lett. **111**, 103605 (2017).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett111_103605."""
        
        # initialize super class
        super().__init__(params, cb_update=cb_update)

        # set attributes
        self.code = 'phys_rev_lett_111_103605'
        self.name = 'Coupled QOM System in Phys. Rev. Lett. 111, 103605'
        # default parameters
        self.params = {
            'E': params.get('E', 320),
            'g': params.get('g', 0.005),
            'gamma': params.get('gamma', 0.005),
            'kappa': params.get('kappa', 0.15),
            'mu': params.get('mu', 0.02),
            'n_b': params.get('n_b', 0),
            'omega_1': params.get('omega_1', 1),
            'omega_2': params.get('omega_2', 1.005)
        }

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
        Delta_0s    = [params[0], params[1]]
        g           = params[3]
        gamma       = params[4]
        kappa       = params[5]
        mu          = params[6]
        omegas      = [params[7], params[8]]
        alphas      = [modes[0], modes[2]]
        betas       = [modes[1], modes[3]]
        tau         = 2 * np.pi / omegas[0]

        # effective values
        Deltas = list()
        Gs = list()
        for i in range(2):
            Deltas.append(Delta_0s[i] + 2 * g * np.real(betas[i]))
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
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.

        c : list
            Constant parameters of the system.
            First two elements contain the laser detunings.
            Next element contains the laser amplitude.
            Next element contains the optomechanical coupling strength.
            Next element contains the mechanical decay rate.
            Next element contains the optical decay rate.
            Next element contains the inter-cavity coupling constant.
            Next two elements contain the mechanical frequencies.
        """

        # extract frequently used variables
        E           = self.params['E']
        g           = self.params['g']
        gamma       = self.params['gamma']
        kappa       = self.params['kappa']
        mu          = self.params['mu']
        n_b         = self.params['n_b']
        omega_1     = self.params['omega_1']
        omega_2     = self.params['omega_2']
        dim         = (2 * self.num_modes, 2 * self.num_modes)
        tau         = 2 * np.pi / omega_1
 
        # initial mode values as 1D list
        modes_0 = np.zeros(self.num_modes, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros(dim, dtype=np.float_)
        for i in range(2):
            corrs_0[4* i + 0][4* i + 0] = 1/2
            corrs_0[4* i + 1][4* i + 1] = 1/2
            corrs_0[4* i + 2][4* i + 2] = (n_b + 1/2)
            corrs_0[4* i + 3][4* i + 3] = (n_b + 1/2)

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
        params = [omega_1, omega_2] + \
            [E, g, gamma, kappa, mu] + \
            [omega_1, omega_2]

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
        Delta_0s    = [params[0], params[1]]
        E           = params[2]
        g           = params[3]
        gamma       = params[4]
        kappa       = params[5]
        mu          = params[6]
        omegas      = [params[7], params[8]]
        alphas      = [modes[0], modes[2]]
        betas       = [modes[1], modes[3]]
        tau         = 2 * np.pi / omegas[0]

        # effective values
        Deltas = list()
        Gs = list()
        for i in range(2):
            Deltas.append(Delta_0s[i] + 2 * g * np.real(betas[i]))
            Gs.append(g * alphas[i])

        # calculate rates
        dalpha_dts = list()
        dbeta_dts = list()
        for i in range(2):
            dalpha_dts.append(((- kappa + 1j * Deltas[i]) * alphas[i] + E * kappa))
            dbeta_dts.append((1j * Gs[i] * np.conjugate(alphas[i]) + (- gamma - 1j * omegas[i]) * betas[i] + 1j * mu * betas[1 - i]))

        # rearrange rates
        mode_rates = [dalpha_dts[0] * tau, dbeta_dts[0] * tau, dalpha_dts[1] * tau, dbeta_dts[1] * tau]

        return mode_rates