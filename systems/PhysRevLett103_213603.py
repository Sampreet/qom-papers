#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the systems described in Phys. Rev. Lett. **103**, 213603 (2009).

# Gently Modulating Optomechanical Systems

==============  =====================================================
Author          Affiliation
==============  =====================================================
A. Mari         Institute of Physics and Astronomy, Univerity of Potsdam, D-14476 Postdam, Germany
J. Eisert       Institute of Advanced Study Berlin, D-14193 Berlin, Germany
==============  =====================================================

## Objectives
* Modulated entanglement dynamics between optics and mechanics.

## Novelty
* Driving with mildly amplitude-modulated light.
* No classical feedback or squeezed-light driving.

## Assumptions
* Linearized description of the modes.

## Results
* High degrees of squeezing below the vacuum noise level.
"""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-05-16'
__updated__ = '2021-08-30'
__version__ = '0.8.0'

# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import SOSMSystem

class PhysRevLett103_213603(SOSMSystem):
    """Class to simulate the gently modulated QOM system in Phys. Rev. Lett. **103**, 213603 (2009).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett103_213603."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'phys_rev_lett_103_213603'
        self.name = 'Gently Modulated QOM System in Phys. Rev. Lett. 103, 213603'
        # default parameters
        self.params = {
            'F': params.get('F', 1.4e4),
            'L': params.get('L', 25e-3),
            'lambda_l': params.get('lambda_l', 1064e-9),
            'm': params.get('m', 150e-12),
            'omega_m': params.get('omega_m', 2e6 * np.pi),
            'P_0': params.get('P_0', 10e-3),
            'P_1': params.get('P_1', 2e-3),
            'Q': params.get('Q', 1e6),
            'T': params.get('T', 0.1)
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
        Delta_0 = params[0]
        G_0     = params[3]
        gamma_m = params[4]
        kappa   = params[5]
        omega_m = params[6]
        Omega   = params[7]
        alpha   = modes[0]
        beta    = modes[1]
        tau     = 2 * np.pi / Omega

        # effective values
        Delta = Delta_0 - np.sqrt(2) * G_0 * np.real(beta)
        G = np.sqrt(2) * G_0 * alpha

        # initialize drift matrix
        if self.A is None or np.shape(self.A) != (4, 4):
            self.A = np.zeros([4, 4], dtype=np.float_)
        # optical position quadrature
        self.A[0][0] = - kappa 
        self.A[0][1] = Delta 
        self.A[0][2] = - np.imag(G) 
        # optical momentum quadrature
        self.A[1][0] = - Delta
        self.A[1][1] = - kappa
        self.A[1][2] = np.real(G)
        # mechanical position quadrature
        self.A[2][3] = omega_m
        # mechanical momentum quadrature
        self.A[3][0] = np.real(G)
        self.A[3][1] = np.imag(G)
        self.A[3][2] = - omega_m
        self.A[3][3] = - gamma_m

        return self.A * tau

    def get_ivc(self):
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.

        c : list
            Constant parameters of the system.
            First element contains the laser detuning.
            Next two elements contain the laser base and modulation amplitudes.
            Next element contains the optomechanical coupling constant.
            Next two elements contain the mechanical and optical decay rates.
            Next two elements contain the mechanical and modulation frequencies.
        """
        
        # extract frequently used variables
        F       = self.params['F']
        L       = self.params['L']
        lambda_l= self.params['lambda_l']
        m       = self.params['m']
        omega_m = self.params['omega_m']
        P_0     = self.params['P_0']
        P_1     = self.params['P_1']
        Q       = self.params['Q']
        T       = self.params['T']

        # laser detuning
        Delta_0 = omega_m
        # mechanical decay rate
        gamma_m = omega_m / Q
        # optical decay rate
        kappa = np.pi * sc.c / (2 * F * L)
        # laser frequency
        omega_l = 2 * np.pi * sc.c / lambda_l
        # cavity frequency
        omega_c = Delta_0 + omega_l
        # coupling strength
        G_0 = np.sqrt(sc.hbar / (m * omega_m)) * omega_c / L
        # modulation frequency
        Omega = 2 * omega_m
        # time frame
        tau = 2 * np.pi / Omega

        # thermal phonon number
        if T != 0.0 or T != 0:
            n_th = 1 / (np.exp(sc.hbar * omega_m / (sc.k * T)) - 1)
        else:
            n_th = 0
 
        # initial mode values as 1D list
        modes_0 = np.zeros(2, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros([4, 4], dtype=np.float_)
        corrs_0[0][0] = 1/2 
        corrs_0[1][1] = 1/2
        corrs_0[2][2] = (n_th + 1/2)
        corrs_0[3][3] = (n_th + 1/2)

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros([4, 4], dtype=np.float_)
        D[0][0] = kappa * tau
        D[1][1] = kappa * tau
        D[3][3] = gamma_m * (2 * n_th + 1) * tau

        # laser amplitudes
        E_0 = np.sqrt(2 * kappa * P_0 / (sc.hbar * omega_l)) 
        E_1 = np.sqrt(2 * kappa * P_1 / (sc.hbar * omega_l))
        
        # constant parameters
        params = [Delta_0] + \
            [E_0, E_1] + \
            [G_0, gamma_m, kappa] + \
            [omega_m, Omega]

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
        Delta_0 = params[0]
        E_0     = params[1]
        E_1     = params[2]
        G_0     = params[3]
        gamma_m = params[4]
        kappa   = params[5]
        omega_m = params[6]
        Omega   = params[7]
        alpha   = modes[0]
        beta    = modes[1]
        tau     = 2 * np.pi / Omega

        # effective values
        Delta = Delta_0 - np.sqrt(2) * G_0 * np.real(beta)
        G = np.sqrt(2) * G_0 * alpha

        # calculate rates
        dalpha_dt = (- (kappa + 1j * Delta) * alpha + E_0 + E_1 * (np.exp(- 1j * Omega * t * tau) + np.exp(1j * Omega * t * tau)))
        dbeta_dt = (1j * G * np.conjugate(alpha) / 2 - (gamma_m + 1j * omega_m) * beta)

        # arrange rates
        mode_rates = [dalpha_dt * tau, dbeta_dt * tau]

        return mode_rates