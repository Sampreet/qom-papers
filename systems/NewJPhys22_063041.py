#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the systems described in New J. Phys. **22**, 063041 (2020).

# Stationary Quantum Entanglement between a Massive Mechanical Membrane and a Low Frequency LC Circuit

==================  =====================================================
Author              Affiliation
==================  =====================================================
Jie Li              Delft University of Technology, 2628CJ Delft, The Netherlands
Simon Groblacher    Delft University of Technology, 2628CJ Delft, The Netherlands
==================  =====================================================

## Objectives
* Entangling a massive mechanical oscillator to a macroscopic low-frequency LC circuit.

## Novelty
* LC resonator frequency in the radio-frequency domain close to the mechanical frequency.
* Membrane-LC interaction is quadratic in LC charge and linear in mechanical position.
* DC drive to enhance the electromechanical coupling rate.
* Red-detuned laser to cool the mechanical mode.

## Assumptions
* Nearly resonant mechanical and LC oscillators.
* Linearized description of the modes.

## Results
* Strong optomechanical and electromechanical rates are required.
* Entanglement is robust against temperature.
"""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-05-15'
__updated__ = '2021-08-30'
__version__ = '0.8.0'


# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import SODMSystem

class NewJPhys22_063041(SODMSystem):
    """Class to simulate the OEM system in New J. Phys. **22**, 063041 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for NewJPhys22_063041."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'new_j_phys_22_063041'
        self.name = 'OEM System in New J. Phys. 22, 063041'  
        # default parameters
        self.params = {
            'Delta_norm': params.get('Delta_norm', 1.0),
            'G_norm': params.get('G_norm', 3.0),
            'g_norm': params.get('g_norm', 3.0),
            'gamma_LC_norm': params.get('gamma_LC_norm', 1e-5),
            'gamma_m_norm': params.get('gamma_m_norm', 1e-6),
            'kappa_norm': params.get('kappa_norm', 0.1),
            'omega_LC_norm': params.get('omega_LC_norm', 1.0),
            'omega_m': params.get('omega_m', 2e6 * np.pi),
            'T_LC': params.get('T_LC', 1e-2),
            'T_m': params.get('T_m', 1e-2),
            't_Delta': params.get('t_Delta', 'absolute')
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
        Delta           = params[0]
        G               = params[1]
        g               = params[2]
        gamma_LC        = params[3]
        gamma_m         = params[4]
        kappa           = params[5]
        omega_LC_prime  = params[6]
        omega_m         = params[7]

        # drift matrix
        if self.A is None or np.shape(self.A) != (6, 6):
            self.A = np.zeros([6, 6], dtype=np.float_)
        # optical quadratures
        self.A[0][0] = - kappa
        self.A[0][1] = Delta
        self.A[1][0] = - Delta
        self.A[1][1] = - kappa
        self.A[1][2] = G
        # first mechanical mode quadratures
        self.A[2][3] = omega_m
        self.A[3][0] = G
        self.A[3][2] = - omega_m
        self.A[3][3] = - gamma_m
        self.A[3][4] = - g
        # second mechanical mode quadratures
        self.A[4][5] = omega_LC_prime
        self.A[5][2] = - g
        self.A[5][4] = - omega_LC_prime
        self.A[5][5] = - gamma_LC

        return self.A

    def get_ivc(self):
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.

        c : list
            Constant parameters of the system.
            First element contains the laser detuning.
            Next two elements contain the coupling constants.
            Next two elements contain the LC circuit and mechanical damping rates.
            Next element contains the optical decay rate.
            Next two elements contain the LC circuit and mechanical frequencies.
        """

        # extract frequently used variables
        Delta_norm      = self.params['Delta_norm']
        G_norm          = self.params['G_norm']
        g_norm          = self.params['g_norm']
        gamma_LC_norm   = self.params['gamma_LC_norm']
        gamma_m_norm    = self.params['gamma_m_norm']
        kappa_norm      = self.params['kappa_norm']
        omega_LC_norm   = self.params['omega_LC_norm']
        omega_m         = self.params['omega_m']
        T_LC            = self.params['T_LC']
        T_m             = self.params['T_m']
        t_Delta         = self.params['t_Delta']

        # effective cavity-laser detuning
        if t_Delta == 'relative':
            Delta = Delta_norm * omega_LC_norm * omega_m
        else:
            Delta = Delta_norm * omega_m
        # LC frequency
        omega_LC = omega_LC_norm * omega_m
        # effective LC frequency
        omega_LC_prime = omega_LC

        # update decays
        G       = G_norm * kappa_norm * omega_m
        g       = g_norm * kappa_norm * omega_m
        gamma_LC= gamma_LC_norm * omega_LC
        gamma_m = gamma_m_norm * omega_m
        kappa   = kappa_norm * omega_m
        n_LC    = sc.k * T_LC / sc.hbar / omega_LC
        n_m     = sc.k * T_m / sc.hbar / omega_m

        # effective cavity-laser detuning
        if t_Delta == 'relative':
            Delta = Delta_norm * omega_LC_norm * omega_m
        else:
            Delta = Delta_norm * omega_m
        # LC frequency
        omega_LC = omega_LC_norm * omega_m
        # effective LC frequency
        omega_LC_prime = omega_LC
 
        # initial mode values as 1D list
        modes_0 = np.zeros(3, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros([6, 6], dtype=np.float_)
        corrs_0[0][0] = 1/2 
        corrs_0[1][1] = 1/2
        corrs_0[2][2] = (n_m + 1/2)
        corrs_0[3][3] = (n_LC + 1/2)

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros([6, 6], dtype=np.float_)
        D[0][0] = kappa
        D[1][1] = kappa
        D[3][3] = gamma_m * (2 * n_m + 1)
        D[5][5] = gamma_LC * (2 * n_LC + 1)
        
        # constant parameters
        params = [Delta] + \
            [G, g] + \
            [gamma_LC, gamma_m, kappa] + \
            [omega_LC_prime, omega_m]

        # all constants
        c = D.flatten().tolist() + params

        return iv, c