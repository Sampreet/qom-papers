#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the systems described in Phys. Rev. A **101**, 053836 (2020).

# Strong Mechanical Squeezing in a Standard Optomechanical System by Pump Modulation

==============  =====================================================
Author          Affiliation
==============  =====================================================
C.-H. Bai       School of Physics, Harbin Institute of Technology, Harbin, Heilongjiang 150001, China
D.-Y. Wang      School of Physics, Harbin Institute of Technology, Harbin, Heilongjiang 150001, China
S. Zhang        School of Physics, Harbin Institute of Technology, Harbin, Heilongjiang 150001, China
                Department of Physics, College of Science, Yanbian University, Yanji, Jilin 133002, China
S. Liu          School of Physics, Harbin Institute of Technology, Harbin, Heilongjiang 150001, China
H.-F. Wang      Department of Physics, College of Science, Yanbian University, Yanji, Jilin 133002, China
==============  =====================================================

## Objectives
* Manipulating a micromechanical oscillator into a quantum squeezed stateusing single-tone pump.

## Novelty
* A single-tone driving field.
* Optimizing the ratio of the effective optomechanical coupling sideband strengths.

## Assumptions

## Results
* Cooling down of Bogoliobov mode of the mechanical mode.
* Strong mechanical sqeezing far surpassing the SQL.
* Robust against mechanical thermal noises and survives at high bath temperatures.
"""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-10-22'
__updated__ = '2022-01-09'
__version__ = '0.8.2'

# dependencies
import numpy as np

# qom modules
from qom.systems import SOSMSystem

class PhysRevA101_053836(SOSMSystem):
    """Class to simulate the single-tone modulated QOM system in Phys. Rev. A **101**, 053836 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA101_053836."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'phys_rev_a_101_053836'
        self.name = 'Single-tone Modulated QOM System in Phys. Rev. A 101, 053836'
        # default parameters
        self.params = {
            'Delta_a': params.get('Delta_a', 1.0),
            'G_0': params.get('G_0', 0.1),
            'G_p1': params.get('G_p1', 0.05),
            'G_m1': params.get('G_m1', 0.01),
            'gamma_m': params.get('gamma_m', 1e-6),
            'kappa': params.get('kappa', 0.1),
            'omega_m': params.get('omega_m', 1.0),
            'Omega': params.get('Omega', 2.0),
            'n_a': params.get('n_a', 0),
            'n_m': params.get('n_m', 10)
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
        Delta_a = params[0]
        G_0     = params[1]
        G_p1    = params[2]
        G_m1    = params[3]
        gamma_m = params[4]
        kappa   = params[5]
        omega_m = params[6]
        Omega   = params[7]

        # effective coupling strength
        G = G_m1 * np.exp(1j * Omega * t) + G_0 + G_p1 * np.exp(-1j * Omega * t)

        # initialize drift matrix
        if self.A is None or np.shape(self.A) != (4, 4):
            self.A = np.zeros([4, 4], dtype=np.float_)
        # optical position quadrature
        self.A[0][0] = - kappa / 2
        self.A[0][1] = Delta_a
        self.A[0][2] = - 2 * np.imag(G) 
        # optical momentum quadrature
        self.A[1][0] = - Delta_a
        self.A[1][1] = - kappa / 2
        self.A[1][2] = 2 * np.real(G)
        # mechanical position quadrature
        self.A[2][2] = - gamma_m / 2
        self.A[2][3] = omega_m
        # mechanical momentum quadrature
        self.A[3][0] = 2 * np.real(G)
        self.A[3][1] = 2 * np.imag(G)
        self.A[3][2] = - omega_m
        self.A[3][3] = - gamma_m / 2

        return self.A

    def get_ivc(self):
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.

        c : list
            Constant parameters of the system.
            First (2*2)^2 elements contain the noise matrix.
            First element contains the effective laser detuning.
            Next two elements contain the effective base and modulated coupling strengths.
            Next two elements contain the mechanical and optical decay rates.
            Next two elements contain the mechanical and modulation frequencies.
        """

        # extract frequently used variables
        Delta_a = self.params['Delta_a']
        G_0     = self.params['G_0']
        G_p1    = self.params['G_p1']
        G_m1    = self.params['G_m1']
        gamma_m = self.params['gamma_m']
        kappa   = self.params['kappa']
        omega_m = self.params['omega_m']
        Omega   = self.params['Omega']
        n_a     = self.params['n_a']
        n_m     = self.params['n_m']
 
        # initial mode values as 1D list
        modes_0 = np.zeros(2, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros([4, 4], dtype=np.float_)
        corrs_0[0][0] = n_a + 1 / 2 
        corrs_0[1][1] = n_a + 1 / 2
        corrs_0[2][2] = n_m + 1 / 2
        corrs_0[3][3] = n_m + 1 / 2

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros([4, 4], dtype=np.float_)
        D[0][0] = kappa * (2 * n_a + 1) / 2
        D[1][1] = kappa * (2 * n_a + 1) / 2
        D[2][2] = gamma_m * (2 * n_m + 1) / 2
        D[3][3] = gamma_m * (2 * n_m + 1) / 2
        
        # constant parameters
        params = [Delta_a] + \
            [G_0, G_p1, G_m1] + \
            [gamma_m, kappa] + \
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
        Delta_a     = params[0]
        G_0         = params[1]
        G_p1        = params[2]
        G_m1        = params[3]
        gamma_m     = params[4]
        kappa       = params[5]
        omega_m     = params[6]
        Omega       = params[7]
        alpha       = modes[0]
        beta        = modes[1]

        # effective coupling strength
        G = G_m1 * np.exp(1j * Omega * t) + G_0 + G_p1 * np.exp(-1j * Omega * t)

        # arrange rates
        mode_rates = [1, 1]

        return mode_rates

class PhysRevA101_053836_RWA(SOSMSystem):
    """Class to simulate the single-tone modulated QOM system in Phys. Rev. A **101**, 053836 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA101_053836."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'phys_rev_a_101_053836_rwa'
        self.name = 'Single-tone Modulated QOM System in Phys. Rev. A 101, 053836 with RWA'
        # default parameters
        self.params = {
            'Delta_a': params.get('Delta_a', 1.0),
            'G_0': params.get('G_0', 0.1),
            'G_p1': params.get('G_p1', 0.05),
            'G_m1': params.get('G_m1', 0.01),
            'gamma_m': params.get('gamma_m', 1e-6),
            'kappa': params.get('kappa', 0.1),
            'omega_m': params.get('omega_m', 1.0),
            'Omega': params.get('Omega', 2.0),
            'n_a': params.get('n_a', 0),
            'n_m': params.get('n_m', 10)
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
        G_0     = params[1]
        G_p1    = params[2]
        gamma_m = params[4]
        kappa   = params[5]

        # effective coupling strength
        G_p = G_0 + G_p1
        G_m = G_0 - G_p1

        # initialize drift matrix
        if self.A is None or np.shape(self.A) != (4, 4):
            self.A = np.zeros([4, 4], dtype=np.float_)
        # optical position quadrature
        self.A[0][0] = - kappa / 2
        self.A[0][3] = - G_m
        # optical momentum quadrature
        self.A[1][1] = - kappa / 2
        self.A[1][2] = G_p
        # mechanical position quadrature
        self.A[2][1] = - G_m
        self.A[2][2] = - gamma_m / 2
        # mechanical momentum quadrature
        self.A[3][0] = G_p
        self.A[3][3] = - gamma_m / 2

        return self.A

    def get_ivc(self):
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.

        c : list
            Constant parameters of the system.
            First (2*2)^2 elements contain the noise matrix.
            First element contains the effective laser detuning.
            Next two elements contain the effective base and modulated coupling strengths.
            Next two elements contain the mechanical and optical decay rates.
            Next two elements contain the mechanical and modulation frequencies.
        """

        # extract frequently used variables
        Delta_a = self.params['Delta_a']
        G_0     = self.params['G_0']
        G_p1    = self.params['G_p1']
        G_m1    = self.params['G_m1']
        gamma_m = self.params['gamma_m']
        kappa   = self.params['kappa']
        omega_m = self.params['omega_m']
        Omega   = self.params['Omega']
        n_a     = self.params['n_a']
        n_m     = self.params['n_m']
 
        # initial mode values as 1D list
        modes_0 = np.zeros(2, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros([4, 4], dtype=np.float_)
        corrs_0[0][0] = n_a + 1 / 2 
        corrs_0[1][1] = n_a + 1 / 2
        corrs_0[2][2] = n_m + 1 / 2
        corrs_0[3][3] = n_m + 1 / 2

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros([4, 4], dtype=np.float_)
        D[0][0] = kappa * (2 * n_a + 1) / 2
        D[1][1] = kappa * (2 * n_a + 1) / 2
        D[2][2] = gamma_m * (2 * n_m + 1) / 2
        D[3][3] = gamma_m * (2 * n_m + 1) / 2
        
        # constant parameters
        params = [Delta_a] + \
            [G_0, G_p1, G_m1] + \
            [gamma_m, kappa] + \
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

        # arrange rates
        mode_rates = [1, 1]

        return mode_rates