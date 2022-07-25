#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the hybrid OEM system in New J. Phys. **22**, 063041 (2020).

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
__updated__ = '2022-07-24'
__version__ = '0.8.5'


# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import SODMSystem

class NewJPhys_22_063041(SODMSystem):
    r"""Class to simulate the hybrid OEM system in New J. Phys. **22**, 063041 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        Delta_norm      (float) normalized laser detuning :math:`\Delta / \omega_{m}` when ``t_Delta`` is "absolute" and :math:`\Delta / \omega_{LC}` when ``t_Delta`` is "relative". Default is :math:`1.0`.
        G_norm          (float) normalized effective optomechanical coupling strength :math:`G / \kappa`. Default is :math:`3.0`.
        g_norm          (float) normalized effective electromechanical coupling strength :math:`g / \kappa`. Default is :math:`3.0`.
        gamma_LC_norm   (float) normalized LC circuit damping rate :math:`\gamma_{LC} / \omega_{LC}`. Default is :math:`10^{-5}`.
        gamma_m_norm    (float) normalized mechanical damping rate :math:`\gamma_{m} / \omega_{m}`. Default is :math:`10^{-6}`.
        kappa_norm      (float) normalized optical decay rate :math:`\kappa / \omega_{m}`. Default is :math:`0.1`.
        omega_LC_norm   (float) normalized LC circuit frequency :math:`\omega_{LC} / \omega_{m}`. Default is :math:`1.0`.
        omega_m         (float) mechanical frequency :math:`\omega_{m}`. 
        T_LC            (float) LC circuit bath temperature :math:`T_{LC}`. Default is :math:`10^{-2}`.
        T_m             (float) mechanical bath temperature :math:`T_{m}`. Default is :math:`10^{-2}`.
        t_Delta         (float) type of :math:`\Delta`. Options are "absolute" when normalized by :math:`\omega_{m}` and "relative" when normalized by :math:`\omega_{LC}`. Default is "absolute". 
        ============    ========================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default system parameters
    system_defaults = {
        'Delta_norm'    : 1.0,
        'G_norm'        : 3.0,
        'g_norm'        : 3.0,
        'gamma_LC_norm' : 1e-5,
        'gamma_m_norm'  : 1e-6,
        'kappa_norm'    : 0.1,
        'omega_LC_norm' : 1.0,
        'omega_m'       : 2 * np.pi * 1e6,
        'T_LC'          : 1e-2,
        'T_m'           : 1e-2,
        't_Delta'       : 'absolute'
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for NewJPhys_22_063041."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'newjphys_22_063041'
        self.name = 'Hybrid System in New J. Phys. 22, 063041'  

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
        Delta, G, g = params[0:3]
        gamma_LC,  gamma_m, kappa = params[3:6]
        omega_LC_prime, omega_m = params[6:8]

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
        r"""Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.
            First element contains the optical mode amplitude.
            Next element contains the mechanical mode amplitude.
            Next :math:`4 \times 3^{2}` elements contain the correlations.

        c : list
            Constants of the IVP.
            First :math:`4 \times 3^{2}` elements contain the noise matrix.
            Rest of the elements contain the system parameters ``params`` in the following order:
            ========    =============================================
            index       parameter
            ========    =============================================
            0           laser detuning :math:`\Delta`.
            1           effective optomechanical coupling strength :math:`G`.
            2           effective electromechanical coupling strength :math:`g`.
            3           LC circuit damping rate :math:`\gamma_{LC}`.
            4           mechanical damping rate :math:`\gamma_{m}`.
            5           optical decay rate :math:`\kappa`.
            6           effective LC circuit frequency :math:`\omega_{LC}^{\prime}`.
            7           mechanical frequency :math:`\omega_{m}`.
            ========    =============================================
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
        params = [Delta, G, g] + \
            [gamma_LC, gamma_m, kappa] + \
            [omega_LC_prime, omega_m]

        # all constants
        c = D.flatten().tolist() + params

        return iv, c