#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the moveable end-mirror QOM system in Phys. Rev. Lett. **98**, 030405 (2007)."""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.0'
__created__ = '2022-07-31'
__updated__ = '2023-07-07'

# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import BaseSystem

class PhysRevLett_98_030405(BaseSystem):
    r"""Class to simulate the moveable end-mirror QOM system in Phys. Rev. Lett. **98**, 030405 (2007).
    
    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        Delta_norm  (*float*) normalized effective detuning of the cavity from the laser :math:`\Delta / \omega_{m}`. Default is :math:`1.0`.
        F           (*float*) Finessed of the cavity :math:`\mathcal{F}`. Default is :math:`1.07 \times 10^{4}`.
        gamma_m     (*float*) mechanical damping rate :math:`gamma_{m}` in Hertz. Default is :math:`2 \pi \times 100` Hz.
        L           (*float*) length of the cavity :math:`L` in metres. Default is :math:`10^{-3}` m.
        lamb        (*float*) wavelength of the laser :math:`lambda` in metres. Default is :math:`810 \times 10^{-9}` m.
        m           (*float*) mass of the mirror :math:`m` in kilograms. Default is :math:`5 \times 10^{12}` kg.
        omega_m     (*float*) mechanical frequency :math:`\omega_{m}` in Hertz.  Default is :math:`2 \pi \times 10^{7}` Hz.
        P           (*float*) power of the laser :math:`P` in Watts. Default is :math:`50 \times 10^{-3}` W.
        T           (*float*) temperature of the bath :math:`T` in Kelvins. Default is :math:`400 \times 10^{-3}` K.
        ========    ============================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'Delta_norm': 1.0,
        'F'         : 1.07e4,
        'gamma_m'   : 2.0 * np.pi * 1e2,
        'L'         : 1e-3,
        'lamb'      : 810e-9,
        'm'         : 5e-12,
        'omega_m'   : 2.0 * np.pi * 10e6,
        'P'         : 50e-3,
        'T'         : 400e-3
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_98_030405."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevLett_98_030405',
            desc='Moveable End-mirror System in Phys. Rev. Lett. 98, 030405',
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
        Delta, _, G_0, kappa    = c
        alpha, _                = modes

        # effective values
        G = np.sqrt(2) * G_0 * alpha
        
        # X quadratures
        self.A[0][0]    = - kappa
        self.A[0][1]    = Delta
        # Y quadratures
        self.A[1][0]    = - Delta
        self.A[1][1]    = - kappa
        self.A[1][2]    = G
        # Q quadratures
        self.A[2][3]    = self.params['omega_m']
        # P quadratures
        self.A[3][0]    = G
        self.A[3][2]    = - self.params['omega_m']
        self.A[3][3]    = - self.params['gamma_m']

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
        T       = self.params['T']
        kappa   = c[3]
        n_bar   = 0.0 if T == 0.0 else 1.0 / (np.exp(sc.hbar * self.params['omega_m'] / sc.k / T) - 1)

        # update drift matrix
        self.D[0][0]    = kappa
        self.D[1][1]    = kappa
        self.D[3][3]    = self.params['gamma_m'] * (2.0 * n_bar + 1.0)

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
            ========    =====================================================
            index       parameter
            ========    =====================================================
            0           effective detuning :math:`\Delta`.
            1           laser amplitude :math:`E`.
            2           coupling constant :math:`G_{0}`.
            3           optical decay rate :math:`\kappa`.
            ========    =====================================================
        """

        # extract frequently used variables
        Delta_norm  = self.params['Delta_norm']
        F           = self.params['F']
        gamma_m     = self.params['gamma_m']
        L           = self.params['L']
        lamb        = self.params['lamb']
        m           = self.params['m']
        omega_m     = self.params['omega_m']
        P           = self.params['P']
        T           = self.params['T']
        Delta       = Delta_norm * omega_m
        kappa       = np.pi * sc.c / L / F
        omega_l     = 2.0 * np.pi * sc.c / lamb
        E           = np.sqrt(2.0 * P * kappa / sc.hbar / omega_l)
        omega_c     = omega_l + Delta
        G_0         = omega_c / L * np.sqrt(sc.hbar / m / omega_m)
        n_bar       = 0.0 if T == 0.0 else 1.0 / (np.exp(sc.hbar * omega_m / sc.k / T) - 1.0)

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = 0.5 
        iv_corrs[1][1]  = 0.5
        iv_corrs[2][2]  = n_bar + 0.5
        iv_corrs[3][3]  = n_bar + 0.5
        
        # constant parameters
        c = np.array([Delta, E, G_0, kappa], dtype=np.float_)

        return None, iv_corrs, c

    def get_modes_steady_state(self, c):
        """Method to obtain the steady state modes.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns 
        -------
        Modes : *numpy.ndarray*
            Steady state modes with shape `(dim, num_modes)`.
        """

        # extract frequently used variables
        Delta, E, _, kappa = c

        return np.array([[
            np.abs(E / (kappa + 1.0j * Delta)),
            1.0
        ]], dtype=np.complex_)
