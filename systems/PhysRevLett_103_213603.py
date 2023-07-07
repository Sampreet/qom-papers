#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the modulated QOM system in Phys. Rev. Lett. **103**, 213603 (2009).

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
__toolbox__ = 'qom-v1.0.0'
__created__ = '2021-05-16'
__updated__ = '2023-07-07'

# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import BaseSystem

class PhysRevLett_103_213603(BaseSystem):
    r"""Class to simulate the modulated QOM system in Phys. Rev. Lett. **103**, 213603 (2009).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ========================================================
        key         meaning
        ========    ========================================================
        F           (*float*) finesse of the cavity :math:`\mathcal{F}`. Default is :math:`1.4 \times 10^{4}`.
        L           (*float*) length of the cavity :math:`L` in meters. Default is :math:`25 \times 10^{-3}`.
        lambda_l    (*float*) wavelength of the laser light :math:`\lambda_{l}` in meters. Default is :math:`1064 \times 10^{-9}`.
        m           (*float*) mass of the mechanical mirror :math:`m` in kilograms. Default is :math:`150 \times 10^{-12}`.
        omega_m     (*float*) mechanical freqency :math:`\omega_{m}`. Default is :math:`2 \pi \times 10^{6}`.
        Ps          (*float*) base and sideband input powers :math:`[ P_{0}, P_{1} ]` in Watts. Default is :math:`[ 10^{-2}, 2 \times 10^{-3} ]`.
        Q           (*float*) quality factor of the mechanical mirror :math:`Q`. Default is :math:`10^{6}`.
        T           (*float*) temperature of the mechanical bath :math:`T` in Kelvins. Default is :math:`0.1`.
        ========    ========================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'F'         : 1.4e4,
        'L'         : 25e-3,
        'lambda_l'  : 1064e-9,
        'm'         : 150e-12,
        'omega_m'   : 2.0 * np.pi * 1e6,
        'Ps'        : [10e-3, 2e-3],
        'Q'         : 1e6,
        'T'         : 0.1
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_103_213603."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevLett_103_213603',
            desc='Modulated System in Phys. Rev. Lett. 103, 213603',
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
        Delta_0, _, _, G_0      = c[0:4]
        gamma_m, kappa, Omega   = c[4:7]
        alpha, beta             = modes

        # effective values
        Delta   = Delta_0 - np.sqrt(2.0) * G_0 * np.real(beta)
        G       = np.sqrt(2.0) * G_0 * alpha

        # optical position quadrature
        self.A[0][0]    = - kappa 
        self.A[0][1]    = Delta 
        self.A[0][2]    = - np.imag(G) 
        # optical momentum quadrature
        self.A[1][0]    = - Delta
        self.A[1][1]    = - kappa
        self.A[1][2]    = np.real(G)
        # mechanical position quadrature
        self.A[2][3]    = self.params['omega_m']
        # mechanical momentum quadrature
        self.A[3][0]    = np.real(G)
        self.A[3][1]    = np.imag(G)
        self.A[3][2]    = - self.params['omega_m']
        self.A[3][3]    = - gamma_m
        # normalize
        self.A *= 2.0 * np.pi / Omega

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
        gamma_m, kappa, Omega = c[4:7]

        # update noise matrix
        self.D[0][0]    = kappa
        self.D[1][1]    = kappa
        self.D[3][3]    = gamma_m
        # normalize
        self.D *= 2.0 * np.pi / Omega
        
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
            0           laser detuning :math:`\Delta_{0}`.
            1           base laser amplitude :math:`E_{0}`.
            2           sideband laser amplitude :math:`E_{1}`.
            3           coupling strength :math:`G_{0}`.
            4           mechanical damping rate :math:`\gamma_{m}`.
            5           optical decay rate :math:`\kappa`.
            6           modulation frequency :math:`\Omega`.
            ========    =============================================
        """
        
        # extract frequently used variables
        F           = self.params['F']
        L           = self.params['L']
        lambda_l    = self.params['lambda_l']
        m           = self.params['m']
        omega_m     = self.params['omega_m']
        P_0, P_1    = self.params['Ps']
        Q           = self.params['Q']
        T           = self.params['T']

        # laser detuning
        Delta_0 = omega_m
        # mechanical decay rate
        gamma_m = omega_m / Q
        # optical decay rate
        kappa = np.pi * sc.c / (2.0 * F * L)
        # laser frequency
        omega_l = 2.0 * np.pi * sc.c / lambda_l
        # cavity frequency
        omega_c = Delta_0 + omega_l
        # coupling strength
        G_0 = np.sqrt(sc.hbar / (m * omega_m)) * omega_c / L
        # modulation frequency
        Omega = 2.0 * omega_m

        # thermal phonon number
        n_th = 0.0 if T == 0.0 else 1.0 / (np.exp(sc.hbar * omega_m / (sc.k * T)) - 1.0)

        # laser amplitudes
        E_0 = np.sqrt(2 * kappa * P_0 / (sc.hbar * omega_l)) 
        E_1 = np.sqrt(2 * kappa * P_1 / (sc.hbar * omega_l))
 
        # initial values of the modes
        iv_modes = np.zeros(self.num_modes, dtype=np.complex_)

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = 0.5 
        iv_corrs[1][1]  = 0.5
        iv_corrs[2][2]  = n_th + 0.5
        iv_corrs[3][3]  = n_th + 0.5
        
        # derived constants
        c = np.array([Delta_0, E_0, E_1, G_0, gamma_m, kappa, Omega], dtype=np.float_)

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
        Delta_0, E_0, E_1, G_0  = c[0:4]
        gamma_m, kappa, Omega   = c[4:7]
        alpha, beta             = modes
        tau                     = 2 * np.pi / Omega

        # effective values
        Delta = Delta_0 - np.sqrt(2) * G_0 * np.real(beta)
        G = np.sqrt(2) * G_0 * alpha

        # calculate rates
        dalpha_dt = (- (kappa + 1j * Delta) * alpha + E_0 + E_1 * (np.exp(- 1j * Omega * t * tau) + np.exp(1j * Omega * t * tau)))
        dbeta_dt = (1j * G * np.conjugate(alpha) / 2 - (gamma_m + 1j * self.params['omega_m']) * beta)

        # arrange rates, normalize and return
        return np.array([dalpha_dt, dbeta_dt], dtype=np.complex_) * tau