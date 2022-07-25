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
__created__ = '2021-05-16'
__updated__ = '2022-07-25'
__version__ = '0.8.5'

# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import SOSMSystem

class PhysRevLett_103_213603(SOSMSystem):
    r"""Class to simulate the modulated QOM system in Phys. Rev. Lett. **103**, 213603 (2009).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ========================================================
        key         meaning
        ========    ========================================================
        F           (float) finesse of the cavity :math:`\mathcal{F}`. Default is :math:`1.4 \times 10^{4}`.
        L           (float) length of the cavity :math:`L` in meters. Default is :math:`25 \times 10^{-3}`.
        lambda_l    (float) wavelength of the laser light :math:`\lambda_{l}` in meters. Default is :math:`1064 \times 10^{-9}`.
        m           (float) mass of the mechanical mirror :math:`m` in kilograms. Default is :math:`150 \times 10^{-12}`.
        omega_m     (float) mechanical freqency :math:`\omega_{m}`. Default is :math:`2 \pi \times 10^{6}`.
        Ps          (float) base and sideband input powers :math:`[ P_{0}, P_{1} ]` in Watts. Default is :math:`[ 10^{-2}, 2 \times 10^{-3} ]`.
        Q           (float) quality factor of the mechanical mirror :math:`Q`. Default is :math:`10^{6}`.
        T           (float) temperature of the mechanical bath :math:`T` in Kelvins. Default is :math:`0.1`.
        ========    ========================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default system parameters
    system_defaults = {
        'F'         : 1.4e4,
        'L'         : 25e-3,
        'lambda_l'  : 1064e-9,
        'm'         : 150e-12,
        'omega_m'   : 2 * np.pi * 1e6,
        'Ps'        : [10e-3, 2e-3],
        'Q'         : 1e6,
        'T'         : 0.1
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevLett_103_213603."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'physrevlett_103_213603'
        self.name = 'Modulated System in Phys. Rev. Lett. 103, 213603'

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
        Delta_0, E_0, E_1 = params[0:3]
        G_0, gamma_m, kappa   = params[3:6]
        Omega, omega_m = params[6:8]
        alpha, beta = modes
        tau = 2 * np.pi / Omega

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
        r"""Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.
            First element contains the optical mode amplitude.
            Next element contains the mechanical mode amplitude.
            Next :math:`4 \times 2^{2}` elements contain the correlations.

        c : list
            Constants of the IVP.
            First :math:`4 \times 2^{2}` elements contain the noise matrix.
            The elements contain the system parameters ``params`` in the following order:
            ========    =============================================
            index       parameter
            ========    =============================================
            0           laser detuning :math:`\Delta_{0}`.
            1           base laser amplitude :math:`E_{0}`.
            2           sideband laser amplitude :math:`E_{1}`.
            4           mechanical damping rate :math:`\gamma_{m}`.
            5           optical decay rate :math:`\kappa`.
            6           modulation frequency :math:`\Omega`.
            7           mechanical frequency :math:`\omega_{m}`.
            ========    =============================================
        """
        
        # extract frequently used variables
        F       = self.params['F']
        L       = self.params['L']
        lambda_l= self.params['lambda_l']
        m       = self.params['m']
        omega_m = self.params['omega_m']
        P_0, P_1= self.params['Ps']
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
        params = [Delta_0, E_0, E_1] + \
            [G_0, gamma_m, kappa] + \
            [Omega, omega_m]

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
        Delta_0, E_0, E_1 = params[0:3]
        G_0, gamma_m, kappa   = params[3:6]
        Omega, omega_m = params[6:8]
        alpha, beta = modes
        tau = 2 * np.pi / Omega

        # effective values
        Delta = Delta_0 - np.sqrt(2) * G_0 * np.real(beta)
        G = np.sqrt(2) * G_0 * alpha

        # calculate rates
        dalpha_dt = (- (kappa + 1j * Delta) * alpha + E_0 + E_1 * (np.exp(- 1j * Omega * t * tau) + np.exp(1j * Omega * t * tau)))
        dbeta_dt = (1j * G * np.conjugate(alpha) / 2 - (gamma_m + 1j * omega_m) * beta)

        # arrange rates
        mode_rates = [dalpha_dt * tau, dbeta_dt * tau]

        return mode_rates