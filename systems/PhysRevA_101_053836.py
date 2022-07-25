#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Module to simulate the single-tone modulated QOM system in Phys. Rev. A **101**, 053836 (2020).

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
__updated__ = '2022-07-25'
__version__ = '0.8.5'

# dependencies
import numpy as np

# qom modules
from qom.systems import SOSMSystem

class PhysRevA_101_053836(SOSMSystem):
    r"""Class to simulate the single-tone modulated QOM system in Phys. Rev. A **101**, 053836 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        Delta_a_norm    (float) normalized effective detuning of the cavity from the laser :math:`\Delta_{a} / \omega_{m}`. Default is :math:`1.0`.
        G_norms         (list) normalized base and sideband amplitudes of the effective coupling strength :math:`[ G_{0} / \omega_{m}, G_{-1} / \omega_{m}, G_{+1} / \omega_{m} ]`. Default is :math:`[ 0.1, 0.01, 0.05 ]`.
        gamma_m_norm    (float) normalized mechanical damping rate :math:`\gamma_{m} / \omega_{m}`. Default is :math:`10^{-6}`.
        kappa_norm      (float) normalized optical decay rate :math:`\kappa / \omega_{m}`. Default is :math:`0.1`.
        ns              (list) quanta of thermal photons and phonons :math:`[ n_{a}, n_{b} ]`. Default is :math:`[ 0.0, 10.0 ]`.
        Omega_norm      (float) normalized modulation frequency :math:`\Omega / \omega_{m}`. Default is :math:`2.0`.
        t_rwa           (bool) option to work under RWA. Default is `True`.
        ============    ========================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    # default system parameters
    system_defaults = {
        'Delta_a_norm'  : 1.0,
        'G_norms'       : [0.1, 0.01, 0.05],
        'gamma_m_norm'  : 1e-6,
        'kappa_norm'    : 0.1,
        'ns'            : [0.0, 10.0],
        'Omega_norm'    : 2.0,
        't_rwa'         : True
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA_101_053836."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'physreva_101_053836'
        self.name = 'Single-tone Modulated System in Phys. Rev. A 101, 053836'

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
        Delta_a_norm, G_0_norm, G_m1_norm, G_p1_norm = params[0:4]
        gamma_m_norm, kappa_norm, Omega_norm = params[4:7]
        t_rwa = params[7]

        # initialize drift matrix
        if self.A is None or np.shape(self.A) != (4, 4):
            self.A = np.zeros([4, 4], dtype=np.float_)

        # with RWA
        if t_rwa:
            # substituted expresssions
            G_m_norm = G_0_norm - G_p1_norm
            G_p_norm = G_0_norm + G_p1_norm

            # optical position quadrature
            self.A[0][0] = - kappa_norm / 2
            self.A[0][3] = - G_m_norm
            # optical momentum quadrature
            self.A[1][1] = - kappa_norm / 2
            self.A[1][2] = G_p_norm
            # mechanical position quadrature
            self.A[2][1] = - G_m_norm
            self.A[2][2] = - gamma_m_norm / 2
            # mechanical momentum quadrature
            self.A[3][0] = G_p_norm
            self.A[3][3] = - gamma_m_norm / 2

        # without RWA
        else:
            # effective coupling strength
            G_norm = G_0_norm + G_m1_norm * np.exp(1j * Omega_norm * t) + G_p1_norm * np.exp(-1j * Omega_norm * t)

            # optical position quadrature
            self.A[0][0] = - kappa_norm / 2
            self.A[0][1] = Delta_a_norm
            self.A[0][2] = - 2 * np.imag(G_norm) 
            # optical momentum quadrature
            self.A[1][0] = - Delta_a_norm
            self.A[1][1] = - kappa_norm / 2
            self.A[1][2] = 2 * np.real(G_norm)
            # mechanical position quadrature
            self.A[2][2] = - gamma_m_norm / 2
            self.A[2][3] = 1.0
            # mechanical momentum quadrature
            self.A[3][0] = 2 * np.real(G_norm)
            self.A[3][1] = 2 * np.imag(G_norm)
            self.A[3][2] = - 1.0
            self.A[3][3] = - gamma_m_norm / 2

        return self.A

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
            0           normalized effective detuning :math:`\Delta_{a} / \omega_{m}`.
            1           normalized base amplitude of the effective coupling strength :math:`G_{0} / \omega_{m}`.
            2           normalized minus sideband amplitude of the effective coupling strength :math:`G_{-1} / \omega_{m}`.
            3           normalized minus sideband amplitude of the effective coupling strength :math:`G_{+1} / \omega_{m}`.
            4           normalized mechanical damping rate :math:`\gamma_{m} / \omega_{m}`.
            5           normalized optical decay rate :math:`\kappa / \omega_{m}`.
            6           normalized modulation frequency :math:`\Omega / \omega_{m}`.
            7           option to work under RWA ``t_rwa``.
            ========    =============================================
        """

        # extract frequently used variables
        Delta_a_norm= self.params['Delta_a_norm']
        G_norms     = self.params['G_norms']
        gamma_m_norm= self.params['gamma_m_norm']
        kappa_norm  = self.params['kappa_norm']
        n_a, n_m    = self.params['ns']
        Omega_norm  = self.params['Omega_norm']
        t_rwa       = self.params['t_rwa']
 
        # initial mode values as 1D list
        modes_0 = np.zeros(2, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros([4, 4], dtype=np.float_)
        corrs_0[0][0] = n_a + 0.5 
        corrs_0[1][1] = n_a + 0.5
        corrs_0[2][2] = n_m + 0.5
        corrs_0[3][3] = n_m + 0.5

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros([4, 4], dtype=np.float_)
        D[0][0] = kappa_norm * (n_a + 0.5)
        D[1][1] = kappa_norm * (n_a + 0.5)
        D[2][2] = gamma_m_norm * (n_m + 0.5)
        D[3][3] = gamma_m_norm * (n_m + 0.5)
        
        # constant parameters
        params = [Delta_a_norm] + G_norms + \
            [gamma_m_norm, kappa_norm, Omega_norm] + \
            [t_rwa]

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

        return [1, 1]