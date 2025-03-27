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

__authors__ = ["Sampreet Kalita"]
__toolbox__ = "qom-v1.1.0"
__created__ = "2021-10-22"
__updated__ = "2025-03-11"
__all__     = ['PhysRevA_101_053836']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevA_101_053836(BaseSystem):
    r"""Class to simulate the single-tone modulated QOM system in Phys. Rev. A **101**, 053836 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ========================================================
        key             meaning
        ============    ========================================================
        Delta_a_norm    (*float*) normalized effective detuning of the cavity from the laser :math:`\Delta_{a} / \omega_{m}`. Default is :math:`1.0`.
        G_norms         (*list*) normalized base and sideband amplitudes of the effective coupling strength :math:`[ G_{0} / \omega_{m}, G_{-1} / \omega_{m}, G_{+1} / \omega_{m} ]`. Default is :math:`[ 0.1, 0.01, 0.05 ]`.
        gamma_m_norm    (*float*) normalized mechanical damping rate :math:`\gamma_{m} / \omega_{m}`. Default is :math:`10^{-6}`.
        kappa_norm      (*float*) normalized optical decay rate :math:`\kappa / \omega_{m}`. Default is :math:`0.1`.
        ns              (*list*) quanta of thermal photons and phonons :math:`[ n_{a}, n_{b} ]`. Default is :math:`[ 0.0, 10.0 ]`.
        Omega_norm      (*float*) normalized modulation frequency :math:`\Omega / \omega_{m}`. Default is :math:`2.0`.
        t_rwa           (*bool*) option to work under RWA. Default is `True`.
        ============    ========================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
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
        super().__init__(
            params=params,
            name='PhysRevA_101_053836',
            desc="Single-tone Modulated System in Phys. Rev. A 101, 053836",
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
        G_0_norm, G_m1_norm, G_p1_norm  = self.params['G_norms']

        # with RWA
        if self.params['t_rwa']:
            # substituted expresssions
            G_m_norm    = G_0_norm - G_p1_norm
            G_p_norm    = G_0_norm + G_p1_norm

            # optical position quadrature
            self.A[0][0]    = - self.params['kappa_norm'] / 2.0
            self.A[0][3]    = - G_m_norm
            # optical momentum quadrature
            self.A[1][1]    = - self.params['kappa_norm'] / 2.0
            self.A[1][2]    = G_p_norm
            # mechanical position quadrature
            self.A[2][1]    = - G_m_norm
            self.A[2][2]    = - self.params['gamma_m_norm'] / 2.0
            # mechanical momentum quadrature
            self.A[3][0]    = G_p_norm
            self.A[3][3]    = - self.params['gamma_m_norm'] / 2.0

        # without RWA
        else:
            # effective coupling strength
            G_norm = G_0_norm + G_m1_norm * np.exp(1.0j * self.params['Omega_norm'] * t) + G_p1_norm * np.exp(- 1.0j * self.params['Omega_norm'] * t)

            # optical position quadrature
            self.A[0][0]    = - self.params['kappa_norm'] / 2.0
            self.A[0][1]    = self.params['Delta_a_norm']
            self.A[0][2]    = - 2.0 * np.imag(G_norm) 
            # optical momentum quadrature
            self.A[1][0]    = - self.params['Delta_a_norm']
            self.A[1][1]    = - self.params['kappa_norm'] / 2.0
            self.A[1][2]    = 2.0 * np.real(G_norm)
            # mechanical position quadrature
            self.A[2][2]    = - self.params['gamma_m_norm'] / 2.0
            self.A[2][3]    = 1.0
            # mechanical momentum quadrature
            self.A[3][0]    = 2.0 * np.real(G_norm)
            self.A[3][1]    = 2.0 * np.imag(G_norm)
            self.A[3][2]    = - 1.0
            self.A[3][3]    = - self.params['gamma_m_norm'] / 2.0

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
        gamma_m_norm    = self.params['gamma_m_norm']
        kappa_norm      = self.params['kappa_norm']
        n_a, n_m        = self.params['ns']

        # update noise matrix
        self.D[0][0]    = kappa_norm * (n_a + 0.5)
        self.D[1][1]    = kappa_norm * (n_a + 0.5)
        self.D[2][2]    = gamma_m_norm * (n_m + 0.5)
        self.D[3][3]    = gamma_m_norm * (n_m + 0.5)

        return self.D

    def get_ivc(self):
        """Method to obtain the initial values of the modes, correlations and derived constants and controls.
        
        Returns
        -------
        iv_modes : *numpy.ndarray*
            Initial values of the classical modes.
        iv_corrs : *numpy.ndarray*
            Initial values of the quantum correlations.
        c : *numpy.ndarray*
            Derived constants and controls.
        """

        # extract frequently used variables
        n_a, n_m    = self.params['ns']

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float64)
        iv_corrs[0][0]  = n_a + 0.5 
        iv_corrs[1][1]  = n_a + 0.5
        iv_corrs[2][2]  = n_m + 0.5
        iv_corrs[3][3]  = n_m + 0.5

        return None, iv_corrs, None