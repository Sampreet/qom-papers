#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the PT-symmetric QOM systems in Phys. Rev. A **100**, 063846 (2019)."""

__authors__ = ['Sampreet Kalita']
__toolbox__ = 'qom-v1.0.2'
__created__ = '2023-09-13'
__updated__ = '2024-06-23'
__all__     = ['PhysRevA_100_063846_00', 'PhysRevA_100_063846_01']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevA_100_063846_00(BaseSystem):
    r"""Class to simulate the Bipartite PT-symmetric QOM system in Phys. Rev. A **100**, 063846 (2019).
    
    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        gamma_norm  (*float*) normalized mechanical detuning rate :math:`\gamma / \Gamma`. Default is :math:`10^{-3}`.
        J_norm      (*float*) normalized coupling strength :math:`J / \Gamma`. Default is :math:`0.5`.
        n_th        (*float*) average thermal phonon occupancy :math:`n_{th}`. Default is :math:`0.0`.
        is_noisy    (*bool*) option to add or remove noises. Default is `True`.
        ========    ============================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA_100_063846_00."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevA_100_063846_00',
            desc='Bipartite PT-symmetric QOM system in Phys. Rev. A 100, 063846',
            num_modes=2,
            cb_update=cb_update
        )

        # attributes
        self.is_A_constant = True

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
        
        # update drift matrix
        self.A[0][0]    = 0.5 - self.params['gamma_norm'] / 2.0
        self.A[0][3]    = - self.params['J_norm']
        self.A[1][1]    = 0.5 - self.params['gamma_norm'] / 2.0
        self.A[1][2]    = self.params['J_norm']
        self.A[2][1]    = - self.params['J_norm']
        self.A[2][2]    = - 0.5 - self.params['gamma_norm'] / 2.0
        self.A[3][0]    = self.params['J_norm']
        self.A[3][3]    = - 0.5 - self.params['gamma_norm'] / 2.0

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
        if not self.params['is_noisy']:
            return self.D

        # update noise matrix
        self.D[0][0]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[1][1]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[2][2]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[3][3]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)

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

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = 0.5 * np.cosh(2.0)
        iv_corrs[0][2]  = 0.5 * np.sinh(2.0)
        iv_corrs[1][1]  = 0.5 * np.cosh(2.0)
        iv_corrs[1][3]  = - 0.5 * np.sinh(2.0)
        iv_corrs[2][0]  = 0.5 * np.sinh(2.0)
        iv_corrs[2][2]  = 0.5 * np.cosh(2.0)
        iv_corrs[3][1]  = - 0.5 * np.sinh(2.0)
        iv_corrs[3][3]  = 0.5 * np.cosh(2.0)

        return None, iv_corrs, None
    
    def get_omega_norms(self, c):
        """Method to obtain the normalized eigenfrequencies.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns 
        -------
        omega_norms : *numpy.ndarray*
            Normalized eigenfrequencies.
        """

        # complex value
        temp = np.sqrt(self.params['J_norm']**2 - 0.25 + 0j)

        return np.array([temp, -temp], dtype=np.complex_)

class PhysRevA_100_063846_01(BaseSystem):
    r"""Class to simulate the Tripartite PT-symmetric QOM system in Phys. Rev. A **100**, 063846 (2019).
    
    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        gamma_norm  (*float*) normalized mechanical detuning rate :math:`\gamma / \Gamma`. Default is :math:`10^{-3}`.
        J_norm      (*float*) normalized coupling strength :math:`J / \Gamma`. Default is :math:`0.5`.
        n_th        (*float*) average thermal phonon occupancy :math:`n_{th}`. Default is :math:`0.0`.
        is_noisy    (*bool*) option to add or remove noises. Default is `True`.
        ========    ============================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'gamma_norm': 1e-3,
        'J_norm'    : 0.5,
        'n_th'      : 0.0,
        'is_noisy'  : True
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA_100_063846_01."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevA_100_063846_01',
            desc='Tripartite PT-symmetric QOM system in Phys. Rev. A 100, 063846',
            num_modes=3,
            cb_update=cb_update
        )

        # attributes
        self.is_A_constant = True

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
        
        # update drift matrix
        self.A[0][0]    = 0.5 - self.params['gamma_norm'] / 2.0
        self.A[0][3]    = - self.params['J_norm']
        self.A[1][1]    = 0.5 - self.params['gamma_norm'] / 2.0
        self.A[1][2]    = self.params['J_norm']
        self.A[2][1]    = - self.params['J_norm']
        self.A[2][2]    = - self.params['gamma_norm'] / 2.0
        self.A[2][5]    = - self.params['J_norm']
        self.A[3][0]    = self.params['J_norm']
        self.A[3][3]    = - self.params['gamma_norm'] / 2.0
        self.A[3][4]    = self.params['J_norm']
        self.A[4][3]    = - self.params['J_norm']
        self.A[4][4]    = - 0.5 - self.params['gamma_norm'] / 2.0
        self.A[5][2]    = self.params['J_norm']
        self.A[5][5]    = - 0.5 - self.params['gamma_norm'] / 2.0

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
        if not self.params['is_noisy']:
            return self.D

        # update noise matrix
        self.D[0][0]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[1][1]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[2][2]    = self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[3][3]    = self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[4][4]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)
        self.D[5][5]    = 0.5 + self.params['gamma_norm'] * (self.params['n_th'] + 0.5)

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

        # frequently used variables
        v_00 = np.exp(2.0) / 3.0 + 2 * np.exp(- 2.0) / 3.0
        v_02 = np.exp(2.0) / 3.0 - np.exp(- 2.0) / 3.0
        v_11 = np.exp(- 2.0) / 3.0 + 2 * np.exp(2.0) / 3.0
        v_13 = np.exp(- 2.0) / 3.0 - np.exp(2.0) / 3.0

        # initial values of the correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = 0.5 * v_00
        iv_corrs[0][2]  = 0.5 * v_02
        iv_corrs[0][4]  = 0.5 * v_02
        iv_corrs[1][1]  = 0.5 * v_11
        iv_corrs[1][3]  = 0.5 * v_13
        iv_corrs[1][5]  = 0.5 * v_13
        iv_corrs[2][0]  = 0.5 * v_02
        iv_corrs[2][2]  = 0.5 * v_00
        iv_corrs[2][4]  = 0.5 * v_02
        iv_corrs[3][1]  = 0.5 * v_13
        iv_corrs[3][3]  = 0.5 * v_11
        iv_corrs[3][5]  = 0.5 * v_13
        iv_corrs[4][0]  = 0.5 * v_02
        iv_corrs[4][2]  = 0.5 * v_02
        iv_corrs[4][4]  = 0.5 * v_00
        iv_corrs[5][1]  = 0.5 * v_13
        iv_corrs[5][3]  = 0.5 * v_13
        iv_corrs[5][5]  = 0.5 * v_11

        return None, iv_corrs, None
    
    def get_omega_norms(self, c):
        """Method to obtain the normalized eigenfrequencies.
        
        Parameters
        ----------
        c : *numpy.ndarray*
            Derived constants and controls.
        
        Returns 
        -------
        omega_norms : *numpy.ndarray*
            Normalized eigenfrequencies.
        """

        return np.roots([
            1,
            0.0,
            0.25 - 2 * self.params['J_norm']**2,
            0.0
        ])
