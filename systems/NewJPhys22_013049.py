#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the systems described in New. J. Phys. **22**, 013049 (2020)."""

__authors__ = ['Sampreet Kalita']
__created__ = '2021-07-27'
__updated__ = '2021-09-03'
__version__ = '0.8.0'

# dependencies
import numpy as np

# qom modules
from qom.systems import SOSMSystem

class NewJPhys22_013049(SOSMSystem):
    """Class to simulate the OM system in New. J. Phys. **22**, 013049 (2020).

    Parameters
    ----------
    params : dict
        Parameters for the system.
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.
    """

    def __init__(self, params, cb_update=None):
        """Class constructor for NewJPhys22_013049."""
        
        # initialize super class
        super().__init__(params=params, cb_update=cb_update)

        # set attributes
        self.code = 'new_j_phys_22_013049'
        self.name = 'OM System New. J. Phys. **22**, 013049 (2020)'  
        # default parameters
        self.params = {
            'Delta_norm': params.get('Delta_norm', 0.0),
            'gamma_norm': params.get('gamma_norm', 1e-4),
            'kappa_norm': params.get('kappa_norm', 0.1),
            'P': params.get('P', 0.1),
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
            Constant parameters.
        t : float
            Time at which the drift matrix is calculated.
        
        Returns
        -------
        A : list
            Drift matrix.
        """
        
        # extract frequently used variables
        Delta_norm  = params[0]
        gamma_norm  = params[1]
        kappa_norm  = params[2]
        P           = params[3]
        alpha, beta = [modes[0], modes[1]]
        dim         = (2 * self.num_modes, 2 * self.num_modes)
        
        # drift matrix
        if self.A is None or np.shape(self.A) != dim:
            self.A = np.zeros(dim, dtype=np.float_)
        # optical mode
        self.A[0][0] = - kappa_norm / 2
        self.A[0][1] = - Delta_norm - 2 * np.real(beta)
        self.A[0][2] = - 2 * np.imag(alpha)
        self.A[1][0] = Delta_norm + 2 * np.real(beta)
        self.A[1][1] = - kappa_norm / 2
        self.A[1][2] = 2 * np.real(alpha)
        # mechanical mode
        self.A[2][2] = - gamma_norm / 2
        self.A[2][3] = 1.0
        self.A[3][0] = P * np.real(alpha)
        self.A[3][1] = P * np.imag(alpha)
        self.A[3][2] = - 1.0
        self.A[3][3] = - gamma_norm / 2

        return self.A

    def get_jac(self, modes, params, t):
        """Function to obtain the jacobian matrix.

        Parameters
        ----------
        modes : list
            Values of the modes.
        params : list
            Constant parameters.
        t : float
            Time at which the rates are calculated.
        
        Returns
        -------
        jac : list
            Jacobian matrix.
        """

        # extract frequently used variables
        Delta_norm  = params[0]
        gamma_norm  = params[1]
        kappa_norm  = params[2]
        P           = params[3]
        alpha, beta = [modes[0], modes[1]]

        # optical mode
        jac = [ [1j * (Delta_norm + 2 * np.real(beta)) - kappa_norm / 2, 1j * alpha],
                [1j * P / 2 * np.conjugate(alpha), - 1j - gamma_norm / 2]]

        return np.array(jac)

    def get_ivc(self):
        """Function to obtain the initial values and constants required for the IVP.
        
        Returns
        -------
        iv : list
            Initial values of variables.
            First element contains the optical mode.
            Next element contains the mechanical mode.
            Next (2*2)^2 elements contain the correlations.

        c : list
            Constant parameters.
            First (2*2)^2 elements contain the noise matrix.
            Next element contains the normalized laser detuning.
            Next element contains the normalized mechanical decay rate.
            Next element contains the normalized optical decay rate.
            Next element contains the normalized coupling.
        """

        # extract frequently used variables
        Delta_norm  = self.params['Delta_norm']
        gamma_norm  = self.params['gamma_norm']
        kappa_norm  = self.params['kappa_norm']
        P           = self.params['P']
        dim         = (2 * self.num_modes, 2 * self.num_modes)
 
        # initial mode values as 1D list
        modes_0 = np.zeros(self.num_modes, dtype=np.complex_).tolist()

        # initial quadrature correlations
        corrs_0 = np.zeros(dim, dtype=np.float_)
        corrs_0[0][0] = 0.5
        corrs_0[1][1] = 0.5
        corrs_0[2][2] = 0.5
        corrs_0[3][3] = 0.5

        # convert to 1D list and concatenate all variables
        iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

        # noise correlation matrix
        D = np.zeros(dim, dtype=np.float_)
        for i in range(2):
            D[0][0] = kappa_norm
            D[1][1] = kappa_norm
            D[2][2] = gamma_norm
            D[3][3] = gamma_norm
        
        # constant parameters
        params = [Delta_norm, gamma_norm, kappa_norm, P]

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
        Delta_norm  = params[0]
        gamma_norm  = params[1]
        kappa_norm  = params[2]
        P           = params[3]
        alpha, beta = [modes[0], modes[1]]

        # initialize lists
        dalpha_dt = (1j * Delta_norm - kappa_norm / 2) * alpha + 1j * alpha * (beta + np.conjugate(beta)) + 1 / 2
        dbeta_dt = (- 1j - gamma_norm / 2) * beta + 1j * P / 2 * np.conjugate(alpha) * alpha

        # rearrange per system
        mode_rates = [dalpha_dt, dbeta_dt]

        return mode_rates

    def get_oss_args(self, params):
        r"""Method to obtain the arguments required to calculate the optical steady state.
        
        Parameters
        ----------
        params : list
            Constant parameters of the system.
        
        Returns 
        -------
        A_l : float
            Amplitude of the laser.
        Delta : float
            Effective detuning if method is "basic", else detuning of the laser.
        kappa : float
            Optical decay rate.
        C : float
            Coefficient of :math:`|\alpha_{s}|^{2}`.
        """

        # extract frequently used variables
        Delta_norm  = params[0]
        gamma_norm  = params[1]
        kappa_norm  = params[2]
        P           = params[3]
        
        # drive amplitude
        A_l = 0.5
        # Coefficient of the mean optical occupancies
        C = 4 * P / (gamma_norm**2 + 4)

        return 0.5, Delta_norm, kappa_norm, C

    def get_ss_modes(self, params):
        """Method to obtain the steady state optical and mechanical mode apmlitudes.
        
        Parameters
        ----------
        params : list
            Constant parameters of the system.
        
        Returns 
        -------
        modes : list
            Optical and mechanical mode amplitudes.
        """
        
        # frequently used variables
        Delta_norm, gamma_norm, kappa_norm, P = params

        # optical mode occupancy
        N_o, _ = self.get_mean_optical_occupancies()
        n_o = N_o[0] if len(N_o) == 1 else N_o[1]

        # mechanical mode position
        beta_real = 2 * P / (gamma_norm**2 + 4) * n_o

        # mode amplitudes
        alpha = 1 / (kappa_norm - 2j * (Delta_norm + 2 * beta_real))
        beta = P * n_o * (2 + 1j * gamma_norm) / (gamma_norm**2 + 4)

        # list of mode amplitudes
        modes = [alpha, beta]

        return modes


