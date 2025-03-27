#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate the cavity soliton system in Phys. Rev. A **100**, 053814 (2019)."""

__authors__ = ["Sampreet Kalita"]
__toolbox__ = "qom-v1.1.0"
__created__ = "2022-07-25"
__updated__ = "2025-03-11"
__all__     = ['PhysRevA_100_053814']

# dependencies
import numpy as np

# qom modules
from qom.systems import BaseSystem

class PhysRevA_100_053814(BaseSystem):
    r"""Class to simulate the cavity soliton system in Phys. Rev. A **100**, 053814 (2019).

    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ========    ============================================================
        key         meaning
        ========    ============================================================
        N           (*int*) number of divisions in fast time :math:`N`. Default is :math:`2^{12}`.
        Delta       (*float*) normalized detuing of the cavity :math:`\Delta`. Default is :math:`3.0`.
        S_0         (*float*) normalized value of the source term :math:`S_{0}`. Default is :math:`\sqrt{3.5}`.
        tau_max     (*float*) maximum value of fast time :math:`\tau_{max}`. Default is :math:`100.0`.
        t_alphas    (*str*) type of initial occupancies `t_alphas`. Options are `"sech"` for sec-hyperbolic function and `"none"` for zero initial occupancies. Default is `"sech"`.
        ========    ============================================================
    cb_update : *callable*, optional
        Callback function to update status and progress, formatted as `cb_update(status, progress, reset)`, where `status` is a string, `progress` is a float and `reset` is a boolean.
    """

    # default parameters of the system
    system_defaults = {
        'N'         : 2**12,
        'Delta'     : 3.0,
        'S_0'       : np.sqrt(3.5),
        'tau_max'   : 100.0,
        't_alphas'  : 'sech'
    }

    def __init__(self, params, cb_update=None):
        """Class constructor for PhysRevA_100_053814."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='PhysRevA_100_053814',
            desc="Cavity Soliton System in Phys. Rev. A 100, 053814",
            num_modes=2 * params.get('N', self.system_defaults['N']),
            cb_update=cb_update
        )

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
        N           = self.params['N']
        tau_max     = self.params['tau_max']
        t_alphas    = self.params['t_alphas']
 
        # initial values of the modes
        iv_modes = np.zeros(self.num_modes, dtype=np.complex128)
        if t_alphas == 'sech':
            iv_modes[::2] = 0.5 + 1.0 / np.cosh([(i - N / 2.0) * 2.0 * tau_max / N for i in range(N)])

        return iv_modes, None, None

    def get_coeffs_dispersion(self, modes, c, t):
        r"""Method to get the coefficients of :math:``( i \omega )^{j}`` in the dispersions.
        
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
        coeffs : *numpy.ndarray*
            Coefficients in the dispersion operator.
        """

        # extract frequently used variables
        N       = self.params['N']
        delta_N = N / 2.0 / self.params['tau_max']

        # return coefficients
        return np.array([0.0j, 0.0j, - 1.0j * (-1.0) * delta_N**2], dtype=np.complex128)

    def get_nonlinearities(self, modes, c, t):
        """Method to get the nonlinearities.
        
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
        nonlinearities : *numpy.ndarray*
            Nonlinearities.
        """

        # return nonlinearities
        return - 1.0 + 1.0j * (np.abs(modes[::2])**2 - self.params['Delta'])

    def get_sources(self, modes, c, t):
        """Method to get the source terms.
        
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
        sources : *numpy.ndarray*
            Source terms.
        """

        return self.params['S_0']


