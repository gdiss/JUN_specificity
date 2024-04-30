
"""
Custom function
"""

import torch

def TwoStateFractionFoldedMod(
    X = None,
    trainable_parameters = {}):
    """
    1-dimensional nonlinear transformation relating modified Gibbs free energy of folding to fraction of molecules folded.

    :param X: list of tensors (required).
    :param trainable_parameters: dictionary of global parameter names (optional).
    :returns: fraction of molecules folded tensor.
    """  
    if X is None:
        return {}
    else:
        return torch.pow(1+torch.exp(X[0]*torch.abs(X[1])), -1)
