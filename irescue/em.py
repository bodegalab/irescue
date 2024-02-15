import numpy as np

def e_step(matrix, counts):
    """
    Performs E-step of EM algorithm: proportionally assigns reads to features
    based on relative feature abundances.
    """
    colsums = (matrix * counts).sum(axis=1)[:, np.newaxis]
    out = matrix / colsums * counts
    return(out)

def m_step(matrix):
    """
    Performs M-step of EM algorithm: calculates feature abundances from read
    counts proportionally distributed to features.
    """
    counts = matrix.sum(axis=0) / matrix.sum()
    return(counts)

def run_em(matrix, cycles=100):
    """
    Run Expectation-Maximization (EM) algorithm to redistribute read counts
    across a set of features.

    Parameters
    ----------
    matrix : array
        Reads-features compatibility matrix.
    cycles : int, optional
        Number of EM cycles.

    Returns
    -------
    out : list
        Optimized relative feature abundances.
    """
    
    # calculate initial estimation of relative abundance.
    # (let the sum of counts of features be 1,
    # will be multiplied by the real UMI count later)
    nFeatures = matrix.shape[1]
    counts = np.array([1 / nFeatures] * nFeatures)

    # run EM for n cycles
    for _ in range(cycles):
        e_matrix = e_step(matrix=matrix, counts=counts)
        counts = m_step(matrix=e_matrix)

    return(counts)