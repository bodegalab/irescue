import numpy as np

def e_step(matrix, counts):
    """
    Performs E-step of EM algorithm: proportionally assigns reads to features
    based on relative feature abundances.
    """
    rows, cols = matrix.shape
    out = np.empty(matrix.shape, dtype=float)
    for row in range(rows):
        for col in range(cols):
            # each matrix value is divided by the sum of its column and
            # multiplied by the current feature value
            out[row,col] = (
                matrix[row,col]
                / sum([counts[i] for i in np.where(matrix[row] > 0)[0]])
                * counts[col]
            )
    return(out)

def m_step(matrix):
    """
    Performs M-step of EM algorithm: calculates feature abundances from read
    counts proportionally distributed to features.
    """
    nFeatures = matrix.shape[1]
    counts = list()
    for i in range(nFeatures):
        counts.append(matrix[:, i].sum() / matrix.sum())
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
    # let the sum of counts of features be 1
    # (it will be multiplied by the real count later)
    nFeatures = matrix.shape[1]
    counts = [1 / nFeatures] * nFeatures

    # run EM for n cycles
    for _ in range(cycles):
        e_matrix = e_step(matrix=matrix, counts=counts)
        counts = m_step(matrix=e_matrix)

    return(counts)