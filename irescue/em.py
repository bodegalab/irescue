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

def log_likelihood(matrix, counts):
    """
    Compute log-likelihood of data.
    """
    likelihoods = (matrix * counts).sum(axis=1)
    log_likelihood = np.sum(np.log(likelihoods + np.finfo(float).eps))
    return log_likelihood

def run_em(matrix, cycles=100, tolerance=1e-5):
    """
    Run Expectation-Maximization (EM) algorithm to redistribute read counts
    across a set of features.

    Parameters
    ----------
    matrix : array
        Reads-features compatibility matrix.
    cycles : int, optional
        Number of EM cycles.
    tolerance : float
        Tolerance threshold of log-likelihood difference to infer convergence.

    Returns
    -------
    out : list
        Optimized relative feature abundances.
    cycle : int
        Number of EM cycles.
    converged : bool
        Indicates if convergence has been reached before cycles threshold.
    """

    # calculate initial estimation of relative abundance.
    # (let the sum of counts of features be 1,
    # will be multiplied by the real UMI count later)
    nFeatures = matrix.shape[1]
    counts = np.array([1 / nFeatures] * nFeatures)

    # Initial log-likelihood
    prev_loglik = log_likelihood(matrix, counts)

    converged = False
    curr_cycle = 0

    # Run EM iterations
    while curr_cycle < cycles:
        curr_cycle += 1
        e_matrix = e_step(matrix=matrix, counts=counts)
        counts = m_step(matrix=e_matrix)

        # Compute the new log-likelihood
        loglik = log_likelihood(matrix, counts)

        # Check for convergence
        loglikdiff = loglik - prev_loglik
        if np.abs(loglikdiff) < tolerance:
            converged = True
            break

        prev_loglik = loglik

    return counts, (curr_cycle, converged, loglik, loglikdiff)
