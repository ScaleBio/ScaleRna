"""
Function definitions for cell finder algorithm
"""
import numpy as np
import pandas as pd
from scipy.stats import dirichlet_multinomial, false_discovery_control
from scipy.sparse import sparray
from scipy.linalg import lstsq
from scipy.optimize import minimize

def interpolate_unseen_species(mtx: sparray, smoothed_probabilities: np.ndarray) -> np.ndarray:
    """
    (Subroutine of the Good-Turing algorithm used in construct_ambient_profile). For each gene not seen at all among
    ambient barcodes, interpolates the smoothed probability of encountering that gene among a count vector randomly
    sampled from the ambient profile.

    Args:
        mtx: Discretized STARsolo count matrix
        smoothed_probabilities: Array of smoothed probabilities for each gene

    Returns:
        An array of posterior expectations for the proportion of counts (summing to 1) assigned to each gene in the
        ambient profile
    """
    unseen_smoothed_probability_total = 1 - np.sum(smoothed_probabilities)
    unseen_gene_indices = np.where(smoothed_probabilities == 0.0)[0]
    unseen_gene_counts = mtx[unseen_gene_indices, :].sum(axis = 1).A
    unseen_genes_total_count = np.sum(unseen_gene_counts)
    unseen_smoothed_probabilities = np.multiply(unseen_gene_counts, unseen_smoothed_probability_total / unseen_genes_total_count)
    posterior_expectations = smoothed_probabilities
    posterior_expectations[unseen_gene_indices] = unseen_smoothed_probabilities.reshape(unseen_gene_indices.shape)

    return posterior_expectations

def construct_ambient_profile(discrete_mtx: sparray, ambient_barcode_indices: np.ndarray) -> np.ndarray:
    """
    (Subroutine of call_cells). Constructs a proportion vector representing the ambient profile.
    Applies the Good-Turing algorithm to the summed unique transcript count vector of all ambient barcodes
    to obtain a posterior expectation for the representation of each gene among those barcodes.

    Args:
        discrete_mtx: Discretized STARsolo count matrix
        ambient_barcode_indices: Array of indices for ambient barcodes

    Returns:
        An array of posterior expectations for the proportion of unique transcript counts (summing to 1)
        assigned to each gene in the ambient profile
    """
    # Sum counts of each gene across ambient barcodes
    gene_count_vector = discrete_mtx[:, ambient_barcode_indices].sum(axis = 1).A

    # Add a pseudo-feature of 1 to the summed counts to avoid zero probabilities being returned by the Good-Turing algorithm
    gene_count_vector = np.append(gene_count_vector, 1)

    # Dictionary mapping summed counts to the number of barcodes with those counts
    gene_count_counts = pd.value_counts(gene_count_vector).to_dict()
    del gene_count_counts[0.0]

    # Obtain smoothed probabilities for each gene across the ambient profile
    smoothed_probabilities = smooth_probabilities(gene_count_vector, gene_count_counts)

    # Remove pseudo-feature added above
    smoothed_probabilities = smoothed_probabilities[0:-1]

    # Interpolate smoothed probabilities for genes with zero counts in the ambient profile
    posterior_expectations = interpolate_unseen_species(discrete_mtx, smoothed_probabilities)

    return posterior_expectations

def scale_ambient_profile(discrete_mtx: sparray, ambient_profile: np.ndarray, ambient_barcode_indices: np.ndarray) -> tuple[np.ndarray, float]:
    """
    (Subroutine of test_ambiguous_barcodes). Scales the ambient profile by the maximum likelihood estimate of the scaling factor alpha. 

    Args:
        discrete_mtx: Discretized STARsolo count matrix
        ambient_profile: Array of posterior expectations for the proportion of counts (summing to 1) assigned to each gene in the ambient profile 
        ambient_barcode_indices: Array of indices for ambient barcodes

    Returns:
        The scaled ambient_profile, along with estimated alpha scaling factor by which the ambient profile was scaled
    """
    ambient_profile_count_vector = np.concatenate(discrete_mtx[:, ambient_barcode_indices].sum(axis = 1).A).astype(int)
    ambient_profile_total_counts = np.sum(ambient_profile_count_vector)

    # Estimate alpha by globally minimizing the negative log-likelihood function for the Dirichlet multinomial distribution using the limited-memory BFGS algorithm
    def negative_log_likelihood(alpha):
        return (-1) * dirichlet_multinomial.logpmf(x = ambient_profile_count_vector, alpha = np.multiply(ambient_profile, alpha), n = ambient_profile_total_counts)
    alpha = minimize(fun = negative_log_likelihood, x0 = 0.01, method = 'L-BFGS-B').x[0]

    # Scale the ambient profile by this estimate of alpha
    scaled_ambient_profile = np.multiply(ambient_profile, alpha)

    return scaled_ambient_profile, alpha

def rescue_cells(discrete_mtx: sparray, scaled_ambient_profile: np.ndarray, alpha: float, ambiguous_barcode_indices: np.ndarray) -> np.ndarray:
    """
    (Subroutine of test_ambiguous_barcodes). Tests each ambiguous barcode for deviation
    from the ambient profile. 10,000 count vectors are sampled from the ambient profile
    under the Dirichlet multinomial distribution, and the likelihood of each sampled
    ambient count vector is compared with the likelihood of each ambiguous barcode under
    the same distribution. A Monte Carlo p-value is computed for each ambiguous barcode
    over the 10,000 iterations.

    Args:
        discrete_mtx: Discretized STARsolo count matrix
        scaled_ambient_profile: Array of posterior expectations for the proportion of counts (summing to 1) assigned to each gene in the ambient profile 
        alpha: The estimated alpha scaling factor by which the ambient profile was scaled
        ambiguous_barcode_indices: Array of indices for ambiguous barcodes
        
    Returns:
        An array of Monte Carlo p-values associated with each ambiguous barcode.
        Each p-value in the array expresses the fraction of sampled 
        ambient count vectors whose likelihoods are <= the likelihood of the given 
        ambiguous barcode under the ambient model.
    """
    # Set a seed to ensure results here are reproducible for the same inputs
    np.random.seed(42)

    # Initialize array of p-values
    monte_carlo_p_values = np.array([1.0] * len(ambiguous_barcode_indices))

    # Obtain number of unique unique transcript counts among the ambiguous barcodes
    ambiguous_mtx = discrete_mtx[:, ambiguous_barcode_indices]
    total_count_by_ambiguous_barcode = np.asarray(ambiguous_mtx.sum(axis = 0)[0])[0]
    unique_total_counts_by_ambiguous_barcode = np.array(sorted(pd.value_counts(total_count_by_ambiguous_barcode).keys()))
    unique_total_count = int(unique_total_counts_by_ambiguous_barcode[0])

    # Wrapper function for the Dirichlet multinomial logPMF 
    def dirichlet_multinomial_logpmf(x):
        return dirichlet_multinomial.logpmf(x, alpha = scaled_ambient_profile, n = np.sum(x))
    
    # Calculate likelihood of each ambiguous barcode, splitting for memory efficiency
    ambiguous_barcode_likelihoods = []
    splits_of_ambiguous_barcodes = np.array_split(ambiguous_barcode_indices, np.ceil(len(ambiguous_barcode_indices) / 100000))
    for split_of_ambiguous_barcodes in splits_of_ambiguous_barcodes:
        split_mtx = discrete_mtx[:, split_of_ambiguous_barcodes].toarray()
        split_ambiguous_barcode_likelihoods = np.apply_along_axis(func1d = dirichlet_multinomial_logpmf, axis = 0, arr = split_mtx)
        ambiguous_barcode_likelihoods = np.hstack((ambiguous_barcode_likelihoods, split_ambiguous_barcode_likelihoods))

    # Sample 10,000 vectors from the ambient profile, each with a starting sum of the lowest total unique transcript count above minUTC
    dirichlet_probability_vectors = np.random.dirichlet(alpha = scaled_ambient_profile, size = 10000)
    def sample_from_dirichlet_prior(dirichlet_probability_vector):
        return np.random.multinomial(n = unique_total_count, pvals = dirichlet_probability_vector, size = 1)[0]
    ambient_vectors = np.apply_along_axis(func1d = sample_from_dirichlet_prior, axis = 1, arr = dirichlet_probability_vectors)

    # Compare the likelihoods of each of these 10,000 vectors to the likelihoods of each ambiguous barcode
    ambient_vector_likelihoods = np.apply_along_axis(func1d = dirichlet_multinomial_logpmf, axis = 1, arr = ambient_vectors)
    while unique_total_count <= np.max(unique_total_counts_by_ambiguous_barcode):
        if unique_total_count in total_count_by_ambiguous_barcode:
            for i in np.where(total_count_by_ambiguous_barcode == unique_total_count)[0]:
                R_1 = np.sum(ambient_vector_likelihoods <= ambiguous_barcode_likelihoods[i])
                if np.isnan(R_1):
                    R_1 = 0
                monte_carlo_p_values[i] = (monte_carlo_p_values[i] + R_1) / 10001
        unique_total_count += 1
        index_of_gene_to_increment = int(np.where(np.random.multinomial(n = 1, pvals = np.random.dirichlet(alpha = scaled_ambient_profile)) == 1)[0][0])
        alpha_of_gene_incremented = scaled_ambient_profile[index_of_gene_to_increment]

        # Update likelihoods for each next highest total unique transcript count by multiplying by closed form expression (see algorithm description)
        for ambient_vector_index, ambient_vector_likelihood in enumerate(ambient_vector_likelihoods):
            incremented_gene_count = ambient_vectors[ambient_vector_index][index_of_gene_to_increment] + 1
            ambient_vectors[ambient_vector_index][index_of_gene_to_increment] = incremented_gene_count
            ambient_vector_likelihoods[ambient_vector_index] = ambient_vector_likelihood + np.log((unique_total_count / (unique_total_count + alpha - 1)) * ((incremented_gene_count + alpha_of_gene_incremented - 1) / incremented_gene_count))
        
    return monte_carlo_p_values   

def test_ambiguous_barcodes(discrete_mtx: sparray, ambient_profile: np.ndarray, ambient_barcode_indices: np.ndarray, ambiguous_barcode_indices: np.ndarray) -> np.ndarray:
    """
    (Subroutine of call_cells). Tests all ambiguous barcodes with unique 
    transcript counts >= minUTC and < UTC for deviation from the ambient 
    profile. If barcodes deviate at the specified FDR from the ambient profile 
    under the Dirichlet multinomial distribution, they are rescued as cells; 

    Args:
        discrete_mtx: Discretized STARsolo count matrix
        ambient_profile: Array of posterior expectations for the proportion of counts (summing to 1) assigned to each gene in the ambient profile 
        ambient_barcode_indices: Array of indices for ambient barcodes 
        ambiguous_barcode_indices: Array of indices for ambiguous barcodes 

    Returns:
        An array of Monte Carlo p-values indicating each ambiguous barcode's degree of deviation from the ambient profile,
        corrected for multiple testing by the Benjamini-Hochberg procedure
    """
    FDRs = []

    if len(ambiguous_barcode_indices) > 0:

        # Scale the ambient profile by the estimated scaling factor alpha for the Dirichlet multinomial distribution
        scaled_ambient_profile, alpha = scale_ambient_profile(discrete_mtx, ambient_profile, ambient_barcode_indices)

        # Test each ambiguous barcode for deviation from the ambient profile, obtaining a Monte Carlo p-value for each
        monte_carlo_p_values = rescue_cells(discrete_mtx, scaled_ambient_profile, alpha, ambiguous_barcode_indices)

        # Adjust p-values for multiple testing using the Benjamini-Hochberg procedure
        FDRs = false_discovery_control(monte_carlo_p_values, method = 'bh')
    
    return FDRs

def interpolate_gene_counts(gene_count_counts: dict[float:int]) -> dict[int:float]:
    """
    (Subroutine of the Good-Turing algorithm used in construct_ambient_profile).
    For each count b, obtains the linear interpolation of b given a, c,
    where a is the greatest count < b, and c is the smallest count > b.

    Args: 
        gene_count_counts: Dictionary mapping counts to the number of genes with that count

    Returns: 
        A dictionary mapping counts to their linear interpolations
    """
    linear_interpolations = {}
    gene_counts = sorted(gene_count_counts.keys())
    for (b_index, b_value) in enumerate(gene_counts):
        if b_index == 0:
            a = 0
        else:
            a = gene_counts[b_index - 1]
        if b_index == len(gene_counts) - 1:
            c = 2 * b_value - a
        else:
            c = gene_counts[b_index + 1]
        linear_interpolations[b_value] = 2 * gene_count_counts[b_value] / float(c - a)

    return linear_interpolations

def smooth_probabilities(gene_count_vector: np.ndarray, gene_count_counts: dict[float:int]) -> np.ndarray:
    """
    (Subroutine of the Good-Turing algorithm used in construct_ambient_profile). For each gene, calculates the smoothed probability of encountering that gene among ambient barcodes given the contingent possibility of encountering genes not seen at all among those barcodes.

    Args:
        gene_count_vector: Array of counts per gene
        gene_count_counts: Dictionary mapping counts to the number of genes with those counts

    Returns:
        An array of smoothed probabilities for each gene
    """
    # Linearly interpolate counts for each gene
    interpolations = interpolate_gene_counts(gene_count_counts)

    # Perform loglinear least-squares regression
    a, b = lstsq(np.c_[np.log(np.array(list(interpolations.keys()))), (1,) * len(np.array(list(interpolations.keys())))], np.log(np.array(list(interpolations.values()))))[0]

    # Perform loglinear smoothing
    smoothed_gene_counts = {}
    use_smoothing = False
    for gene_count in sorted(gene_count_counts.keys()):
        y = float(gene_count + 1) * np.exp(a * np.log(gene_count + 1) + b) / np.exp(a * np.log(gene_count) + b)
        next_higher_gene_count = gene_count + 1
        if next_higher_gene_count not in gene_count_counts:
            use_smoothing = True
        if use_smoothing:
            smoothed_gene_counts[gene_count] = y
            continue
        turing_estimate = (float(next_higher_gene_count) * gene_count_counts[next_higher_gene_count]) / gene_count_counts[gene_count]
        Nr = float(gene_count_counts[gene_count])
        Nr_next_higher = float(gene_count_counts[next_higher_gene_count])
        confidence_interval_width = 0.95 * np.sqrt(float(next_higher_gene_count)**2 * (Nr_next_higher / Nr**2) * (1.0 + (Nr_next_higher / Nr)))
        if abs(turing_estimate - y) > confidence_interval_width:
            smoothed_gene_counts[gene_count] = turing_estimate
        else:
            smoothed_gene_counts[gene_count] = y
        use_smoothing = True

    # Normalize 
    total_smoothing = 0.0
    for gene_count, smoothed_value in smoothed_gene_counts.items():
        total_smoothing += gene_count_counts[gene_count] * smoothed_value  
    p0 = gene_count_counts[1.0] / sum(gene_count_vector)  
    smoothed_probabilities = np.array([((1.0 - p0) * (smoothed_gene_counts[i] / total_smoothing)) if i in smoothed_gene_counts else i * 0.0 for i in gene_count_vector])

    return smoothed_probabilities
