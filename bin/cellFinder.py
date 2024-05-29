"""
Implementation of a statistical method to classify cell and background barcodes, following Lun et al. (2019)
"""
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import dirichlet_multinomial, false_discovery_control
from scipy import stats
from scipy.sparse import sparray
from scipy.linalg import lstsq
from scipy import optimize

# Set a seed to ensure results here are reproducible for the same inputs
RNG = np.random.default_rng(42)

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
    unseen_gene_counts = mtx[unseen_gene_indices, :].sum(axis=1)
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
    gene_count_vector = discrete_mtx[:, ambient_barcode_indices].sum(axis=1)

    # Add a pseudo-feature of 1 to the summed counts to avoid zero probabilities being returned by the Good-Turing algorithm
    gene_count_vector = np.append(gene_count_vector, 1)

    # Dictionary mapping summed counts to the number of barcodes with those counts
    gene_count_counts = pd.Series(gene_count_vector).value_counts().to_dict()
    del gene_count_counts[0.0]

    # Obtain smoothed probabilities for each gene across the ambient profile
    smoothed_probabilities = smooth_probabilities(gene_count_vector, gene_count_counts)

    # Remove pseudo-feature added above
    smoothed_probabilities = smoothed_probabilities[0:-1]

    # Interpolate smoothed probabilities for genes with zero counts in the ambient profile
    posterior_expectations = interpolate_unseen_species(discrete_mtx, smoothed_probabilities)

    return posterior_expectations

def estimate_alpha(discrete_mtx: sparray, ambient_profile: np.ndarray, ambient_barcode_indices: np.ndarray) -> float:
    """
    ML estimate of the overdisperion (dirichlet scaling factor alpha) from a sample of ambient_barcodes. 

    Args:
        discrete_mtx: Discretized count matrix
        ambient_profile: Array of posterior expectations for the proportion of counts (summing to 1) assigned to each gene in the ambient profile 
        ambient_barcode_indices: Array of indices for ambient barcodes

    Returns:
        Estimated alpha scaling factor
    """
    def negative_log_likelihood(alpha, x, ns):
        lls = dirichlet_multinomial.logpmf(x=x, alpha=alpha*ambient_profile, n=ns)
        return -1 * stats.trim_mean(lls, proportiontocut=0.001)
    min_count = 10 # Exclude barcodes with lower total counts
    barcode_inds = np.setdiff1d(ambient_barcode_indices, np.where(discrete_mtx.sum(axis = 0) < min_count))
    num_vecs = min(1000, len(barcode_inds)) # Number of ambient barcodes to sample for optimization
    sampled_vecs = RNG.choice(barcode_inds, num_vecs, replace=False)
    vec = discrete_mtx[:, sampled_vecs].toarray().T
    ns = vec.sum(1)
    res = optimize.minimize_scalar(method='bounded', fun = negative_log_likelihood, args=(vec, ns), bounds=(1, 10000), options={'xatol':0.001})
    if not res.success: 
        print(f"Alpha (overdispersion) optimization failed: {res.message}", file=sys.stderr)
        return 3000 # Default alpha
    alpha = res.x
    print(alpha, negative_log_likelihood(alpha, vec, ns))
    return alpha

def rescue_cells(discrete_mtx: sparray, scaled_ambient_profile: np.ndarray, alpha: float, ambiguous_barcode_indices: np.ndarray, mcvecs: int = 10000) -> np.ndarray:
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
        mcvecs: Number of Monte Carlo vectors to sample from the ambient profile
        
    Returns:
        An array of Monte Carlo p-values associated with each ambiguous barcode.
        Each p-value in the array expresses the fraction of sampled 
        ambient count vectors whose likelihoods are <= the likelihood of the given 
        ambiguous barcode under the ambient model.
    """
    # Initialize array of p-values
    monte_carlo_p_values = np.array([1.0] * len(ambiguous_barcode_indices))

    # Obtain number of unique unique transcript counts among the ambiguous barcodes
    ambiguous_mtx = discrete_mtx[:, ambiguous_barcode_indices]
    total_count_by_ambiguous_barcode = pd.Series(ambiguous_mtx.sum(axis=0))
    unique_total_counts_by_ambiguous_barcode = np.array(sorted(total_count_by_ambiguous_barcode.value_counts().keys()))
    unique_total_count = int(unique_total_counts_by_ambiguous_barcode[0])

    # Calculate likelihood of all ambiguous barcode, split into dense matrix chunks for memory efficiency
    ambiguous_barcode_likelihoods = []
    chunk_size = 10000
    barcode_chunks = np.array_split(ambiguous_barcode_indices, np.ceil(len(ambiguous_barcode_indices) / chunk_size))
    for chunk in barcode_chunks:
        split_mtx = discrete_mtx[:, chunk].toarray().transpose()
        lls = dirichlet_multinomial.logpmf(split_mtx, alpha=np.array([scaled_ambient_profile]), n=split_mtx.sum(1))
        ambiguous_barcode_likelihoods.append(lls)
    ambiguous_barcode_likelihoods = np.concatenate(ambiguous_barcode_likelihoods)

    # Sample 10,000 vectors from the ambient profile, each with a starting sum of the lowest total unique transcript count above minUTC
    dirichlet_probability_vectors = RNG.dirichlet(alpha = scaled_ambient_profile, size = mcvecs)
    ambient_vectors = RNG.multinomial(unique_total_count, dirichlet_probability_vectors)

    # Compare the likelihoods of each of these 10,000 vectors to the likelihoods of each ambiguous barcode
    ambient_vector_likelihoods = stats.dirichlet_multinomial.logpmf(ambient_vectors, [scaled_ambient_profile], ambient_vectors.sum(1))
    max_count = int(np.max(unique_total_counts_by_ambiguous_barcode))
    # Faster to call random once here for all increments for each step and random vector below 
    incrGenes = np.array([RNG.choice(len(scaled_ambient_profile), max_count+1, p=dirichlet_probability_vectors[i]) for i in range(mcvecs)])
    while unique_total_count <= max_count:
        for i in np.where(total_count_by_ambiguous_barcode == unique_total_count)[0]:
            R_1 = np.sum(ambient_vector_likelihoods <= ambiguous_barcode_likelihoods[i])
            monte_carlo_p_values[i] = (R_1+1) / (mcvecs+1)

        unique_total_count += 1
        # Update likelihoods for each next highest total unique transcript count by multiplying by closed form expression (see algorithm description)
        for ambient_vector_index, ambient_vector_likelihood in enumerate(ambient_vector_likelihoods):
            index_of_gene_to_increment = incrGenes[ambient_vector_index][unique_total_count-1]
            alpha_of_gene_incremented = scaled_ambient_profile[index_of_gene_to_increment]
            incremented_gene_count = ambient_vectors[ambient_vector_index][index_of_gene_to_increment] + 1
            ambient_vectors[ambient_vector_index][index_of_gene_to_increment] = incremented_gene_count
            ambient_vector_likelihoods[ambient_vector_index] += np.log((unique_total_count / (unique_total_count + alpha - 1)) * ((incremented_gene_count + alpha_of_gene_incremented - 1) / incremented_gene_count))
    
    return monte_carlo_p_values   

def test_ambiguous_barcodes(discrete_mtx: sparray, ambient_profile: np.ndarray, ambient_barcode_indices: np.ndarray, ambiguous_barcode_indices: np.ndarray, alpha: Optional[float]=0) -> np.ndarray:
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
        if not alpha:
            alpha = estimate_alpha(discrete_mtx, ambient_profile, ambient_barcode_indices)
        scaled_ambient_profile = alpha * ambient_profile

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
