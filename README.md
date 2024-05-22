# EEL_correction
Code used in the effective experiment length correction manuscript

# Description 
R function to adjust condition-specific changes in fitness values, removing any correlation between reference fitness and condition-specific changes in fitness which may have occurred from a different number of library doublings between experiments.

This function breaks the reference and condition-specific changes in fitnesses of a mutant library into bins of a specified size, in order of increasing reference fitness. For each bin, the median of the reference and condition-specific changes in fitnesses is taken. The difference in medians is added to all condition-specific changes in fitnesses within that bin such that the median value is the same as reference for all mutants in that bin. Done over the whole library, this removes the correlation effect which occurs from different effective experimental lengths. The function returns the corrected condition-specific changes in fitnesses and a summary table for each bin.

# Code Usage
The input format is ```correct_l2fc_v3(initial_gammas, other_gammas, window_number=10000, norm_to_constant=NA)```

```initial_gammas``` and ```other_gammas``` are vectors containing the reference (```initial_gammas```) or condition-specific (```other_gammas```) fitnesses for all sgRNAs/mutants, assumed to be in the same order.

```window_number``` is the size of the bins. The number of bins depends on the size of the library. For Tn-seq data, 50 is likely appropriate, but for CRISPRi data, a larger bin size (200,1000) can be used, as there are more growth-deficient strains. default = 1000

```norm_to_constant``` this param depends on how the condition-specific L2FC was calculated. If both the condition-specific and reference fitness L2FCs were calculated compared to a t0 sample, use ```NA``` (default). If only the reference fitness L2FC was calculated compared to a t0, and the condition-specific L2FC was calculated compared to the reference endpoint, this value is the constant value (usually 0). 

For example, say a sick mutant has a reference fitness of -10, calculated by doing L2FC to a t0 sample. In this example we expect the mutant to have the same fitness value in the reference and condition-specific experiments. If the fitness in the reference condition was calculated to a t0 sample as well, it would be expected to have a condition-specific changes in fitness of -10 in the drug condition. However, if the condition-specific fitness was calculated L2FC to the reference, it would expected to be 0. This difference requires two different approaches in how the condition-specific changes in fitness values are corrected depending on the calculation method.

The function outputs a list with two values. First is a matrix summarizing the bins. Tor each bin of specified window length (rows), the columns contain:
1. the median of the reference fitnesses 
2. median of the condition-specific changes in fitnesses
3. the median difference between the reference and condition-specific changes in fitnesses 
4. the number of sgRNAs/mutants which have a non-NA value.
5. median average distance (MAD) of the difference
The second value is a vector containing the corrected condition-specific changes in fitness values, in the same order as the original other_gammas input vector.

In cases where there is more than one experiment to be corrected (ie a large-scale chemical genomics screen with many conditions tested and multiple replicates), it would be prudent to run in a loop against all conditions and consolidate all corrected vectors into a single corrected table.

Psuedocode:
```all_experiments <- read.csv("allexperimentsfitnesses.csv") 
reference <- all_experiments$reference
for (all columns in all_experiments) :
      curr_results <- correct_l2fc_v3(reference, current experiment, window_number = 10000, norm_to_constant = NA)
      curr_corrected <- curr_results[[2]]
      corrected_results[,current experiment] <- curr_corrected

