#' Correction of condition-specific changes in fitness values.
#' 
#' `correct_l2fc_v3` applies a correction to a vector of condition-specific changes in fitnesses which ablates any correlation between condition-specific changes in and reference fitnesses which
#' may have arisen from different effective experimental lengths.
#'
#' This function breaks the reference and condition-specific changes in fitnesses of a mutant library into bins of a specified size, in order of increasing reference fitness. For each bin, 
#' the median of the reference and condition-specific changes in fitnesses is taken. The difference in medians is added to all condition-specific changes in fitnesses within that bin such that the median value
#' is the same as reference for all mutants in that bin. Done over the whole library, this removes the correlation effect which occurs from different effective 
#' experimental lengths. The function returns the corrected condition-specific changes in fitnesses, a summary table for each bin, and two lists of initial values in each bin.
#' 
#' @param initial_gammas a vector containing the reference fitnesses for all sgRNAs/mutants. Order of mutants assumed to be same as other_gammas.
#' 
#' @param other_gammas a vector containing the condition-specific changes in fitnesses for sgRNAs/mutants to be corrected for a single experiment. Order of mutants assumed to be same as initial_gammas.
#' 
#' @param window_number the size of the bins. The number of bins depends on the size of the library. For Tn-seq data, 50 is likely appropriate, but for CRISPRi data,
#' a larger bin size (200,1000) can be used, as there are more growth-deficient strains. default = 1000
#' 
#' @param norm_to_constant this param depends on how the condition-specific L2FC was calculated. If both the condition-specific 
#' and reference fitness L2FCs were calculated compared to a t0 sample, use NA (default). If only the reference fitness L2FC was
#' calculated compared to a t0, and the condition-specific L2FC was calculated compared to the reference, this value is the constant value (usually 0).
#' For example, say a sick mutant has a reference fitness of -10, calculated by doing L2FC to a t0 sample. In this example
#' we expect the mutant to have the same fitness value in the reference and condition-specific experiments.
#' If the fitness in the reference condition was calculated to a t0 sample as well, it would be expected to have a condition-specific changes in fitness of
#' -10 in the drug condition. However, if the condition-specific fitness was calculated L2FC to the reference, it would expected to be 0. This
#' difference requires two different approaches in how the condition-specific changes in fitness values are corrected depending on the calculation method.
#' 
#' @returns a list with two values. 
#' The first value of the list is a matrix: for each bin of specified window length (rows), the columns contain:
#' 1. the median of the reference fitnesses 
#' 2. median of the condition-specific changes in fitnesses
#' 3. the median difference between the reference and condition-specific changes in fitnesses 
#' 4. the number of sgRNAs/mutants which have a non-NA value.
#' 5. median average distance (MAD) of the difference
#' The second value is a vector containing the corrected condition-specific changes in fitness values, in the same order as the original other_gammas 
#' input vector.
#' 
#' @author Horia Todor and Lili Kim
#'

correct_l2fc_v3 <- function(initial_gammas, other_gammas, window_number = 1000, norm_to_constant = NA){
  
  #initialize vector of corrected condition-specific changes in fitnesses to be returned
  corrected_other <- rep(NA, length(other_gammas))
  
  #only select mutant fitnesses which are not NA in either the reference or condition-specific experiments; NA values are unmeasured and preclude reliable correction
  not_na <- which(!is.na(initial_gammas) & !is.na(other_gammas))
  initial_gammas <- initial_gammas[not_na]
  other_gammas <- other_gammas[not_na]
  
  #get the order of input reference fitnesses by lowest to highest for later binning and accessing
  initial_gamma_order <- order(initial_gammas)
  
  #initialize results table based on window and input size
  results_size <- ifelse(length(initial_gammas) %% window_number == 0, 
                         length(initial_gammas) %/% window_number, 
                         length(initial_gammas) %/% window_number+1)
  results <- matrix(NA, results_size, 5)
  colnames(results) <- c("initial median", "other median", "difference", "number of not NA sgRNAs", "other mad")
  
  list_of_bins <- list()
  list_of_bins2 <- list()
  
  #for each bin of specified window size, we report the median reference and condition-specific changes in fitnesses, determine and apply the correction, and report the summary.
  for (i in 1:dim(results)[1]){
    
    #set the current bin
    start <- ((i-1)*window_number)+1
    stop <- min(c(length(initial_gammas), (i*window_number)))
    
    #get the medians of the reference condition-specific changes in fitnesses in the current bin and report in the summary table.
    #The bins are populated by increasing reference fitness 
    results[i,1] <- median(initial_gammas[initial_gamma_order[start:stop]], na.rm = TRUE)
    results[i,2] <- median(other_gammas[initial_gamma_order[start:stop]], na.rm = TRUE)
    list_of_bins[[i]] <- other_gammas[initial_gamma_order[start:stop]]
    list_of_bins2[[i]] <- initial_gammas[initial_gamma_order[start:stop]]
    
    #Here we determine the correction. The method depends on how condition-specific changes in fitness was calculated
    #If condition-specific changes in fitness was calculated by taking L2FC to t0, the correction is the median difference between the reference and condition-specific changes in fitnesses
    if (is.na(norm_to_constant)){
      differences <- other_gammas[initial_gamma_order[start:stop]]-initial_gammas[initial_gamma_order[start:stop]]
      results[i,3] <- median(differences, na.rm = TRUE)
      results[i,5] <- mad(differences, na.rm = TRUE)

    } 
    #If condition-specific changes in fitness was calculated by taking L2FC to reference, the correction is the median difference between condition-specific changes in fitnesses and the expected constant value (0) 
    else {
      differences <- other_gammas[initial_gamma_order[start:stop]]-norm_to_constant
      results[i,3] <- median(differences, na.rm = TRUE)
      results[i,5] <- mad(differences, na.rm = TRUE)
    }
    
    #report the number of non-NA mutant fitness values used
    results[i,4] <- sum(!is.na(differences))

    #apply the correction to the condition-specific changes in fitnesses and add it to the returned corrected vector, which preserves the original ordering of fitnesses.
    corrected_other[not_na[initial_gamma_order[start:stop]]] <- other_gammas[initial_gamma_order[start:stop]] - median(differences, na.rm = TRUE)
    
  }
  
  #at the conclusion of the loop, return the full table and vector of corrected condition-specific changes in fitnesses
  return(list(results, corrected_other, list_of_bins, list_of_bins2))
  
}
