##The Best Alternative: The "Parameter Sweep" (AUC) Method
##If you do not want to defend the arbitrary choice of 5% in your paper, 
##you don't have to. The modern, highly robust solution used in advanced network neuroscience is 
##the Parameter Sweep (Area Under the Curve) method.

##Instead of picking one magic number, you calculate LAM across a range of reasonable thresholds, 
##and then average the results.

##How it works:You run your get_laminarity function at $RR = 0.01$ (1%).
##You run it again at $0.02$, $0.03$, $0.04$, and $0.05$.You plot these 5 LAM scores on a curve 
##and calculate the Area Under the Curve (AUC), or simply take the Mean LAM.By doing this, 
##you are telling the reviewers: "Subject 5 has a higher Laminarity than Subject 1 
##regardless of whether we define a recurrence as 1%, 3%, or 5% of the network. 
##Their temporal stability is fundamentally higher across all scales."


get_laminarity_sweep <- function(coords_mat, min_vert_line = 2) {
  # Define our sweep range: 1% to 5%
  rr_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  lam_scores <- c()
  
  for (target_RR in rr_thresholds) {
    # Call your original get_laminarity function here for each RR
    current_lam <- get_laminarity(coords_mat, target_RR, min_vert_line)
    lam_scores <- c(lam_scores, current_lam)
  }
  
  # Return the average Laminarity across the whole sweep
  mean_lam <- mean(lam_scores, na.rm = TRUE)
  return(mean_lam)
}