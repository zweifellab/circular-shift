# ======================
# == Useful functions ==
# ======================

# This function will plot a single trial (or average of trials). This function
# is useful when trying to visualize what is happening in one instance.
plot_one_trial <- function(slice, show_tone = TRUE, caption = NA,
                           pre_window_start = NA, pre_window_end = NA,
                           post_window_start = NA, post_window_end = NA,
                           chop_left = NA, chop_right = NA,
                           plot_pre_avg = FALSE, plot_post_avg = FALSE,
                           start_tone_line_index = 301,
                           end_tone_line_index = 400,
                           stimulus_index = 431) {
  data <- data.frame(time = 1:length(slice),
                     dFF = slice)
  lines <- data.frame(num = c(start_tone_line_index, end_tone_line_index, 
                              stimulus_index),
                      Time = c('start tone', 'end tone', 'stimulus'),
                      color = c('red', 'blue', 'purple'), 
                      stringsAsFactors = FALSE)
  
  cap <- ''
  if (!is.na(caption)) {
    cap <- caption
  }
  
  return_plot <- ggplot() +
    geom_line(data = data, mapping = aes(x = time, y = dFF)) +
    ggtitle(cap) +
    geom_vline(data = lines, mapping = aes(xintercept = num, color = Time))
  #geom_vline(data = lines, aes(xintercept = 301, color = 'red')) +
  #geom_vline(data = lines, aes(xintercept = 400, color = 'blue')) + 
  #geom_vline(aes(xintercept = 431, color = 'purple')) +
  #scale_color_discrete(labels = c('start tone', 'end tone', 'food'))
  
  # Add on a pretone window
  if (!is.na(pre_window_start) & !is.na(pre_window_end)) {
    return_plot <- return_plot +
      annotate("rect", xmin = pre_window_start, xmax = pre_window_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.2, fill = "green")
    
    # Plot the pretone average
    if (plot_pre_avg) {
      pretone_mean <- mean(data$dFF[pre_window_start:pre_window_end])
      return_plot <- return_plot +
        geom_segment(aes(x = pre_window_start, y = pretone_mean,
                         xend = pre_window_end, yend = pretone_mean))
    }
  }
  
  # Add on a posttone window
  if (!is.na(post_window_start) & !is.na(post_window_end)) {
    return_plot <- return_plot +
      annotate("rect", xmin = post_window_start, xmax = post_window_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.2, fill = "blue")
    
    # Plot the posttone average
    if (plot_post_avg) {
      posttone_mean <- mean(data$dFF[post_window_start:post_window_end])
      return_plot <- return_plot +
        geom_segment(aes(x = post_window_start, y = posttone_mean,
                         xend = post_window_end, yend = posttone_mean))
    }
  }
  
  # Crop
  if (!is.na(chop_left) & !is.na(chop_right)) {
    return_plot <- return_plot +
      coord_cartesian(xlim=c(chop_left, chop_right))
  }
  
  return_plot
}

# Wrap a neuron
wrap_data <- function(original_data, t) {
  if (t == 1) {
    return(original_data)
  } else {
    return(c(original_data[t:length(original_data)], original_data[1:(t-1)]))
  }
}

# ============================================
# == Algorithm for assessing responsiveness ==
# ============================================

responsiveness_randomization <- function(Y.in, pre_event, post_event, 
                                         wrapping_method = "new",
                                         wrapping_period = NULL, B = 100) {
  
  # Make sure wrapping method is valid
  if (!(wrapping_method %in% c("new", "old"))) {
    stop("Wrapping method must be \"old\" or \"new\".")
  }
  
  # By default, wrap on all the data
  if (missing(wrapping_period)) {
    wrapping_period <- 1:dim(Y.in)[3]
  }
  
  # Only care about data lying within the wrapping period (after subsetting,
  # we don't have to think about this anymore)
  Y <- Y.in[, , wrapping_period]
  
  n.c <- dim(Y)[1]
  n.tr <- dim(Y)[2]
  timepoints <- dim(Y)[3]
  
  # Stored results
  WSRT_obs <- rep(NA, n.c)
  WSRT_randomization <- array(numeric(), c(n.c, B))
  
  # NEW METHOD
  # Use all neurons to get null distributions
  if (wrapping_method == "new") {
    
    # Create a repository of test statistics for null data
    Wib.tilde <- rep(0, B)
    for (b in 1:B) {
      i.star <- sample(1:n.c, size = 1)
      for (j in 1:n.tr) {
        wrap.q <- sample(1:timepoints, size = 1)
        Yij.tilde <- wrap_data(Y[i.star, j, ], wrap.q)
        Yij.tilde.pre <- Yij.tilde[pre_event]
        Yij.tilde.post <- Yij.tilde[post_event]
        Wib.tilde[b] <- Wib.tilde[b] + wilcox.test(Yij.tilde.post, 
                                                   Yij.tilde.pre, paired = FALSE)$statistic
      }
    }
    # Store it all
    for (i in 1:n.c) {
      WSRT_randomization[i, ] <- Wib.tilde
    }
  }
  # OLD METHOD
  # Wrap each neuron individually to get null distributions
  else {
    # Loop through neurons
    for (i in 1:n.c) {
      
      Wib.tilde <- rep(0, B)
      # Loop through randomizations
      for (b in 1:B) {
        wrap.t <- sample(1:timepoints, size = 1)
        # Loop through trials
        for (j in 1:n.tr) {
          Yij.tilde <- wrap_data(Y[i, j, ], wrap.t)
          Yij.tilde.pre <- Yij.tilde[pre_event]
          Yij.tilde.post <- Yij.tilde[post_event]
          Wib.tilde[b] <- Wib.tilde[b] + wilcox.test(Yij.tilde.post, 
                                                     Yij.tilde.pre, paired = FALSE)$statistic
        }
      }
      
      # Store randomization WSRT
      WSRT_randomization[i, ] <- Wib.tilde
    }
  }
  
  # Observed test statistics method stays the same
  for (i in 1:n.c) {
    Wi.tilde <- 0
    for (j in 1:n.tr) {
      Yij.pre <- Y[i, j, pre_event]
      Yij.post <- Y[i, j, post_event]
      Wi.tilde <- Wi.tilde + wilcox.test(Yij.post, Yij.pre, paired = FALSE)$statistic
    }
    
    # Store observed WSRT
    WSRT_obs[i] <- Wi.tilde
  }
  
  return(list(Wi_obs = WSRT_obs, Wi_rand = WSRT_randomization))
}

# Function to get p-values based upon the WSRT statistics returned from the
# algorithm above.
alg1_pvals <- function(Wi_obs, Wi_rand) {
  n.c <- length(Wi_obs)
  p.vals <- rep(NA, n.c)
  resp.type <- rep('X', n.c)
  for (i in 1:n.c) {
    # p.val.1 <- mean(Wi_rand[i, ] >= Wi_obs[i])
    # p.val.2 <- mean(Wi_rand[i, ] <= Wi_obs[i])
    p.val.1 <- (sum(Wi_rand[i, ] >= Wi_obs[i]) + 1) / (length(Wi_rand[i, ]) + 1)
    p.val.2 <- (sum(Wi_rand[i, ] <= Wi_obs[i]) + 1) / (length(Wi_rand[i, ]) + 1)
    
    p.vals[i] <- 2*min(p.val.1, p.val.2)
    
    if (p.val.1 < p.val.2) {
      resp.type[i] <- '+'
    } else {
      resp.type[i] <- '-'
    }
  }
  
  return(list(p.vals = p.vals, resp.type = resp.type))
}