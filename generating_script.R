# =====================================================================
# == This document generates the results in the analysis for Mi-Seon ==
# -- and Larry.                                                      ==
# == Author: Ethan Ancell                                            ==
# =====================================================================

# =====================================================================
# == NOTE: This file also implements the better version of wrapping. ==
# =====================================================================

library(tidyverse)
library(gridExtra)
library(R.utils)

# This file contains some helpful extra functions such as
#     (1) plotting functions
#     (2) functions to do wrapping of neurons
#     (3) old algorithm for assessing responsiveness of neuron
#     (4) new algorithm for assessing responsiveness of neuron
source('auxilary_functions.R')

# ===============
# == Load data ==
# ===============

d10 <- readRDS('data/appet_d10_matrix.Rds') # Appetitive conditioning day 10 data. 872 neurons, 20 trials, 720 datapoint per trial.
d11 <- readRDS('data/averse_d11_matrix.Rds') # Aversive conditioning. 519 neurons, 10 trials, 598 datapoint per trial.
d12 <- readRDS('data/tone_d12_matrix.Rds') # Post test. Alternating CS-food and CS-shock. 20 trials total, 10 trial each.  

# Reorder data dimensions to make consistent with 
# paper notation. With this new index ordering, the dimensions are
# (NEURONS, TRIAL, TIME POINT)
d10 <- wrap.array(d10, list(3,1,2))
d11 <- wrap.array(d11, list(3,1,2))
d12 <- wrap.array(d12, list(3,1,2))

# Split day 12 into the two types of tones
d12_f <- d12[, seq(1, 20, by = 2), ] # CS-food
d12_s <- d12[, seq(2, 20, by = 2), ] # CS-shock

# ======================================================
# == An example of what plotting data might look like ==
# ======================================================

# Plot 10th neuron in 5th trial in day 10 data
plot_one_trial(d10[10, 5, ], pre_window_start = 300, pre_window_end = 400)

# =============================
# == Run the algorithm above ==
# =============================

# Day 10 tone
 set.seed(1)
 res_d10_WSRT <- responsiveness_randomization(d10, pre_event = 100:300,
                                              post_event = 301:400,
                                              wrapping_method = "new",
                                              wrapping_period = 1:400, B = 500)
 # Day 10 food
 set.seed(2)
 res_d10_food_WSRT <- responsiveness_randomization(d10, pre_event = 250:430,
                                                   post_event = 431:600,
                                                   wrapping_method = "new",
                                                   wrapping_period = 1:dim(d10)[3],
                                                   B = 500)
 # Day 11 tone
 set.seed(3)
 res_d11_WSRT <- responsiveness_randomization(d11, pre_event = 100:300, post_event = 301:400,
                                              wrapping_method = "new",
                                              wrapping_period = 1:400, B = 500)
 # Day 11 shock
 set.seed(4)
 res_d11_shock_WSRT <- responsiveness_randomization(d11, pre_event = 250:430,
                                                    post_event = 431:480,
                                                    wrapping_method = "new",
                                                    wrapping_period = 1:dim(d11)[3], B = 500)
 # Day 12 food tone
 set.seed(5)
 res_d12_f_WSRT <- responsiveness_randomization(d12_f, pre_event = 100:300,
                                                post_event = 301:400,
                                                wrapping_method = "new",
                                                wrapping_period = 1:dim(d12_f)[3],
                                                B = 500)
 # Day 12 shock tone
 set.seed(6)
 res_d12_s_WSRT <- responsiveness_randomization(d12_s, pre_event = 100:300,
                                                post_event = 301:400,
                                                wrapping_method = "new",
                                                wrapping_period = 1:dim(d12_f)[3], B = 500)
# Save results so we don't have to run it again
 saveRDS(res_d10_WSRT, file = 'savedata/NEW_res_d10_WSRT.RDS')
 saveRDS(res_d10_food_WSRT, file = 'savedata/NEW_res_d10_food_WSRT.RDS')
 saveRDS(res_d11_WSRT, file = 'savedata/NEW_res_d11_WSRT.RDS')
 saveRDS(res_d11_shock_WSRT, file = 'savedata/NEW_res_d11_shock_WSRT.RDS')
 saveRDS(res_d12_f_WSRT, file = 'savedata/NEW_res_d12_f_WSRT.RDS')
 saveRDS(res_d12_s_WSRT, file = 'savedata/NEW_res_d12_s_WSRT.RDS')

# Load results
res_d10_WSRT <- readRDS(file = 'savedata/NEW_res_d10_WSRT.RDS')
res_d10_food_WSRT <- readRDS(file = 'savedata/NEW_res_d10_food_WSRT.RDS')
res_d11_WSRT <- readRDS(file = 'savedata/NEW_res_d11_WSRT.RDS')
res_d11_shock_WSRT <- readRDS(file = 'savedata/NEW_res_d11_shock_WSRT.RDS')
res_d12_f_WSRT <- readRDS(file = 'savedata/NEW_res_d12_f_WSRT.RDS')
res_d12_s_WSRT <- readRDS(file = 'savedata/NEW_res_d12_s_WSRT.RDS')

# Get p-values
res_d10_all <- alg1_pvals(res_d10_WSRT[['Wi_obs']], res_d10_WSRT[['Wi_rand']])
res_d10 <- res_d10_all[['p.vals']]

res_d10_food_all <- alg1_pvals(res_d10_food_WSRT[['Wi_obs']], res_d10_food_WSRT[['Wi_rand']])
res_d10_food <- res_d10_food_all[['p.vals']]

res_d11_all <- alg1_pvals(res_d11_WSRT[['Wi_obs']], res_d11_WSRT[['Wi_rand']])
res_d11 <- res_d11_all[['p.vals']]

res_d11_shock_all <- alg1_pvals(res_d11_shock_WSRT[['Wi_obs']], res_d11_shock_WSRT[['Wi_rand']])
res_d11_shock <- res_d11_shock_all[['p.vals']]

res_d12_f_all <- alg1_pvals(res_d12_f_WSRT[['Wi_obs']], res_d12_f_WSRT[['Wi_rand']])
res_d12_f <- res_d12_f_all[['p.vals']]

res_d12_s_all <- alg1_pvals(res_d12_s_WSRT[['Wi_obs']], res_d12_s_WSRT[['Wi_rand']])
res_d12_s <- res_d12_s_all[['p.vals']]

# ===============================
# == Find the significant ones ==
# ===============================

p_threshold <- 0.05

# Code positive response +1, inhibited is -1, non responsive is 0
coded_responses <- rep(0, dim(d10)[1])
coded_responses[res_d10_all$p.vals < p_threshold & res_d10_all$resp.type == '+'] <- 1
coded_responses[res_d10_all$p.vals < p_threshold & res_d10_all$resp.type == '-'] <- -1

# Get indices of positive responsive
indices_pos <- which(coded_responses == 1)
indices_neg <- which(coded_responses == -1)

# ====================================================================
# == Save a copy of the p-values as well as the response directions ==
# ====================================================================

# New wrapping
d10_df <- data.frame(neuron_index = 1:dim(d10)[1],
                     p_value_tone = res_d10,
                     direction_tone = res_d10_all[['resp.type']],
                     p_value_food = res_d10_food,
                     direction_food = res_d10_food_all[['resp.type']])
write.csv(d10_df, 'results/new_wrapping/results_d10.csv')

d11_df <- data.frame(neuron_index = 1:dim(d11)[1],
                     p_value_tone = res_d11,
                     direction_tone = res_d11_all[['resp.type']],
                     p_value_shock = res_d11_shock,
                     direction_shock = res_d11_shock_all[['resp.type']])
write.csv(d11_df, 'results/new_wrapping/results_d11.csv')

d12_df <- data.frame(neuron_index = 1:dim(d12)[1],
                     p_value_food_tone = res_d12_f,
                     direction_food_tone = res_d12_f_all[['resp.type']],
                     p_value_shock_tone = res_d12_s,
                     direction_shock_tone = res_d12_s_all[['resp.type']])
write.csv(d12_df, 'results/new_wrapping/results_d12.csv')


# ========================================
# == Evaluate responsiveness of neurons ==
# ========================================

# What proportion of neurons respond?
mean(res_d10 < 0.05) # Tone d10
mean(res_d10_food < 0.05) # Food d10
mean(res_d11 < 0.05) # Tone d11
mean(res_d11_shock < 0.05) # Shock d11
mean(res_d12_f < 0.05) # Food tone d12
mean(res_d12_s < 0.05) # Shock tone d12

# What proportion of these are positive and negative responses?
sum((res_d10_all[['resp.type']] == '+')[res_d10 < 0.05]) / length(res_d10)
sum((res_d10_all[['resp.type']] == '-')[res_d10 < 0.05]) / length(res_d10)

sum((res_d10_food_all[['resp.type']] == '+')[res_d10_food < 0.05]) / length(res_d10_food)
sum((res_d10_food_all[['resp.type']] == '-')[res_d10_food < 0.05]) / length(res_d10_food)

sum((res_d11_all[['resp.type']] == '+')[res_d11 < 0.05]) / length(res_d11)
sum((res_d11_all[['resp.type']] == '-')[res_d11 < 0.05]) / length(res_d11)

sum((res_d11_shock_all[['resp.type']] == '+')[res_d11_shock < 0.05]) / length(res_d11_shock)
sum((res_d11_shock_all[['resp.type']] == '-')[res_d11_shock < 0.05]) / length(res_d11_shock)

sum((res_d12_f_all[['resp.type']] == '+')[res_d12_f < 0.05]) / length(res_d12_f)
sum((res_d12_f_all[['resp.type']] == '-')[res_d12_f < 0.05]) / length(res_d12_f)

sum((res_d12_s_all[['resp.type']] == '+')[res_d12_s < 0.05]) / length(res_d12_s)
sum((res_d12_s_all[['resp.type']] == '-')[res_d12_s < 0.05]) / length(res_d12_s)
