# EKG_EEG

This repository contains some scripts to analyze simultaneous EEG and ECG recordings.

To detrend the Plethysmography it uses the code available here: https://github.com/cvxgrp/l1_tf/tree/master/matlab
Functions for calculating mutual information and other information theoretic quantities using a parametric Gaussian copula: https://github.com/robince/gcmi/tree/master/matlab

Overview of the scripts:

   extract_HR_and_epochs_DEAP.m              detect R peaks and extract EEG epochs (1 trial)
   loop_check_HR_DEAP.m                      visualize the peaks
   loop_extract_HR_and_epochs_DEAP.m         detect R peaks and extract EEG epochs (loop over trials and subjects)
   calc_temp_inf_HEP_DEAP.m                  calculate the information theoretic quantities + create plots
   calc_temp_inf_HEP_DEAP_test_cluster.m     test for significance using cluster-based permutation test (multiple comparisons corrected)
   calc_temp_inf_HEP_DEAP_test_pixel.m       test for significance using permutation test at pixel level (uncorrected)
