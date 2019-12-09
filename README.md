# Human confidence judgments reflect reliability-based hierarchical integration of contextual information

This repository contains the Matlab code to carry out the data analysis for the manuscript entitled "Human confidence judgments reflect reliability-based hierarchical integration of contextual information" by P. Schustek, A. Hyafil and R. Moreno-Bote, published in 2019 at Nature Communications ([link](https://www.nature.com/articles/s41467-019-13472-z)).


## Structure
```
./illustrations 
	/exp1 	
	/exp2
./models	# compute CVLL
	/exp1 	
	/exp2
./results
	/exp1 	
		/si
	/exp2
		/si
./src 		# function repository
```

## Usage
The following table indicates which Matlab script can be used to produce the results. Please note that model-based results require executing the respective models or some fits before.

| Figure | Script | 
|---------------|---------------|
|Fig. 1a | /illustrations/exp1/exp1_task_schema.m |
|Fig. 1b | /illustrations/exp1/exp1_sample_size_dependence.m |
|Fig. 1c | /illustrations/exp1/exp1_confH_mEv_N.m |
|Fig. 1d | /illustrations/exp1/exp1_slopes_psychometric.m |
|Fig. 2a | /results/exp1/basic_psychometric_N.m |
|Fig. 2b | /results/exp1/basic_N_regression.m |
|Fig. 3a-c | /illustrations/exp2/ |
|Fig. 4 | /results/exp2/pattern_overview.m |
|Fig. 5a | /results/exp2/prior_psychometric_context.m |
|Fig. 5b | /results/exp2/prior_opt_inferred_tendency.m |
|Fig. 6a | /results/exp2/prior_N_regression.m |
|Fig. 6b | /results/exp2/prior_modulation_N_current_trial.m |
|Fig. 6c | /results/exp2/prior_modulation_N_previous_trial.m |
|Fig. 7a | /results/exp2/prior_weight_previous_trials.m |
|Fig. 7b | /results/exp2/prior_accum_block_length.m |
|Fig. 8a | /results/exp2/prior_psychometric_context_fitted.m |
|Fig. 8b | /results/exp2/prior_opt_inferred_tendency_fitted.m |
|Fig. 8c | /results/exp2/prior_modulation_N_current_trial_fitted.m |
|Fig. 8d | /results/exp2/prior_modulation_N_previous_trial_fitted.m |
|Fig. 8e | /results/exp2/prior_weight_previous_trials_fitted.m |
|Fig. 8f | /results/exp2/prior_accum_block_length_fitted.m |
|Fig. S1a | /results/exp1/supp_info/si_basic_confidence_calibration.m |
|Fig. S1b | /results/exp2/supp_info/si_prior_confidence_calibration.m |
|Fig. S2 | /results/exp1/basic_N_regression.m |
|Fig. S3ab | /results/exp1/basic_N_regression.m |
|Fig. S3c | /results/exp1/supp_info/si_basic_N_regression_linearity.m |
|Fig. S4 | /results/exp1/supp_info/si_basic_bms_paired_ref_model_comparison.m |
|Fig. S5a | /results/exp1/supp_info/si_basic_sensory_noise.m |
|Fig. S5b | /results/exp2/supp_info/si_prior_sensory_noise.m |
|Fig. S6a | /results/exp1/supp_info/si_performance_trials_basic.m |
|Fig. S6b | /results/exp2/supp_info/si_performance_trials_prior.m |
|Fig. S7 | /results/exp2/supp_info/si_prior_bms_paired_ref_model_comparison.m |
|Fig. S8a | /results/exp2/prior_opt_inferred_tendency.m |
|Fig. S8b | /results/exp2/supp_info/si_prior_modulation_N_previous_message.m |
|Fig. S9 | /results/exp2/supp_info/si_prior_weight_previous_trials_spillover.m |
|Fig. S10a-e | /results/exp2/supp_info/si_prior_modulation_tally.m |
|Fig. S10f | /results/exp2/prior_modulation_N_previous_trial.m |


### Remarks on data availability
The original behavioral data used in this analysis is available as Supplementary Data to the manuscript (see link above). For convenience, the fields in data structure provide responses of the ideal observer model. 

