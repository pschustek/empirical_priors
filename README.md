# Human confidence judgments reflect reliability-based hierarchical integration of contextual information

This repository contains the Matlab code to carry out the data analysis for the manuscript entitled "Human confidence judgments reflect reliability-based hierarchical integration of contextual information" currently available as a [preprint ](https://doi.org/10.1101/425462).


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
The following table indicates which Matlab script can be used to produce the results.

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
|Fig. S1a | /results/exp1/supp_info/si_basic_confidence_calibration.m |
|Fig. S1b | /results/exp1/supp_info/si_prior_confidence_calibration.m |
|Fig. S2 | /results/exp1/supp_info/si_basic_bms_paired_ref_model_comparison.m |
|Fig. S3 | /results/exp2/supp_info/si_prior_bms_paired_ref_model_comparison.m |
|Fig. S6a | /results/exp2/supp_info/si_prior_psychometric_context.m |
|Fig. S6b | /results/exp2/supp_info/si_prior_opt_inferred_tendency.m |
|Fig. S6c | /results/exp2/supp_info/si_prior_modulation_N_current_trial.m |
|Fig. S6d | /results/exp2/supp_info/si_prior_modulation_N_previous_trial.m |
|Fig. S6e | /results/exp2/supp_info/si_prior_weight_previous_trials.m |
|Fig. S6f | /results/exp2/supp_info/si_prior_accum_block_length.m |
Some more code refactoring might be done in the future, depending as well on changes to the analysis.

### Remarks on data availability
The original behavioral data used in this analysis will be released upon final publication of the manuscript. For convenience, the fields in data structure provide responses of the ideal observer model. 

