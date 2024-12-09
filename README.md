# drydown-of-NNsm
Prepare data on an annual basis:
1) Soil moisture: NNsm_year_YEAR.MAT
2) precipitation: GPM_EASE_Daily_6AM_Local_36km_YEAR01_YEAR12.mat
  YEAR:2002-2023

Code processing is divided into three steps:

Step 1: Run NNsm_tau_process_xyw_f1_20022006. m

step 2: Run NNsm_post_multidata_tauL_process_f1.m

Step 3: Run NNsm_stage1memory_time_f1.m

step 4: Run NNSM_soil_thresholds_process.m

  frequency: f1 or f3
  
  year: 20022006, 20072011, 20122016, 20172021, 20222023
