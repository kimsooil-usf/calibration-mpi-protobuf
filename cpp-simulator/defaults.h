//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#ifndef DEFAULTS_H_
#define DEFAULTS_H_
#include <string>

struct defaults{
  std::string NUM_DAYS = "120";
  std::string START_DAY = "0";
  std::string INIT_FRAC_INFECTED = "0.0001";
  std::string INIT_FIXED_NUMBER_INFECTED = "100";
  std::string MEAN_INCUBATION_PERIOD = "4.50";
  std::string MEAN_ASYMPTOMATIC_PERIOD = "0.5";
  std::string MEAN_SYMPTOMATIC_PERIOD = "5";
  std::string SYMPTOMATIC_FRACTION = "0.67";
  std::string MEAN_HOSPITAL_REGULAR_PERIOD = "8";
  std::string MEAN_HOSPITAL_CRITICAL_PERIOD = "8";
  std::string MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD="2*365";//--Loss of immunity Shakir
  std::string COMPLIANCE_PROBABILITY = "1";
  std::string F_KERNEL_A = "10.751";
  std::string F_KERNEL_B = "5.384";
  std::string MEASURES ="1";//1,2,3,4,5 ---Check simulator.cc for explanation. It reduces or increases beta's by 25 percent.
                        //The first wo scenarios change all betas while the later two scenarios change 
                        //all beta's except for home.
  std::string REDUCTION_FACTOR="0.25";//Reduce or increase the betas by this factor to implement social measures
  std::string BETA_H = "0.15";//---->Comment by 1.227
  std::string BETA_W = "0.3";
  std::string BETA_C = "0.15";
  std::string BETA_S = "0.5";
  std::string BETA_TRAVEL = "0.1";
  std::string HD_AREA_FACTOR = "1.5";
  std::string HD_AREA_EXPONENT = "0";
  std::string INTERVENTION = "0";
  std::string CYCLIC_POLICY_TYPE = "0";
  std::string output_dir = "outputs/test_output_timing";
  std::string input_base = "../simulator/input_files";
  std::string attendance_filename = "attendance.json";
  std::string mask_filename = "maskwearing_November10_2022.json";//mask file
  std::string mask_folder = "maskensembles";//mask file
  std::string SEED_HD_AREA_POPULATION = "false";
  std::string SEED_ONLY_NON_COMMUTER = "false";
  std::string SEED_FIXED_NUMBER = "false";
  std::string IGNORE_ATTENDANCE_FILE = "true";
  std::string CALIBRATION_DELAY = "13";
  std::string DAYS_BEFORE_LOCKDOWN = "11";
  std::string FIRST_PERIOD = "30";
  std::string SECOND_PERIOD = "30";
  std::string THIRD_PERIOD = "90";
  std::string FOURTH_PERIOD = "270";
  std::string OE_SECOND_PERIOD = "30";
  std::string LOCKED_COMMUNITY_LEAKAGE = "1.0";
  std::string COMMUNITY_LOCK_THRESHOLD = "0.001";
  std::string LOCKED_NEIGHBORHOOD_LEAKAGE = "1.0";
  std::string NEIGHBORHOOD_LOCK_THRESHOLD = "0.001";
  std::string MASK_ACTIVE = "true";
  std::string MASK_FACTOR = "0.8";
  std::string MASK_START_DELAY = "0";
  std::string USE_AGE_DEPENDENT_MIXING = "false";
  std::string SIGNIFICANT_EIGEN_VALUES = "3.0";
  std::string NUM_AGE_GROUPS = "16";
  std::string CITY_SW_LAT = "27.649319212594804";//"12.8340125";
  std::string CITY_SW_LON = "-82.6280312634527";
  std::string CITY_NE_LAT = "28.171026529715622"; 
  std::string CITY_NE_LON = "-82.05614700746086";
  std::string NBR_CELL_SIZE = "0.178"; // in km
  std::string ENABLE_CONTAINMENT = "false";
  std::string ENABLE_NBR_CELLS = "true";
  std::string ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT = "false";
  std::string WARD_CONTAINMENT_THRESHOLD  = "0";
  std::string intervention_params_filename = "intervention_params.json";
  std::string BETA_PROJECT = ".2";
  std::string BETA_CLASS = "0.2";
  std::string BETA_RANDOM_COMMUNITY = ".2";
  std::string BETA_NBR_CELLS = "0.2";
  std::string ENABLE_TESTING = "false";
  std::string TESTING_PROTOCOL = "0";
  std::string testing_protocol_filename = "testing_protocol.json";
  // sk ///////////////////////////////////
  std::string BETA_COHORT = "0";
  std::string ENABLE_COHORTS = "false";
  std::string COHORT_SIZE = "15";
  std::string FRACTION_IN_TRAINS = "1.0";
  std::string CROWDING_FACTOR = "1.0";
  std::string ISOLATE_COHORTS = "false";
  std::string COHORT_SEVERITY_FRACTION = "0";
  std::string COHORT_STRATEGY = "1";
  std::string STORE_STATE_TIME_STEP = "0";
  std::string LOAD_STATE_TIME_STEP = "0"; // sk: allow load agents state day=0
  std::string ONE_OFF_TRAVELERS_RATIO = "0";
  std::string agent_load_file = "agentStore.pbstore";
  ////////////////////////////////
  std::string VACCINATION_EFFECTIVENESS1 = "0.75";
  std::string VACCINATION_EFFECTIVENESS2 = "0.95";
  std::string VACCINATION_EFFECTIVENESS_WANING = "0.9";
  std::string VACCINATION_EFFECTIVENESS_BOOSTED = "0.99";
  std::string VACCINATION_EFFECTIVENESS_WANING2 = "0.99";
  std::string VACCINATION_EFFECTIVENESS_BOOSTED2 = "0.99";

  // std::string TIME_ALPHA="1216";//30 december 2020
  // std::string TIME_DELTA="1708";//2 May 2021
  // std::string TIME_OMICRON="2564";//2 December 2021
  // std::string TIME_OMICRON_NEW="2772";//23 January 2022
  // std::string TIME_OMICRON_BA4="3192";//8 May 2022
  // std::string TIME_OMICRON_BA5="3192";//8 May 2022
  // std::string TIME_OMICRON_BA6="3472";//17 July 2022
  // std::string TIME_OMICRON_BA7="3640";//28 August 2022  std::string TIME_OMICRON_BA4="3204";//11 May 2022
  // std::string TIME_OMICRON_BA8="3660";//2 September 2022
  // // std::string TIME_OMICRON_BA9="3308";//6 June 2022
  // // std::string TIME_OMICRON_BA10="3308";//6 June 2022

  std::string TIME_ALPHA="1108";//3December 2020//"1000";//30 december 2020
  std::string TIME_DELTA="1680";//25 April 2021//"1828";//1 june 2021
  std::string TIME_OMICRON="2564";//2 December 2021
  std::string TIME_OMICRON_NEW="2772";//23 January 2022
  std::string TIME_OMICRON_BA4="3192";//8 May 2022
  std::string TIME_OMICRON_BA5="3192";//8 May 2022
  std::string TIME_OMICRON_BA6="3472";//17 July 2022
  std::string TIME_OMICRON_BA7="3640";//28 August 2022  std::string TIME_OMICRON_BA4="3204";//11 May 2022
  std::string TIME_OMICRON_BA8="3660";//2 September 2022
  // std::string TIME_OMICRON_BA9="3308";//6 June 2022
  // std::string TIME_OMICRON_BA10="3308";//6 June 2022


  std::string INFECTIOUSNESS_ALPHA = "1.2";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_DELTA = "1.35";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON = "5";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_NEW = "5.2";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA4 = "5.2";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA5 = "5.3";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA6 = "5.3";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA7 = "5.4";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA8 = "5.3";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA9 = "5.3";//---original 3 comment by Shakir
  std::string INFECTIOUSNESS_OMICRON_BA10 = "5.3";//---original 3 comment by Shakir

  // std::string INFECTIOUSNESS_ALPHA = "1.3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_DELTA = "2.6";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_NEW = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA4 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA5 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA6 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA7 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA8 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA9 = "3";//---original 3 comment by Shakir
  // std::string INFECTIOUSNESS_OMICRON_BA10 = "3";//---original 3 comment by Shakir

  std::string VIRULENT_NEW_STRAIN = "1.3";//Governs the probability of hospitalizations and subsequent criticality and death stages
  std::string VIRULENT_NEW_ALPHA="1.35";
  std::string VIRULENT_NEW_DELTA="3";
  std::string VIRULENT_NEW_OMICRON=".1";
  std::string VIRULENT_NEW_OMICRON_NEW=".1";
  std::string VIRULENT_NEW_OMICRON_BA4=".1";
  std::string VIRULENT_NEW_OMICRON_BA5=".1";
  std::string VIRULENT_NEW_OMICRON_BA6=".1";
  std::string VIRULENT_NEW_OMICRON_BA7=".1";
  std::string VIRULENT_NEW_OMICRON_BA8=".1";
  std::string VIRULENT_NEW_OMICRON_BA9=".1";
  std::string VIRULENT_NEW_OMICRON_BA10=".1";

  std::string RELATIVE_HOSPITALIZATION_ALPHA="1.3";
  std::string RELATIVE_HOSPITALIZATION_DELTA="3";
  std::string RELATIVE_HOSPITALIZATION_OMICRON=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_NEW=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_BA4=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_BA5=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_BA6=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_BA7=".1";
  std::string RELATIVE_HOSPITALIZATION_OMICRON_BA8=".1";

  // std::string REINFECTION_ALPHA="0.01";
  // std::string REINFECTION_DELTA="0.05";
  // std::string REINFECTION_OMICRON="0.1";  
  // std::string REINFECTION_OMICRON_NEW="0.2";
  // std::string REINFECTION_OMICRON_BA4="0.2";
  // std::string REINFECTION_OMICRON_BA5="0.2";
  // std::string REINFECTION_OMICRON_BA6="0.2";
  // std::string REINFECTION_OMICRON_BA7="0.2";
  // std::string REINFECTION_OMICRON_BA8="0.2";
  // std::string REINFECTION_OMICRON_BA9="0.2";
  // std::string REINFECTION_OMICRON_BA10="0.2";
  

  std::string FRACTION_NEW_STRAIN = "0.10";
  
  std::string FRACTION_NEW_ALPHA = "0.3";
  std::string FRACTION_NEW_DELTA = "1";
  std::string FRACTION_NEW_OMICRON = "0.40";
  std::string FRACTION_NEW_OMICRON_NEW = "0.30";
  std::string FRACTION_NEW_OMICRON_BA4 = "0.30";
  std::string FRACTION_NEW_OMICRON_BA5 = "0.30";
  std::string FRACTION_NEW_OMICRON_BA6 = "0.30";
  std::string FRACTION_NEW_OMICRON_BA7 = "0.30";  
  std::string FRACTION_NEW_OMICRON_BA8 = "0.30";
  std::string FRACTION_NEW_OMICRON_BA9 = "0.30";
  std::string FRACTION_NEW_OMICRON_BA10 = "0.30";

  // std::string FRACTION_NEW_STRAIN = "0.10";
  
  // std::string FRACTION_NEW_ALPHA = "0.01";
  // std::string FRACTION_NEW_DELTA = "0.10";
  // std::string FRACTION_NEW_OMICRON = "0.20";
  // std::string FRACTION_NEW_OMICRON_NEW = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA4 = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA5 = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA6 = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA7 = "0.30";  
  // std::string FRACTION_NEW_OMICRON_BA8 = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA9 = "0.30";
  // std::string FRACTION_NEW_OMICRON_BA10 = "0.30";


  // std::string FRACTION_SUSCEPTIBLE_ALPHA="1";
  // std::string FRACTION_SUSCEPTIBLE_DELTA="1";
  // std::string FRACTION_SUSCEPTIBLE_OMICRON="1";
  // std::string FRACTION_SUSCEPTIBLE_OMICRON_NEW="1";
  // std::string FRACTION_SUSCEPTIBLE_OMICRON_BA4="1";
  // std::string FRACTION_SUSCEPTIBLE_OMICRON_BA5="1";

} DEFAULTS;

#endif
