//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#include "simulator.h"
#include "outputs.h"
#include "defaults.h"
#include "initializers.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cxxopts.hpp>


int main(int argc, char** argv){
  cxxopts::Options options(argv[0],
			   "Simulate the mean field agent model");

  options.add_options("Basic")
    ("h,help", "display description of program options")
    ("NUM_DAYS", "number of days in the simulation",
     cxxopts::value<count_type>()->default_value(DEFAULTS.NUM_DAYS))
    ("START_DAY", "day of starting the simulator after March 1 2020",
     cxxopts::value<count_type>()->default_value(DEFAULTS.START_DAY))     
    ("output_directory", "output directory",
     cxxopts::value<std::string>()->default_value(DEFAULTS.output_dir))
    ("input_directory", "input directory",
     cxxopts::value<std::string>()->default_value(DEFAULTS.input_base))
    ("PROVIDE_INITIAL_SEED",
     "provide an initial seed to the simulator. If this is not provided, the simulator uses "
     "std::random_device to get the random seed.",
     cxxopts::value<count_type>())
    ("PROVIDE_INITIAL_SEED_GRAPH",
     "provide the initial seed for the interaction graphs. If this is not provided, the simulator uses "
     "std::random_device to get the random seed.",
     cxxopts::value<count_type>())
    ;

  options.add_options("Infection seeding")
    ("SEED_HD_AREA_POPULATION", "seed those living in high-density areas as well",
     cxxopts::value<bool>()->default_value(DEFAULTS.SEED_HD_AREA_POPULATION))
    ("SEED_ONLY_NON_COMMUTER", "seed only those who do not take public transit",
     cxxopts::value<bool>()->default_value(DEFAULTS.SEED_ONLY_NON_COMMUTER))
    ("INIT_FRAC_INFECTED", "initial probability of a person being infected.  If --SEED_FIXED_NUMBER is provided, this is ignored in favour of INIT_FIXED_NUMBER_INFECTED",
     cxxopts::value<double>()->default_value(DEFAULTS.INIT_FRAC_INFECTED))
    ("SEED_FIXED_NUMBER", "seed a fixed number of initial infections.  If this option is provided, INIT_FRAC_INFECTED is ignored in favour of INIT_FIXED_NUMBER_INFECTED",
     cxxopts::value<bool>()->default_value(DEFAULTS.SEED_FIXED_NUMBER))
    ("INIT_FIXED_NUMBER_INFECTED", "initial number of people infected.  If --SEED_FIXED_NUMBER is provided, this supersedes INIT_FRAC_INFECTED",
     cxxopts::value<count_type>()->default_value(DEFAULTS.INIT_FIXED_NUMBER_INFECTED))
    ;

  options.add_options("Disease progression")
    ("MEAN_INCUBATION_PERIOD", "mean incubation period",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_INCUBATION_PERIOD))
    ("MEAN_ASYMPTOMATIC_PERIOD", "mean asymptomatic period",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_ASYMPTOMATIC_PERIOD))
    ("MEAN_SYMPTOMATIC_PERIOD", "mean symptomatic period",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_SYMPTOMATIC_PERIOD))
    ("SYMPTOMATIC_FRACTION", "fraction of people who develop symptoms",
     cxxopts::value<double>()->default_value(DEFAULTS.SYMPTOMATIC_FRACTION))
    ("MEAN_HOSPITAL_REGULAR_PERIOD", "mean period of regular hospitalization",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_HOSPITAL_REGULAR_PERIOD))
    ("MEAN_HOSPITAL_CRITICAL_PERIOD", "mean period of critical care hospitalization",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_HOSPITAL_CRITICAL_PERIOD))
    ("MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD", "Loss of immunity by the recovered",
     cxxopts::value<double>()->default_value(DEFAULTS.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD))     
    ;

  options.add_options("City")
    ("COMPLIANCE_PROBABILITY", "default compliance probability",
     cxxopts::value<double>()->default_value(DEFAULTS.COMPLIANCE_PROBABILITY))
    ("HD_COMPLIANCE_PROBABILITY", "default compliance probability for high-density areas",
     cxxopts::value<double>()->default_value(DEFAULTS.COMPLIANCE_PROBABILITY))      
    ("F_KERNEL_A", "the 'a' parameter in the distance kernel, for distance in km",
     cxxopts::value<double>()->default_value(DEFAULTS.F_KERNEL_A))
    ("F_KERNEL_B", "the 'b' parameter in the distance kernel, for distance in km",
     cxxopts::value<double>()->default_value(DEFAULTS.F_KERNEL_B))
    ("BETA_H", "the beta_home parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_H))
    ("BETA_W", "the beta_workplace parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_W))
    ("BETA_C", "the beta_community parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_C))
    ("BETA_S", "the beta_school parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_S))
    ("BETA_TRAVEL", "the beta_travel parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_TRAVEL))
    ("BETA_CLASS", "the beta_class parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_CLASS))
    ("BETA_PROJECT", "the beta_project parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_PROJECT))
    ("BETA_RANDOM_COMMUNITY", "the beta_random_community parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_RANDOM_COMMUNITY))
    ("BETA_NBR_CELLS", "the beta_travel parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_NBR_CELLS))
    ("BETA_COHORT", "the beta_cohorts parameter",
     cxxopts::value<double>()->default_value(DEFAULTS.BETA_COHORT))
    ("HD_AREA_FACTOR", "multiplicative factor for high density areas",
     cxxopts::value<double>()->default_value(DEFAULTS.HD_AREA_FACTOR))
    ("HD_AREA_EXPONENT", "exponent for community size for high density areas",
     cxxopts::value<double>()->default_value(DEFAULTS.HD_AREA_EXPONENT))
    ;

  options.add_options("Intervention - basic")
    ("INTERVENTION", "index of the intervention",
     cxxopts::value<count_type>()->default_value(DEFAULTS.INTERVENTION))
    ("CALIBRATION_DELAY", "delay observed in calibration",
     cxxopts::value<double>()->default_value(DEFAULTS.CALIBRATION_DELAY))
    ("MEASURES", "Increase or decrease social measure for school office or work/community +1 or -1 to increase or decrease social measures",
     cxxopts::value<int>()->default_value(DEFAULTS.MEASURES))
    ("DAYS_BEFORE_LOCKDOWN", "no intervention period prior to interventions",
     cxxopts::value<double>()->default_value(DEFAULTS.DAYS_BEFORE_LOCKDOWN))
    ("FIRST_PERIOD", "length in days of the first intervention period",
     cxxopts::value<double>()->default_value(DEFAULTS.FIRST_PERIOD))
    ("SECOND_PERIOD", "length in days of the second intervention period",
     cxxopts::value<double>()->default_value(DEFAULTS.SECOND_PERIOD))
    ("THIRD_PERIOD", "length in days of the third intervention period",
     cxxopts::value<double>()->default_value(DEFAULTS.THIRD_PERIOD))
    ("FOURTH_PERIOD", "length in days of the fourth intervention period",
     cxxopts::value<double>()->default_value(DEFAULTS.FOURTH_PERIOD))     
    ("OE_SECOND_PERIOD", "length in days of the second odd-even intervention period",
     cxxopts::value<double>()->default_value(DEFAULTS.OE_SECOND_PERIOD))
    ;

  options.add_options("Intervention - cylic strategy")
    ("CYCLIC_POLICY_TYPE", "whether cyclic policy (only relevant in interventions which enable cyclic work schedule) is implemented based on homes or individuals",
     cxxopts::value<count_type>()->default_value(DEFAULTS.CYCLIC_POLICY_TYPE))
    ;

  options.add_options("Intervention - soft containment zones")
    ("LOCKED_COMMUNITY_LEAKAGE", "minimum community infection leakage under containment",
     cxxopts::value<double>()->default_value(DEFAULTS.LOCKED_COMMUNITY_LEAKAGE))
    ("COMMUNITY_LOCK_THRESHOLD", "hospitalisation fraction in a ward beyond which the ward will be cordoned off.",
     cxxopts::value<double>()->default_value(DEFAULTS.COMMUNITY_LOCK_THRESHOLD))
	("LOCKED_NEIGHBORHOOD_LEAKAGE", "minimum neighborhood cell infection leakage under containment",
	 cxxopts::value<double>()->default_value(DEFAULTS.LOCKED_NEIGHBORHOOD_LEAKAGE))
    ("NEIGHBORHOOD_LOCK_THRESHOLD", "hospitalisation fraction in a neighbourhood cell beyond which the neighbourhood cell will be cordoned off.",
     cxxopts::value<double>()->default_value(DEFAULTS.NEIGHBORHOOD_LOCK_THRESHOLD))
	("ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT", "whether neighborhood soft containment is enabled (false if neighborhood cells are not enabled)",
     cxxopts::value<bool>()->default_value(DEFAULTS.ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT))
	;

  options.add_options("Intervention - neighbourhood containment")
    ("ENABLE_CONTAINMENT", "enable neighborhood containment",
     cxxopts::value<bool>()->default_value(DEFAULTS.ENABLE_CONTAINMENT))
    ("ENABLE_NBR_CELLS", "Enable neighbourhood cells",
     cxxopts::value<bool>()->default_value(DEFAULTS.ENABLE_NBR_CELLS))    
    ("CITY_SW_LAT", "south west latitude boundary of the city",
     cxxopts::value<double>()->default_value(DEFAULTS.CITY_SW_LAT))
    ("CITY_SW_LON", "south west longitude boundary of the city",
     cxxopts::value<double>()->default_value(DEFAULTS.CITY_SW_LON))
    ("CITY_NE_LAT", "north east latitude boundary of the city",
     cxxopts::value<double>()->default_value(DEFAULTS.CITY_NE_LAT))
    ("CITY_NE_LON", "north east longitude boundary of the city",
     cxxopts::value<double>()->default_value(DEFAULTS.CITY_NE_LON))
    ("NBR_CELL_SIZE", "neighbourhood cell size (length of side)",
     cxxopts::value<double>()->default_value(DEFAULTS.NBR_CELL_SIZE))
    ("WARD_CONTAINMENT_THRESHOLD",
	 "threhsold on the number of hospitalized individuals in a ward beyond which the ward becomes a containment zone.",
     cxxopts::value<count_type>()->default_value(DEFAULTS.WARD_CONTAINMENT_THRESHOLD))
    ;

  options.add_options("Age-dependent mixing")
    ("USE_AGE_DEPENDENT_MIXING", "whether age-stratified interactions are enabled",
     cxxopts::value<bool>()->default_value(DEFAULTS.USE_AGE_DEPENDENT_MIXING))
    ("SIGNIFICANT_EIGEN_VALUES",
	 "number of principal components of the age-stratification matrix to use",
     cxxopts::value<double>()->default_value(DEFAULTS.SIGNIFICANT_EIGEN_VALUES))
    ("NUM_AGE_GROUPS", "number of age group buckets",
     cxxopts::value<count_type>()->default_value(DEFAULTS.NUM_AGE_GROUPS))
    ;

  options.add_options("Other")
    ("IGNORE_ATTENDANCE_FILE", "whether to ignore the attendance file",
     cxxopts::value<bool>()->default_value(DEFAULTS.IGNORE_ATTENDANCE_FILE))
    ("attendance_filename", "attendance json filename, relative to input_directory",
     cxxopts::value<std::string>()->default_value(DEFAULTS.attendance_filename))
    ("intervention_filename", "intervention json filename, relative to input_directory",
     cxxopts::value<std::string>()->default_value(DEFAULTS.intervention_params_filename))
    ("MASK_ACTIVE", "whether masks are to be incorporated",
     cxxopts::value<bool>()->default_value(DEFAULTS.MASK_ACTIVE))
    ("MASK_FACTOR", "whether masks are to be incorporated",
     cxxopts::value<double>()->default_value(DEFAULTS.MASK_FACTOR))
    ("MASK_START_DELAY", "days after which masks are enforced",
     cxxopts::value<double>()->default_value(DEFAULTS.MASK_START_DELAY))
    ;

    options.add_options("Testing and contact tracing")
    ("ENABLE_TESTING", "enable testing (contact-tracing) functionality",
     cxxopts::value<bool>()->default_value(DEFAULTS.ENABLE_TESTING))
    ("testing_protocol_filename", "intervention json filename, relative to input_directory",
    cxxopts::value<std::string>()->default_value(DEFAULTS.testing_protocol_filename))
    ("TESTING_PROTOCOL", "index of testing protocol",
     cxxopts::value<count_type>()->default_value(DEFAULTS.TESTING_PROTOCOL))
    ;
	
	options.add_options("New_strain_vaccination_reinfection")
    ("VACCINATION_EFFECTIVENESS1", "vaccination effectiveness 1st dose",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS1))
    ("VACCINATION_EFFECTIVENESS2", "vaccination effectiveness 2nd dose",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS2))
    ("VACCINATION_EFFECTIVENESS_BOOSTED", "vaccination effectiveness boosted",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS_BOOSTED))
    ("VACCINATION_EFFECTIVENESS_WANING", "vaccination effectiveness waning",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS_WANING))
    ("VACCINATION_EFFECTIVENESS_WANING2", "vaccination effectiveness waning2",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS_WANING2))
    ("VACCINATION_EFFECTIVENESS_BOOSTED2", "vaccination effectiveness boosted2",
     cxxopts::value<double>()->default_value(DEFAULTS.VACCINATION_EFFECTIVENESS_BOOSTED2))     
    ("TIME_ALPHA", "Time at which alpha",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_ALPHA))     
    ("TIME_DELTA", "Time at which delta",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_DELTA))
    ("TIME_OMICRON", "Time at which omicron",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_OMICRON))             
    ("TIME_OMICRON_NEW", "Time at which omicron_new",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_OMICRON_NEW))             
    ("TIME_OMICRON_BA4", "Time at which omicron_new",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_OMICRON_BA4))           
    ("TIME_OMICRON_BA5", "Time at which omicron_new",
    cxxopts::value<double>()->default_value(DEFAULTS.TIME_OMICRON_BA5))           

    ("INFECTIOUSNESS_ALPHA", "new alpha strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_ALPHA))
    ("INFECTIOUSNESS_DELTA", "new delta strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_DELTA))    
    ("INFECTIOUSNESS_OMICRON", "new omicron strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_OMICRON))    
    ("INFECTIOUSNESS_OMICRON_NEW", "new omicronNew strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_OMICRON_NEW))                 
    ("INFECTIOUSNESS_OMICRON_BA4", "new omicronNew strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_OMICRON_BA4))                 

  ("INFECTIOUSNESS_OMICRON_BA5", "new omicronNew strain increases infectiousness by this factor",
    cxxopts::value<double>()->default_value(DEFAULTS.INFECTIOUSNESS_OMICRON_BA5))    

    ("VIRULENT_NEW_STRAIN", "new strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_STRAIN))
    ("VIRULENT_NEW_ALPHA", "alpha strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_ALPHA))
    ("VIRULENT_NEW_DELTA", "delta strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_DELTA))
    ("VIRULENT_NEW_OMICRON", "omicron strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_OMICRON))               
    ("VIRULENT_NEW_OMICRON_NEW", "omicron_new strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_OMICRON_NEW))                    
    ("VIRULENT_NEW_OMICRON_BA4", "omicron_new strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_OMICRON_BA4))  
    ("VIRULENT_NEW_OMICRON_BA5", "omicron_new strain is more deadly by this factor",
     cxxopts::value<double>()->default_value(DEFAULTS.VIRULENT_NEW_OMICRON_BA5)) 

    ("REINFECTION_ALPHA", "fraction of recovered people eligible for reinfection by alpha on time_alpha",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_ALPHA))
    ("REINFECTION_DELTA", "fraction of recovered people eligible for reinfection by delta on time_delta",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_DELTA))
    ("REINFECTION_OMICRON", "fraction of recovered people eligible for reinfection by omicron on time_omicron",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_OMICRON))
    ("REINFECTION_OMICRON_NEW", "fraction of recovered people eligible for reinfection by omicron_new on time_omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_OMICRON_NEW))

    ("REINFECTION_OMICRON_BA4", "fraction of recovered people eligible for reinfection by omicron_new on time_omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_OMICRON_BA4))
    ("REINFECTION_OMICRON_BA5", "fraction of recovered people eligible for reinfection by omicron_new on time_omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.REINFECTION_OMICRON_BA5))

    ("FRACTION_NEW_STRAIN", "fraction of infected people infected with alpha",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_STRAIN))
    ("FRACTION_NEW_ALPHA", "fraction of infected people infected with delta",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_DELTA))
    ("FRACTION_NEW_DELTA", "fraction of infected people infected with omicron",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_STRAIN))
    ("FRACTION_NEW_OMICRON", "fraction of infected people infected with omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_OMICRON))
    ("FRACTION_NEW_OMICRON_NEW", "fraction of infected people infected with omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_OMICRON_NEW))     

    ("FRACTION_NEW_OMICRON_BA4", "fraction of infected people infected with omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_OMICRON_BA4))   
    ("FRACTION_NEW_OMICRON_BA5", "fraction of infected people infected with omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_NEW_OMICRON_BA5)) 

    ("FRACTION_SUSCEPTIBLE_ALPHA", "fraction of susceptible people to alpha",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_ALPHA))     
    ("FRACTION_SUSCEPTIBLE_DELTA", "fraction of susceptible people to delta",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_DELTA))          
    ("FRACTION_SUSCEPTIBLE_OMICRON", "fraction of susceptible people to omicron",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_OMICRON))          
    ("FRACTION_SUSCEPTIBLE_OMICRON_NEW", "fraction of susceptible people to omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_OMICRON_NEW))               
    ("FRACTION_SUSCEPTIBLE_OMICRON_BA4", "fraction of susceptible people to omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_OMICRON_BA4))    
    ("FRACTION_SUSCEPTIBLE_OMICRON_BA5", "fraction of susceptible people to omicron_new",
     cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_SUSCEPTIBLE_OMICRON_BA5))         
;
  options.add_options("Cohorts")
    ("ENABLE_COHORTS", "enable cohorts",
     cxxopts::value<bool>()->default_value(DEFAULTS.ENABLE_COHORTS))
    ("ISOLATE_COHORTS", "isolate cohorts",
     cxxopts::value<bool>()->default_value(DEFAULTS.ISOLATE_COHORTS))
    ("CROWDING_FACTOR_COHORTS", "crowding factor in trains in cohorts",
    cxxopts::value<double>()->default_value(DEFAULTS.CROWDING_FACTOR))
    ("COHORT_SIZE", "size of a cohort",
    cxxopts::value<double>()->default_value(DEFAULTS.COHORT_SIZE))
    ("FRACTION_IN_TRAINS_COHORTS", "fraction in trains",
    cxxopts::value<double>()->default_value(DEFAULTS.FRACTION_IN_TRAINS))
    ("COHORT_SEVERITY_FRACTION","fraction of symptomatic who will trigger isolation",
    cxxopts::value<double>()->default_value(DEFAULTS.COHORT_SEVERITY_FRACTION))
    ("COHORT_STRATEGY", "index of the cohort strategy",
     cxxopts::value<count_type>()->default_value(DEFAULTS.COHORT_STRATEGY))
    ("ONE_OFF_TRAVELERS_RATIO", "ratio of one-off travelers",
     cxxopts::value<double>()->default_value(DEFAULTS.ONE_OFF_TRAVELERS_RATIO))
    ;

  options.add_options("Store_load_state")
  ("STORE_STATE_TIME_STEP", "state storage timestep",
   cxxopts::value<count_type>()->default_value(DEFAULTS.STORE_STATE_TIME_STEP))
  ("LOAD_STATE_TIME_STEP", "state load timestep",
   cxxopts::value<count_type>()->default_value(DEFAULTS.LOAD_STATE_TIME_STEP))
  ("agent_load_file", "location of agentStore.pbstore, relative to input_directory",
    cxxopts::value<std::string>()->default_value(DEFAULTS.agent_load_file))
  ;

  auto optvals = options.parse(argc, argv);
  
  if(optvals.count("help")){
    std::cout << options.help({"Basic",
			       "Infection seeding",
			       "Disease progression",
			       "City",
			       "Intervention - basic",
			       "Intervention - cyclic strategy",
			       "Intervention - soft containment zones",
			       "Intervention - neighbourhood containment",
			       "Age-dependent mixing",
			       "Other",
             "Testing and contact tracing",
             "New_strain_vaccination_reinfection"
      }) << std::endl;
    return 0;
  }
  
  //Save options
  GLOBAL.SEED_HD_AREA_POPULATION = optvals["SEED_HD_AREA_POPULATION"].count();
  GLOBAL.SEED_ONLY_NON_COMMUTER = optvals["SEED_ONLY_NON_COMMUTER"].count();
  GLOBAL.SEED_FIXED_NUMBER = optvals["SEED_FIXED_NUMBER"].count();
  GLOBAL.NUM_DAYS = optvals["NUM_DAYS"].as<count_type>();
  GLOBAL.START_DAY = optvals["START_DAY"].as<count_type>();
  
  GLOBAL.INIT_FRAC_INFECTED = optvals["INIT_FRAC_INFECTED"].as<double>();
  GLOBAL.INIT_FIXED_NUMBER_INFECTED = optvals["INIT_FIXED_NUMBER_INFECTED"].as<count_type>();
  GLOBAL.MEAN_INCUBATION_PERIOD = optvals["MEAN_INCUBATION_PERIOD"].as<double>();
  GLOBAL.MEAN_ASYMPTOMATIC_PERIOD = optvals["MEAN_ASYMPTOMATIC_PERIOD"].as<double>();
  GLOBAL.MEAN_SYMPTOMATIC_PERIOD = optvals["MEAN_SYMPTOMATIC_PERIOD"].as<double>();
  GLOBAL.SYMPTOMATIC_FRACTION = optvals["SYMPTOMATIC_FRACTION"].as<double>();
  GLOBAL.MEAN_HOSPITAL_REGULAR_PERIOD = optvals["MEAN_HOSPITAL_REGULAR_PERIOD"].as<double>();
  GLOBAL.MEAN_HOSPITAL_CRITICAL_PERIOD = optvals["MEAN_HOSPITAL_CRITICAL_PERIOD"].as<double>();

  GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD = optvals["MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD"].as<double>();

  GLOBAL.COMPLIANCE_PROBABILITY = optvals["COMPLIANCE_PROBABILITY"].as<double>();
  if(optvals["HD_COMPLIANCE_PROBABILITY"].count()){
	  GLOBAL.HD_COMPLIANCE_PROBABILITY = optvals["HD_COMPLIANCE_PROBABILITY"].as<double>();
  } else {
	//If HD_COMPLIANCE_PROBABILITY is not provided then set it equal to whatever
	//value was provided for COMPLIANCE_PROBABILITY
	GLOBAL.HD_COMPLIANCE_PROBABILITY = GLOBAL.COMPLIANCE_PROBABILITY;
  }
  GLOBAL.F_KERNEL_A = optvals["F_KERNEL_A"].as<double>();
  GLOBAL.F_KERNEL_B = optvals["F_KERNEL_B"].as<double>();

  GLOBAL.BETA_H = optvals["BETA_H"].as<double>();
  GLOBAL.BETA_W = optvals["BETA_W"].as<double>();
  GLOBAL.BETA_C = optvals["BETA_C"].as<double>();
  GLOBAL.BETA_S = optvals["BETA_S"].as<double>();
  GLOBAL.BETA_TRAVEL = optvals["BETA_TRAVEL"].as<double>();
  GLOBAL.BETA_CLASS = optvals["BETA_CLASS"].as<double>();
  GLOBAL.BETA_PROJECT = optvals["BETA_PROJECT"].as<double>();
  GLOBAL.BETA_RANDOM_COMMUNITY = optvals["BETA_RANDOM_COMMUNITY"].as<double>();
  GLOBAL.BETA_NBR_CELLS = optvals["BETA_NBR_CELLS"].as<double>();

  GLOBAL.HD_AREA_FACTOR = optvals["HD_AREA_FACTOR"].as<double>();
  GLOBAL.HD_AREA_EXPONENT = optvals["HD_AREA_EXPONENT"].as<double>();
  
  GLOBAL.INTERVENTION
	= static_cast<Intervention>(optvals["INTERVENTION"].as<count_type>());
  GLOBAL.intervention_filename = optvals["intervention_filename"].as<std::string>();

  GLOBAL.MEASURES = optvals["MEASURES"].as<int>();
  GLOBAL.CALIBRATION_DELAY = optvals["CALIBRATION_DELAY"].as<double>();
  GLOBAL.DAYS_BEFORE_LOCKDOWN = optvals["DAYS_BEFORE_LOCKDOWN"].as<double>();
  GLOBAL.NUM_DAYS_BEFORE_INTERVENTIONS = GLOBAL.CALIBRATION_DELAY + GLOBAL.DAYS_BEFORE_LOCKDOWN;

  GLOBAL.FIRST_PERIOD = optvals["FIRST_PERIOD"].as<double>();
  GLOBAL.SECOND_PERIOD = optvals["SECOND_PERIOD"].as<double>();
  GLOBAL.THIRD_PERIOD = optvals["THIRD_PERIOD"].as<double>();
  GLOBAL.FOURTH_PERIOD = optvals["FOURTH_PERIOD"].as<double>();

  GLOBAL.OE_SECOND_PERIOD = optvals["OE_SECOND_PERIOD"].as<double>();

  GLOBAL.CYCLIC_POLICY_TYPE
	= static_cast<Cycle_Type>(optvals["CYCLIC_POLICY_TYPE"].as<count_type>());
  GLOBAL.CYCLIC_POLICY_START_DAY = GLOBAL.NUM_DAYS_BEFORE_INTERVENTIONS +
	GLOBAL.FIRST_PERIOD + GLOBAL.SECOND_PERIOD;

  std::string output_dir(optvals["output_directory"].as<std::string>());
  GLOBAL.output_path = output_dir; //sk

  GLOBAL.input_base = optvals["input_directory"].as<std::string>();
  if(optvals["attendance_filename"].count()){
    GLOBAL.attendance_filename = optvals["attendance_filename"].as<std::string>();
    GLOBAL.IGNORE_ATTENDANCE_FILE = false;
  }
  //GLOBAL.IGNORE_ATTENDANCE_FILE = optvals["IGNORE_ATTENDANCE_FILE"].count();

  GLOBAL.USE_AGE_DEPENDENT_MIXING = optvals["USE_AGE_DEPENDENT_MIXING"].count();
  GLOBAL.SIGNIFICANT_EIGEN_VALUES = optvals["SIGNIFICANT_EIGEN_VALUES"].as<double>();
  GLOBAL.NUM_AGE_GROUPS = optvals["NUM_AGE_GROUPS"].as<count_type>();
  
  if(optvals["PROVIDE_INITIAL_SEED"].count()){
	//Initial seed was provided
	SEED_RNG_PROVIDED_SEED(optvals["PROVIDE_INITIAL_SEED"].as<count_type>());
  } else {
	SEED_RNG(); //No Initial seed was provided
  }
  if(optvals["PROVIDE_INITIAL_SEED_GRAPH"].count()){
	//Initial seed was provided
	SEED_RNG_GRAPH_PROVIDED_SEED(optvals["PROVIDE_INITIAL_SEED_GRAPH"].as<count_type>());
  } else {
	SEED_RNG_GRAPH(); //No Initial seed was provided
  }

  GLOBAL.LOCKED_COMMUNITY_LEAKAGE = optvals["LOCKED_COMMUNITY_LEAKAGE"].as<double>();
  GLOBAL.COMMUNITY_LOCK_THRESHOLD = optvals["COMMUNITY_LOCK_THRESHOLD"].as<double>();
  GLOBAL.LOCKED_NEIGHBORHOOD_LEAKAGE = optvals["LOCKED_NEIGHBORHOOD_LEAKAGE"].as<double>();
  GLOBAL.NEIGHBORHOOD_LOCK_THRESHOLD = optvals["NEIGHBORHOOD_LOCK_THRESHOLD"].as<double>();

  //Compute parametrs based on options
  GLOBAL.NUM_TIMESTEPS = GLOBAL.NUM_DAYS*GLOBAL.SIM_STEPS_PER_DAY;
  GLOBAL.INCUBATION_PERIOD_SCALE = GLOBAL.MEAN_INCUBATION_PERIOD*GLOBAL.SIM_STEPS_PER_DAY / GLOBAL.INCUBATION_PERIOD_SHAPE;

  GLOBAL.ASYMPTOMATIC_PERIOD = GLOBAL.MEAN_ASYMPTOMATIC_PERIOD*GLOBAL.SIM_STEPS_PER_DAY;
  GLOBAL.SYMPTOMATIC_PERIOD = GLOBAL.MEAN_SYMPTOMATIC_PERIOD*GLOBAL.SIM_STEPS_PER_DAY;
  GLOBAL.HOSPITAL_REGULAR_PERIOD = GLOBAL.MEAN_HOSPITAL_REGULAR_PERIOD*GLOBAL.SIM_STEPS_PER_DAY;
  GLOBAL.HOSPITAL_CRITICAL_PERIOD = GLOBAL.MEAN_HOSPITAL_CRITICAL_PERIOD*GLOBAL.SIM_STEPS_PER_DAY;
  
  GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD*GLOBAL.SIM_STEPS_PER_DAY;

  GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD*GLOBAL.SIM_STEPS_PER_DAY/GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_SHAPE;

  GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED1_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED1*GLOBAL.SIM_STEPS_PER_DAY /GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED1_SHAPE;// 2.29 days;

  GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED2_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED2*GLOBAL.SIM_STEPS_PER_DAY /GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_VACCINATED2_SHAPE;// 2.29 days;;

	GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_WANING_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD_WANING*GLOBAL.SIM_STEPS_PER_DAY /GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_WANING_SHAPE;// 2.29 days;;

	GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED*GLOBAL.SIM_STEPS_PER_DAY /GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED_SHAPE;// 2.29 days;;

	GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED2_SCALE=GLOBAL.MEAN_RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED2*GLOBAL.SIM_STEPS_PER_DAY /GLOBAL.RECOVERED_TO_SUSCEPTIBLE_PERIOD_BOOSTED2_SHAPE;// 2.29 days;;


  GLOBAL.MASK_ACTIVE = optvals["MASK_ACTIVE"].count();
  if(GLOBAL.MASK_ACTIVE){
    GLOBAL.MASK_FACTOR = optvals["MASK_FACTOR"].as<double>();
  }

  GLOBAL.MASK_START_DATE = GLOBAL.CALIBRATION_DELAY + optvals["MASK_START_DELAY"].as<double>(); //masks starts from April 9

  //initialise city bounding box co-ordinates with Hillsborough values. Will read from file later.
  GLOBAL.city_SW.lat = optvals["CITY_SW_LAT"].as<double>();
  GLOBAL.city_SW.lon = optvals["CITY_SW_LON"].as<double>();
  GLOBAL.city_NE.lat = optvals["CITY_NE_LAT"].as<double>();
  GLOBAL.city_NE.lon = optvals["CITY_NE_LON"].as<double>();
  GLOBAL.NBR_CELL_SIZE = optvals["NBR_CELL_SIZE"].as<double>();
  GLOBAL.ENABLE_CONTAINMENT = optvals["ENABLE_CONTAINMENT"].count();
  GLOBAL.ENABLE_NBR_CELLS = optvals["ENABLE_NBR_CELLS"].count();

//  std::cout<<"neghbors enabled or not in default drivesimul?"<<"\t"<<GLOBAL.ENABLE_NBR_CELLS<<std::endl;//----Shakir-Cheking neighbors.

  if(GLOBAL.ENABLE_NBR_CELLS){
	GLOBAL.ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT = optvals["ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT"].count();
  } else {
	GLOBAL.ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT = false;
  }
  GLOBAL.WARD_CONTAINMENT_THRESHOLD = optvals["WARD_CONTAINMENT_THRESHOLD"].as<count_type>();

  GLOBAL.ENABLE_TESTING = optvals["ENABLE_TESTING"].count();
  GLOBAL.ENABLE_NBR_CELLS = GLOBAL.ENABLE_NBR_CELLS || GLOBAL.ENABLE_CONTAINMENT;

  GLOBAL.TESTING_PROTOCOL
	= static_cast<Testing_Protocol>(optvals["TESTING_PROTOCOL"].as<count_type>());
  GLOBAL.testing_protocol_filename = optvals["testing_protocol_filename"].as<std::string>();

  GLOBAL.agent_load_file = optvals["agent_load_file"].as<std::string>();//sk

  GLOBAL.VACCINATION_EFFECTIVENESS1 = optvals["VACCINATION_EFFECTIVENESS1"].as<double>();
  GLOBAL.VACCINATION_EFFECTIVENESS2 = optvals["VACCINATION_EFFECTIVENESS2"].as<double>();
  GLOBAL.VACCINATION_EFFECTIVENESS_WANING = optvals["VACCINATION_EFFECTIVENESS_WANING"].as<double>();
  GLOBAL.VACCINATION_EFFECTIVENESS_BOOSTED = optvals["VACCINATION_EFFECTIVENESS_BOOSTED"].as<double>();
  GLOBAL.VACCINATION_EFFECTIVENESS_WANING2 = optvals["VACCINATION_EFFECTIVENESS_WANING2"].as<double>();
  GLOBAL.VACCINATION_EFFECTIVENESS_BOOSTED2 = optvals["VACCINATION_EFFECTIVENESS_BOOSTED2"].as<double>();

  GLOBAL.TIME_ALPHA = optvals["TIME_ALPHA"].as<double>();
  GLOBAL.TIME_DELTA = optvals["TIME_DELTA"].as<double>();
  GLOBAL.TIME_OMICRON = optvals["TIME_OMICRON"].as<double>();
  GLOBAL.TIME_OMICRON_NEW = optvals["TIME_OMICRON_NEW"].as<double>();
  GLOBAL.TIME_OMICRON_BA4 = optvals["TIME_OMICRON_BA4"].as<double>();
  GLOBAL.TIME_OMICRON_BA5 = optvals["TIME_OMICRON_BA5"].as<double>();

  GLOBAL.INFECTIOUSNESS_ALPHA = optvals["INFECTIOUSNESS_ALPHA"].as<double>();
  GLOBAL.INFECTIOUSNESS_DELTA = optvals["INFECTIOUSNESS_DELTA"].as<double>();
  GLOBAL.INFECTIOUSNESS_OMICRON = optvals["INFECTIOUSNESS_OMICRON"].as<double>();
  GLOBAL.INFECTIOUSNESS_OMICRON_NEW = optvals["INFECTIOUSNESS_OMICRON_NEW"].as<double>();
  GLOBAL.INFECTIOUSNESS_OMICRON_BA4 = optvals["INFECTIOUSNESS_OMICRON_BA4"].as<double>();
  GLOBAL.INFECTIOUSNESS_OMICRON_BA5 = optvals["INFECTIOUSNESS_OMICRON_BA5"].as<double>();

  GLOBAL.VIRULENT_NEW_STRAIN = optvals["VIRULENT_NEW_STRAIN"].as<double>();
  GLOBAL.VIRULENT_NEW_ALPHA = optvals["VIRULENT_NEW_ALPHA"].as<double>();
  GLOBAL.VIRULENT_NEW_DELTA = optvals["VIRULENT_NEW_DELTA"].as<double>();
  GLOBAL.VIRULENT_NEW_OMICRON = optvals["VIRULENT_NEW_OMICRON"].as<double>();
  GLOBAL.VIRULENT_NEW_OMICRON_NEW = optvals["VIRULENT_NEW_OMICRON_NEW"].as<double>();
  GLOBAL.VIRULENT_NEW_OMICRON_BA4 = optvals["VIRULENT_NEW_OMICRON_BA4"].as<double>();
  GLOBAL.VIRULENT_NEW_OMICRON_BA5 = optvals["VIRULENT_NEW_OMICRON_BA5"].as<double>();

  
  GLOBAL.REINFECTION_ALPHA = optvals["REINFECTION_ALPHA"].as<double>();
  GLOBAL.REINFECTION_DELTA = optvals["REINFECTION_DELTA"].as<double>();
  GLOBAL.REINFECTION_OMICRON = optvals["REINFECTION_OMICRON"].as<double>();
  GLOBAL.REINFECTION_OMICRON_NEW = optvals["REINFECTION_OMICRON_NEW"].as<double>();
  GLOBAL.REINFECTION_OMICRON_BA4 = optvals["REINFECTION_OMICRON_BA4"].as<double>();
  GLOBAL.REINFECTION_OMICRON_BA5 = optvals["REINFECTION_OMICRON_BA5"].as<double>();

  GLOBAL.FRACTION_NEW_ALPHA = optvals["FRACTION_NEW_ALPHA"].as<double>();
  GLOBAL.FRACTION_NEW_DELTA = optvals["FRACTION_NEW_DELTA"].as<double>();
  GLOBAL.FRACTION_NEW_OMICRON = optvals["FRACTION_NEW_OMICRON"].as<double>();
  GLOBAL.FRACTION_NEW_OMICRON_NEW = optvals["FRACTION_NEW_OMICRON_NEW"].as<double>();
  GLOBAL.FRACTION_NEW_OMICRON_BA4 = optvals["FRACTION_NEW_OMICRON_BA4"].as<double>();
  GLOBAL.FRACTION_NEW_OMICRON_BA5 = optvals["FRACTION_NEW_OMICRON_BA5"].as<double>();

  GLOBAL.FRACTION_SUSCEPTIBLE_ALPHA = optvals["FRACTION_SUSCEPTIBLE_ALPHA"].as<double>();
  GLOBAL.FRACTION_SUSCEPTIBLE_DELTA = optvals["FRACTION_SUSCEPTIBLE_DELTA"].as<double>();
  GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON = optvals["FRACTION_SUSCEPTIBLE_OMICRON"].as<double>();
  GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_NEW = optvals["FRACTION_SUSCEPTIBLE_OMICRON_NEW"].as<double>();
  GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_BA4 = optvals["FRACTION_SUSCEPTIBLE_OMICRON_BA4"].as<double>();
  GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_BA5 = optvals["FRACTION_SUSCEPTIBLE_OMICRON_BA5"].as<double>();

  GLOBAL.FRACTION_NEW_STRAIN = optvals["FRACTION_NEW_STRAIN"].as<double>();

//printf("I am in drive simul\n");


  if(GLOBAL.input_base != "" && GLOBAL.input_base[GLOBAL.input_base.size() - 1] != '/'){ 
	GLOBAL.input_base += '/';
	//Make sure the path of the input_base
	//directory is terminated by a "/"
  }
  // sk ///////////////////////////////////////////
  //settings for cohorts
 GLOBAL.ENABLE_COHORTS = optvals["ENABLE_COHORTS"].count();
 GLOBAL.ISOLATE_COHORTS = optvals["ISOLATE_COHORTS"].count();
 GLOBAL.BETA_COHORT = optvals["BETA_COHORT"].as<double>();
 GLOBAL.crowding_factor = optvals["CROWDING_FACTOR_COHORTS"].as<double>();
 GLOBAL.taking_train_fraction = optvals["FRACTION_IN_TRAINS_COHORTS"].as<double>();
 GLOBAL.COHORT_SIZE = optvals["COHORT_SIZE"].as<double>();
 GLOBAL.COHORT_SEVERITY_FRACTION = optvals["COHORT_SEVERITY_FRACTION"].as<double>();
 GLOBAL.COHORT_STRATEGY = static_cast<cohort_strategy>(optvals["COHORT_STRATEGY"].as<count_type>());
 GLOBAL.ONE_OFF_TRAVELERS_RATIO = optvals["ONE_OFF_TRAVELERS_RATIO"].as<double>();
 
//  std::cout<<"GLOBAL.ENABLE_COHORTS" << GLOBAL.ENABLE_COHORTS << std::endl;

// store or load state.
 GLOBAL.STORE_STATE_TIME_STEP = optvals["STORE_STATE_TIME_STEP"].as<count_type>();
 GLOBAL.LOAD_STATE_TIME_STEP = optvals["LOAD_STATE_TIME_STEP"].as<count_type>();

//////////////////////////////////////////////
  //Initialize the attendance probability
  initialize_office_attendance();

  //Initialize output folders
  // gnuplot gnuplot(output_dir);

  //Run simulations
  auto plot_data = run_simulation();
  
  
  std::string csv_file_name = "infections_from_new_strain"+std::to_string(GLOBAL.START_DAY)+"_"+std::to_string(GLOBAL.NUM_DAYS)+".csv";
  //std::string csv_file_path = output_dir + "/" + csv_file_name;
  std::string csv_file_path = output_dir + csv_file_name; // make sure --output_directory ends with /
  std::cerr << "drive_simulation: Writing to " << csv_file_path << "\n";   
  output_sensitivities_csv(plot_data.logger1, csv_file_path);


  //Start output
//  output_global_params(output_dir);

  output_csv_files(output_dir, 
  // gnuplot, 
  plot_data);
  return 0;
}
