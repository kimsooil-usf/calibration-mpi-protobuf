//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include "models.h"
#include "initializers.h"
#include "updates.h"
#include "simulator.h"
#include "testing.h"
#include "cohorts.h"

using std::vector;
using std::string;

#if defined DEBUG || defined TIMING
#include <iostream>
#include <cstdlib>
using std::cerr;
#endif

#ifdef TIMING
#include <chrono>

template <class T>
auto duration(const std::chrono::time_point<T>& start, const std::chrono::time_point<T>& end){
  return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}
#endif

// sk: not sure what it is / remove if we don't modelling commuting-by-train
namespace { // Anon namespace for local functions
int num_coaches(const std::unordered_map<count_type, std::vector<train_coach>>& coaches) {
	int i =0;
	for (const auto& c : coaches) {
		i += c.second.size();
	}
	return i;
}
}

plot_data_struct run_simulation(){
#ifdef TIMING
//   cerr << "simulator: starting JSON read\n";
  auto start_time = std::chrono::high_resolution_clock::now();
#endif
//std::cout<<"Enabled neighbors? before intiating"<<"\t"<<GLOBAL.ENABLE_NBR_CELLS<<std::endl;//---delete shakir
 
  auto homes = init_homes();
  //auto mask = init_mask();//---disabled temporariliy --shakir
  
  auto workplaces = init_workplaces();
  auto communities = init_community();
  auto nodes = init_nodes();
  auto nbr_cells = init_nbr_cells();
  auto intv_params = init_intervention_params();
  auto testing_protocol_file_read = init_testing_protocol();
//std::cout<<"\n"<<"comunity size"<<"\t"<<communities.size()<<std::endl;

//std::cout<<"\n"<<"neighbors size after initiating"<<"\t"<<nbr_cells.size()<<std::endl;
	// sk /////////////////////////////
	auto train_loader = init_TrainLoader();
	auto cohorts = make_cohorts(nodes, GLOBAL.COHORT_SIZE, train_loader);
	// std::cout<<"\nNum Cohorts: " << get_num_cohorts(cohorts)<< std::endl;
	auto one_off_cohorts = make_oneoff_cohorts(nodes, train_loader);
	// std::cout<<"One Off Cohorts: " << get_num_cohorts(one_off_cohorts)<< std::endl;
	merge_cohorts(cohorts, one_off_cohorts);
	// std::cout<<"Merged Cohorts: " << get_num_cohorts(cohorts)<< std::endl;
	bool coaches_created = false;
	std::unordered_map<count_type, std::vector<train_coach>> train_coaches_am;
	std::unordered_map<count_type, std::vector<train_coach>> train_coaches_pm;
	///////////////////////////////


  auto community_dist_matrix = compute_community_distances(communities);
  auto community_fk_matrix = compute_community_distances_fkernel(community_dist_matrix);

  svd home_age_matrix,
	school_age_matrix,
	workplace_age_matrix,
	community_age_matrix;

  if(GLOBAL.USE_AGE_DEPENDENT_MIXING){
	home_age_matrix = init_home_age_interaction_matrix();
	school_age_matrix = init_school_age_interaction_matrix();
	workplace_age_matrix = init_workplace_age_interaction_matrix();
	community_age_matrix = init_community_age_interaction_matrix();
  }
  
#ifdef TIMING
    auto end_time = std::chrono::high_resolution_clock::now();
	// cerr << "simulator: time for JSON reads (ms): " << duration(start_time, end_time) << "\n";
	start_time = std::chrono::high_resolution_clock::now();
#endif

  assign_individual_home_community(nodes, homes, workplaces, communities);
  //assign_individual_home_community must be called before assign_homes_nbr_cell
  assign_homes_nbr_cell(homes,nbr_cells);
  assign_individual_projects(workplaces, nodes);
  assign_household_community(communities, nodes, homes);
  assign_household_random_community(homes, communities);
	//   assign_inter_cohort(cohorts, nodes); //TODO[v2]: Use when inter-cohort interactions are enabled // sk: just copied for future

  compute_scale_homes(homes);
  compute_scale_workplaces(workplaces);
  compute_scale_communities(nodes, communities);
  compute_scale_random_community(homes, nodes);
  compute_scale_nbr_cells(nodes, nbr_cells, homes);

	compute_scale_intra_cohorts(cohorts, nodes); // sk: not sure if it always/sometimes being used


  double travel_fraction = 0;
  std::vector<double> travel_fraction_higher;
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);
	travel_fraction_higher.push_back(0);

  //double *travel_fraction_higher=(double*) malloc(sizeof(double)*7);
	// sk /////////////////////////////////////// this block exists only in latest version
	//This needs to be done after the initilization.

	int home_ward_infected[GLOBAL.num_wards];
	//int work_ward_infected[GLOBAL.num_wards];

	for (int nwards = 0; nwards < GLOBAL.num_wards; nwards++){
		home_ward_infected[nwards] = 0;
		//work_ward_infected[nwards] = 0;
	}
	/////////////////////////////////////////////



  //This needs to be done after the initilization.
  plot_data_struct plot_data;
  plot_data.nums =
	{
	 {"num_infected", {}},
	 {"num_exposed", {}},
	 {"num_hospitalised", {}},
	{"num_unvaccinated", {}},
		{"num_vaccinated1", {}},
		{"num_vaccinated2", {}},
		{"num_waning", {}},
		{"num_boosted", {}},
	{"num_boosted2", {}},


		{"num_infected_age_group_1", {}},
		{"num_infected_age_group_2", {}},
		{"num_infected_age_group_3", {}},
		{"num_infected_age_group_4", {}},
		{"num_infected_age_group_5", {}},
		{"num_infected_age_group_6", {}},
		{"num_infected_age_group_7", {}},
		{"num_infected_age_group_8", {}},
		{"num_infected_age_group_9", {}},
		{"num_infected_age_group_10", {}},


		{"num_hospitalised_age_group_1", {}},
		{"num_hospitalised_age_group_2", {}},
		{"num_hospitalised_age_group_3", {}},
		{"num_hospitalised_age_group_4", {}},
		{"num_hospitalised_age_group_5", {}},
		{"num_hospitalised_age_group_6", {}},
		{"num_hospitalised_age_group_7", {}},
		{"num_hospitalised_age_group_8", {}},
		{"num_hospitalised_age_group_9", {}},
		{"num_hospitalised_age_group_10", {}},

		{"num_dead_age_group_1", {}},
		{"num_dead_age_group_2", {}},
		{"num_dead_age_group_3", {}},
		{"num_dead_age_group_4", {}},
		{"num_dead_age_group_5", {}},
		{"num_dead_age_group_6", {}},
		{"num_dead_age_group_7", {}},
		{"num_dead_age_group_8", {}},
		{"num_dead_age_group_9", {}},
		{"num_dead_age_group_10", {}},

		{"num_infected_race_group_1", {}},
		{"num_infected_race_group_2", {}},
		{"num_infected_race_group_3", {}},
		{"num_infected_race_group_4", {}},
		{"num_infected_race_group_5", {}},
		{"num_infected_race_group_6", {}},
		{"num_infected_race_group_7", {}},

		{"num_infected_income_group_1", {}},
		{"num_infected_income_group_2", {}},
		{"num_infected_income_group_3", {}},
		{"num_infected_income_group_4", {}},
		{"num_infected_income_group_5", {}},
		{"num_infected_income_group_6", {}},
		{"num_infected_income_group_7", {}},
		{"num_infected_income_group_8", {}},

		{"num_infected_ethnicity_group_1", {}},
		{"num_infected_ethnicity_group_2", {}},


		{"num_infected_gender_group_1", {}},
		{"num_infected_gender_group_2", {}},

		{"num_hospitalised_race_group_1", {}},
		{"num_hospitalised_race_group_2", {}},
		{"num_hospitalised_race_group_3", {}},
		{"num_hospitalised_race_group_4", {}},
		{"num_hospitalised_race_group_5", {}},
		{"num_hospitalised_race_group_6", {}},
		{"num_hospitalised_race_group_7", {}},

		{"num_hospitalised_income_group_1", {}},
		{"num_hospitalised_income_group_2", {}},
		{"num_hospitalised_income_group_3", {}},
		{"num_hospitalised_income_group_4", {}},
		{"num_hospitalised_income_group_5", {}},
		{"num_hospitalised_income_group_6", {}},
		{"num_hospitalised_income_group_7", {}},
		{"num_hospitalised_income_group_8", {}},

		{"num_hospitalised_ethnicity_group_1", {}},
		{"num_hospitalised_ethnicity_group_2", {}},

		{"num_hospitalised_gender_group_1", {}},
		{"num_hospitalised_gender_group_2", {}},


		{"num_dead_race_group_1", {}},
		{"num_dead_race_group_2", {}},
		{"num_dead_race_group_3", {}},
		{"num_dead_race_group_4", {}},
		{"num_dead_race_group_5", {}},
		{"num_dead_race_group_6", {}},
		{"num_dead_race_group_7", {}},

		{"num_dead_income_group_1", {}},
		{"num_dead_income_group_2", {}},
		{"num_dead_income_group_3", {}},
		{"num_dead_income_group_4", {}},
		{"num_dead_income_group_5", {}},
		{"num_dead_income_group_6", {}},
		{"num_dead_income_group_7", {}},
		{"num_dead_income_group_8", {}},

		{"num_dead_ethnicity_group_1", {}},
		{"num_dead_ethnicity_group_2", {}},

		{"num_dead_gender_group_1", {}},
	{"num_dead_gender_group_2", {}},
	

	 {"num_symptomatic", {}},
	 {"num_critical", {}},
	 {"num_fatalities", {}},
	 {"num_recovered", {}},
	 {"num_affected", {}},
	 {"num_cases", {}},
	 {"num_cumulative_hospitalizations", {}},
	 {"num_cumulative_infective", {}}
	};
  for(auto& elem: plot_data.nums){
	elem.second.reserve(GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS);
  }
  plot_data.nums["csvContent"] = {};
  plot_data.nums["csvContent"].reserve((GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS) * GLOBAL.num_communities);

  plot_data.susceptible_lambdas =
	{
	 {"susceptible_lambda", {}},
	 {"susceptible_lambda_H", {}},
	 {"susceptible_lambda_W", {}},
	 {"susceptible_lambda_C", {}},
	 {"susceptible_lambda_T", {}},
	 {"susceptible_lambda_PROJECT", {}},
	 {"susceptible_lambda_NBR_CELL", {}},
	 {"susceptible_lambda_RANDOM_COMMUNITY", {}}
	};

  plot_data.total_lambda_fractions =
	{
	 {"total_fraction_lambda_H", {}},
	 {"total_fraction_lambda_W", {}},
	 {"total_fraction_lambda_C", {}},
	 {"total_fraction_lambda_T", {}},
	 {"total_fraction_lambda_PROJECT", {}},
	 {"total_fraction_lambda_NBR_CELL", {}},
	 {"total_fraction_lambda_RANDOM_COMMUNITY", {}}
	};

  plot_data.mean_lambda_fractions =
	{
	 {"mean_fraction_lambda_H", {}},
	 {"mean_fraction_lambda_W", {}},
	 {"mean_fraction_lambda_C", {}},
	 {"mean_fraction_lambda_T", {}},
	 {"mean_fraction_lambda_PROJECT", {}},
	 {"mean_fraction_lambda_NBR_CELL", {}},
	 {"mean_fraction_lambda_RANDOM_COMMUNITY", {}}
	};

  plot_data.cumulative_mean_lambda_fractions = 
	{
	 {"cumulative_mean_fraction_lambda_H", {}},
	 {"cumulative_mean_fraction_lambda_W", {}},
	 {"cumulative_mean_fraction_lambda_C", {}},
	 {"cumulative_mean_fraction_lambda_T", {}},
	 {"cumulative_mean_fraction_lambda_PROJECT", {}},
	 {"cumulative_mean_fraction_lambda_NBR_CELL", {}},
	 {"cumulative_mean_fraction_lambda_RANDOM_COMMUNITY", {}}
	};

  plot_data.quarantined_stats =
	{
	 {"quarantined_stats", {}},
	 {"curtailment_stats", {}}
	};
	   
  plot_data.disease_label_stats = 
	{
		{"disease_label_stats", {}},
	};
  for(auto& elem: plot_data.susceptible_lambdas){
	elem.second.reserve(GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS);
  }

#ifdef TIMING
  end_time = std::chrono::high_resolution_clock::now();
//   cerr << "simulator: time for setup after JSON reads (ms): " << duration(start_time, end_time) << "\n";

//   cerr << "simulator: starting simulation\n";
  start_time = std::chrono::high_resolution_clock::now();
#endif
  lambda_incoming_data total_lambda_fraction_data;
  lambda_incoming_data mean_lambda_fraction_data;
  lambda_incoming_data cumulative_mean_lambda_fraction_data;
  count_type num_cases = 0; // Total number of agents who have progessed to symptomatic so far
  count_type quarantined_num_cases = 0;
  count_type num_cumulative_hospitalizations = 0; //Total number of agents who have had to go to the hospital so far
  count_type num_cumulative_infective = 0; //Total number of people who have progressed to the infective state so far

  count_type num_total_infections = 0;
  //Total number of individuals who have become infected via transmission so far
  //This does not included the initially seeded infections

  std::vector<long double> infections_by_new_infectives(GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS, 0);
  //For keeping track of infections ascribed to agents that became infective at
  //each time
  
  std::vector<count_type> new_strain_candidates;
  new_strain_candidates.reserve(GLOBAL.num_people);
  std::string logging1 = "Time_step,total_new_infections,new_strain0_infections,new_strain1_infections,new_strain2_infections,new_strain3_infections,new_strain4_infections,new_strain5_infections,new_strain6_infections,active_patients,active_patients_new_strain2,active_patients_new_strain1,active_patients_new_strain2,active_patients_new_strain3,active_patients_new_strain4,active_patients_new_strain5,active_patients_new_strain6" ;
  std::vector<std::string> logger1;
  logger1.push_back(logging1);



//----age group less than 20
  std::vector<count_type> new_vaccinated1_candidates_LE20;
  new_vaccinated1_candidates_LE20.reserve(GLOBAL.num_people);

  std::vector<count_type> new_vaccinated2_candidates_LE20;
  new_vaccinated2_candidates_LE20.reserve(GLOBAL.num_people);
  

  std::vector<count_type> new_waning_candidates_LE20;
  new_waning_candidates_LE20.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted_candidates_LE20;
  new_boosted_candidates_LE20.reserve(GLOBAL.num_people);

  std::vector<count_type> new_waning2_candidates_LE20;
  new_waning_candidates_LE20.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted2_candidates_LE20;
  new_boosted_candidates_LE20.reserve(GLOBAL.num_people);


//age group between 30 and 60
  std::vector<count_type> new_vaccinated1_candidates_30_to_60;
  new_vaccinated1_candidates_30_to_60.reserve(GLOBAL.num_people);

  std::vector<count_type> new_vaccinated2_candidates_30_to_60;
  new_vaccinated2_candidates_30_to_60.reserve(GLOBAL.num_people);
  

  std::vector<count_type> new_waning_candidates_30_to_60;
  new_waning_candidates_30_to_60.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted_candidates_30_to_60;
  new_boosted_candidates_30_to_60.reserve(GLOBAL.num_people);

  std::vector<count_type> new_waning2_candidates_30_to_60;
  new_waning_candidates_30_to_60.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted2_candidates_30_to_60;
  new_boosted_candidates_30_to_60.reserve(GLOBAL.num_people);

//age group above 60
  std::vector<count_type> new_vaccinated1_candidates_GT60;
  new_vaccinated1_candidates_GT60.reserve(GLOBAL.num_people);

  std::vector<count_type> new_vaccinated2_candidates_GT60;
  new_vaccinated2_candidates_GT60.reserve(GLOBAL.num_people);
  

  std::vector<count_type> new_waning_candidates_GT60;
  new_waning_candidates_GT60.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted_candidates_GT60;
  new_boosted_candidates_GT60.reserve(GLOBAL.num_people);

  std::vector<count_type> new_waning2_candidates_GT60;
  new_waning_candidates_GT60.reserve(GLOBAL.num_people);


  std::vector<count_type> new_boosted2_candidates_GT60;
  new_boosted_candidates_GT60.reserve(GLOBAL.num_people);

  //----------------Age group wise vaccination variable declared


  const auto NUM_PEOPLE = GLOBAL.num_people;

//---Vaccination timings---------------//

	double nVacc1=10;//---->rateVacc1=x*time_step+b. The number of people vaccinated changes with time.
	double nVacc2=5;//---->These numbers need to be fixed after discussion.
	double nWaning=3;//---->these numbers need to be fixed after discussion with Ken.
	double nBoosted=10,z;
	double nboosted2=200;//---->these numbers need to be  fixed after discussion with Ken.
	double vaccFn,vaccFn2,vaccFn3,vaccFn4,Time_VaccStart,TimeVacc2Start,TimeBoostStart,TimeVaccConst,TimeVacc2Const;
	double TimeBoostStart1,TimeBoostConst,TimeBoost2Start,TimeBoost2Const;


	Time_VaccStart=306*GLOBAL.SIM_STEPS_PER_DAY;//Jan 1 2021 is 306 days away from March 1, 2020.
	TimeVacc2Start=(12+306)*GLOBAL.SIM_STEPS_PER_DAY;//Jan 1 2021 is 12+ 306 days away from March 1, 2020.
	TimeBoostStart=518*GLOBAL.SIM_STEPS_PER_DAY;//1 august 2021 is 518 days away from March 1, 2020.
	TimeBoostStart1=654*GLOBAL.SIM_STEPS_PER_DAY;//15 December 2021 is 654 days away from 1 March 2020;Between 1 august and december 15 2021 second doses are 1455 per day. after december 15 2021 it is the third polynomial;
	
	TimeVaccConst=760*GLOBAL.SIM_STEPS_PER_DAY;//Maarch 31, 2022 is 760 days away from march 1 2020.

	TimeVacc2Const=760*GLOBAL.SIM_STEPS_PER_DAY;//Maarch 31, 2022 is 760 days away from march 1 2020. (days since vacc start is 760-306=454 days)

	TimeBoostConst=679*GLOBAL.SIM_STEPS_PER_DAY;//January 9, 2022 is 679 days away from march 1 2020. after this time third doses are a constant;
	
	TimeBoost2Const=810*GLOBAL.SIM_STEPS_PER_DAY;//May 10, 2022 is 800 days from 1st march 2022.

//---Vaccination timings Shakir------------//


//   for(count_type time_step = GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY; time_step < GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS; ++time_step){
// #ifdef DEBUG
// 	auto start_time_timestep = std::chrono::high_resolution_clock::now();
// #endif
	// sk /////////////////////////////////////
	count_type time_step_start = 0;
	// std::cout << "LOAD_STATE_TIME_STEP " << GLOBAL.LOAD_STATE_TIME_STEP
	// 					<< "\nSTORE_STATE_TIME_STEP " << GLOBAL.STORE_STATE_TIME_STEP
	// 					<< "\n\n" << std::endl;
	time_step_start = GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY;
	count_type time_step_end = GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS;
	#ifdef ENABLE_PROTO
	if (GLOBAL.LOAD_STATE_TIME_STEP > 0) {
		if (loadAgentsInfo(nodes)) {
			time_step_start = GLOBAL.LOAD_STATE_TIME_STEP;
			std::cout << " Loaded state at "+std::to_string(time_step_start) << std::endl;
		} else {
			std::cout << "ERROR: Loading state failed" << std::endl;
			return plot_data;
		}
	}
	else {
		std::cout << "NOT loading state...(Start from Day 0)" << std::endl;
	}
	#endif

	//for (count_type time_step = time_step_start; time_step < GLOBAL.NUM_TIMESTEPS; ++time_step){ // sk
	/////////////////////////////////////////////
	
  for(count_type time_step = time_step_start; time_step < time_step_end; ++time_step){

int new_day=time_step%GLOBAL.SIM_STEPS_PER_DAY;
//----Time depdendent mask factor
//GLOBAL.MASK_FACTOR=0.8;
//GLOBAL.MASK_ACTIVE=bernoulli(mask[int(time_step/GLOBAL.SIM_STEPS_PER_DAY)].maskcompliance);//put the time dependent function here indicating mask active or not.

//GLOBAL.COMPLIANCE_PROBABILITY=0.;//mask[int(time_step/GLOBAL.SIM_STEPS_PER_DAY)].maskcompliance;
//-------------age group less than 20----------------------------//
		new_vaccinated1_candidates_LE20.clear();
	
		new_vaccinated2_candidates_LE20.clear();

		new_waning_candidates_LE20.clear();

		new_boosted_candidates_LE20.clear();

		new_waning2_candidates_LE20.clear();

		new_boosted2_candidates_LE20.clear();
//-------------age group 30 to 60--------------------------------//
		new_vaccinated1_candidates_30_to_60.clear();
	
		new_vaccinated2_candidates_30_to_60.clear();

		new_waning_candidates_30_to_60.clear();

		new_boosted_candidates_30_to_60.clear();		

		new_waning2_candidates_30_to_60.clear();

		new_boosted2_candidates_30_to_60.clear();	
//-------------age group above 60--------------------------------//

		new_vaccinated1_candidates_GT60.clear();
	
		new_vaccinated2_candidates_GT60.clear();

		new_waning_candidates_GT60.clear();

		new_boosted_candidates_GT60.clear();

		new_waning2_candidates_GT60.clear();

		new_boosted2_candidates_GT60.clear();
//--------------------------------------------------------------//		
	
	total_lambda_fraction_data.set_zero();
	mean_lambda_fraction_data.set_zero();

	count_type num_new_infections = 0;
	count_type num_new_strain0_infections = 0,num_new_strain1_infections = 0,num_new_strain2_infections = 0,num_new_strain3_infections = 0,num_new_strain4_infections = 0,num_new_strain5_infections = 0,num_new_strain6_infections = 0;

	// sk //////////////////////
		#ifdef ENABLE_PROTO
		if (GLOBAL.STORE_STATE_TIME_STEP > 0 && time_step == GLOBAL.STORE_STATE_TIME_STEP) {
			storeAgentsInfo(nodes);
			std::cout << "Stored state at "+GLOBAL.STORE_STATE_TIME_STEP << std::endl;
		}
		#endif
	//////////////////////////////
    if (new_day==0){
      for(count_type j = 0; j < GLOBAL.num_people; ++j){
        nodes[j].attending =
          bernoulli(std::min(communities[nodes[j].community].w_c,
							 nodes[j].neighborhood_access_factor)
                    * nodes[j].get_attendance_probability(time_step));
		nodes[j].forced_to_take_train = GLOBAL.TRAINS_RUNNING
		  && bernoulli(GLOBAL.FRACTION_FORCED_TO_TAKE_TRAIN);
		  //std::cout<<"node is forced to take train"<<"\t"<<nodes[j].forced_to_take_train<<std::endl;
      }
    }

          
	//#pragma omp parallel for
	//
	// Since update_infection uses a random number generator with
	// global state, parallelizing this loop is not straightforward.
	// Puttting the generator in a critical section can keep it
	// correct, but slows down the code too much.

//----Number to be vaccinated is calculated here:Shakir------//
	if(time_step<Time_VaccStart){
		vaccFn=0;
	}
	if(time_step>=Time_VaccStart && time_step<=TimeVaccConst)
	{
	double tv0=time_step/GLOBAL.SIM_STEPS_PER_DAY-306;
	z=((tv0)-259)/149.4;

	vaccFn = (730.4*pow(z,5)-919.8*pow(z,4)-1855*pow(z,3)+2389*pow(z,2)-902.5*z+1398.0);
	}	

	if(time_step>=TimeVaccConst){
		vaccFn=225.0;
	}

	if(time_step<TimeVacc2Start){
		vaccFn2=0;

	}
	if(time_step>=TimeVacc2Start && time_step<=TimeVacc2Const)
	{
	double tv1=time_step/GLOBAL.SIM_STEPS_PER_DAY-306-12;	
	double z1=((tv1)-259)/149.4;
	vaccFn2 = (436.5*pow(z1,5)-908.4*pow(z1,4)-555*pow(z1,3)+1957*pow(z1,2)-1711*z1+1476);

	}	

	if(time_step>TimeVacc2Const){
		vaccFn2=225;
	}
	

	if(time_step<TimeBoostStart){
		vaccFn3=0;
	}
	if(time_step>=TimeBoostStart && time_step<=TimeBoostStart1){
	double tv2=time_step/GLOBAL.SIM_STEPS_PER_DAY-518;
	double z2=((tv2)-83.5)/48.06;
	vaccFn3 = (191.6*pow(z2,5)-156.2*pow(z2,4)-809.3*pow(z2,3)+958*pow(z2,2)-92.61*z2+291.2);

	}
	if(time_step>TimeBoostStart1 && time_step<=TimeBoostConst){
		vaccFn3=1455.0;
	}

	if(time_step>=TimeBoostConst){
		vaccFn3=225;
	}

	if(vaccFn2<0)
	vaccFn2=0;

	if(time_step<TimeBoost2Const){
		vaccFn4=0;
	}
	if(time_step>=TimeBoost2Const){
		vaccFn4=225;
	}

//----Number to be vaccinated is calculated ends here: Shakir------//

// if t < 0 (January 1st, 2021), daily first doses = 0/day
// if 0 < t < 457,  daily doses = first_polynomial(t)
// if t >  457 (March 31 2022), daily first doses = 225/day.

// Second Dose:
// if t < 12 (January 12th 2021), daily second doses = 0/day
// if 12 < t < 456, daily doses = second_polynomial(t) t = 12 = January 12th 2021
// if t > 457 (March 31 2022), daily second doses = 225/day

// Third Dose (This one is a bit different due to a dump of vaccinations in the dataset)
// if t < August 1st, 2021, daily third doses = 0/day
// if August 1st 2021 < t < December 15th 2021, daily third doses = 1455/day
// if t = 0 = December 15th 2021, daily third doses = third_polynomial(t)
// if t >138 (January 9th, 2022) daily third doses = 225/day

count_type DAY_84=84;
count_type FIVE_MONTH=5*30*GLOBAL.SIM_STEPS_PER_DAY;
count_type SIX_MONTH=6*30*GLOBAL.SIM_STEPS_PER_DAY;

	for(count_type j = 0; j < NUM_PEOPLE; ++j){
		int age_index=nodes[j].age_index;
		nodes[j].new_vaccinated1=false;
		nodes[j].new_vaccinated2=false;
		nodes[j].new_waning=false;
		nodes[j].new_boosted=false;
		nodes[j].new_boosted2=false;
// int get_age_index(int age){
//   // Determine age category of individual.
//   if(age < 10) {
//     return 0;
//   } else if(age < 20) {
//     return 1;
//   } else if(age < 30) {
//     return 2;
//   } else if(age < 40) {
//     return 3;
//   } else if(age < 50) {
//     return 4;
//   } else if(age < 60) {
//     return 5;
//   } else if(age < 70) {
//     return 6;
//   } else if(age < 80) {
//     return 7;
//   } else {
//     return 8;
//   }
// }//----->Hint for age_index variable


//-----Hint for age dependent vaccinations-----//

//      if age >= 65:
//            Proportion of vaccine doses= 0.5
//      If age >= 30 and age < 65:
//            Proportion of vaccine doses=0.3
//      if age >= 18 and age < 30:
//             Proportion of vaccine doses= .2


			//---Vaccinating according to multiple doses scheme.
		if (time_step>=Time_VaccStart  && new_day==0)// && nodes[j].vaccinated1==false)//when start date is 1st March 2020 and vaccination start is 1 Jan 2021.
		{

         if(age_index<=1 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery != Progression::vaccinated1) )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_vaccinated1_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here

				}

			}
		    if(age_index>=2 && age_index<=5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery != Progression::vaccinated1) )){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_vaccinated1_candidates_30_to_60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
			}
		    
		    if(age_index>5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery != Progression::vaccinated1) )){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_vaccinated1_candidates_GT60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
			}
//-------------------New Vaccinated2----------------------------------------------------//
		if(age_index<=1 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated1) && time_step-nodes[j].time_at_vaccine1>=DAY_84 )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_vaccinated2_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here

				}

			}
		    if(age_index>=2 && age_index<=5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated1) && time_step-nodes[j].time_at_vaccine1>=DAY_84)){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_vaccinated2_candidates_30_to_60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
			}
		    
		    if(age_index>5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated1) && time_step-nodes[j].time_at_vaccine1>=DAY_84)){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_vaccinated2_candidates_GT60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
					
			}

//-------------------Waning 1---------------------------------------------------//

		if((nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated2) && time_step-nodes[j].time_at_vaccine2>=FIVE_MONTH )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_waning_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here
				}

			}

//------------------Boosted candidates--------------------------------------------//

		if(age_index<=1 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated2) && time_step-nodes[j].time_at_vaccine2>=SIX_MONTH )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_boosted_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here

				}

			}
		    if(age_index>=2 && age_index<=5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated2) && time_step-nodes[j].time_at_vaccine2>=SIX_MONTH)){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_boosted_candidates_30_to_60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
			}
		    
		    if(age_index>5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::vaccinated2) && time_step-nodes[j].time_at_vaccine2>=SIX_MONTH){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_boosted_candidates_GT60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
					
			}
//-------------------Waning and boosted 2-------------------------------------//

//-------------------Waning 2---------------------------------------------------//

		if((nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::boosted) && time_step-nodes[j].time_at_boosted>=FIVE_MONTH )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_waning2_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here
				}

			}

//------------------Boosted2 candidates--------------------------------------------//

		if(age_index<=1 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::boosted) && time_step-nodes[j].time_at_boosted>=SIX_MONTH )){//----->age between 1 and 10
			bool vacc1=false;
			vacc1=true;
				if(vacc1)
				{			
				new_boosted_candidates_LE20.push_back(j);//j-th individula is vaccinated. Choose age distribution here

				}

			}
		    if(age_index>=2 && age_index<=5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::boosted) && time_step-nodes[j].time_at_boosted>=SIX_MONTH)){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_boosted_candidates_30_to_60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
			}
		    
		    if(age_index>5 && (nodes[j].infection_status==Progression::susceptible || (nodes[j].infection_status==Progression::recovered &&nodes[j].state_before_recovery == Progression::boosted) && time_step-nodes[j].time_at_boosted>=SIX_MONTH)){//----->age between 1 and 10
				bool vacc1=false;
				//numvac1+=1;

				vacc1=true;
					if(vacc1)
					{			
					new_boosted_candidates_GT60.push_back(j);//j-th individula is vaccinated. Choose age distribution here

					}
					
					
			}

		}

//----vaccine paraemeters need work: to be updated according to new data and functional forms as and when they become available


	  auto node_update_status = update_infection(nodes[j], time_step); 
	  nodes[j].psi_T = psi_T(nodes[j], time_step);


//1000/10^6=1/1000;

//---Setting the infectiousness of new strains

	  if(time_step==GLOBAL.TIME_ALPHA){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_ALPHA);
                   if(is_new_strain){
                    nodes[j].new_strain = 1;//0,1,2,3,4
		    // nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_ALPHA;
		    nodes[j].infectiousness = nodes[j].infectiousness_original* GLOBAL.INFECTIOUSNESS_ALPHA;
			
                                  }
                              }
                    }

	  if(time_step==GLOBAL.TIME_DELTA){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_DELTA);
                   if(is_new_strain){
                    nodes[j].new_strain = 2;
		    //nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_DELTA;
		    nodes[j].infectiousness = nodes[j].infectiousness_original* GLOBAL.INFECTIOUSNESS_DELTA;

                                  }
                              }
                    }

	  if(time_step==GLOBAL.TIME_OMICRON){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_OMICRON);
                   if(is_new_strain){
                    nodes[j].new_strain = 3;
		    //nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_OMICRON;
			nodes[j].infectiousness =nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON;			
			
                                  }
                              }
                    }
	  if(time_step==GLOBAL.TIME_OMICRON_NEW){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_OMICRON_NEW);
                   if(is_new_strain){
                    nodes[j].new_strain = 4;
		    //nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_OMICRON_NEW;

		    nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_NEW;
                                  }
                              }
                    }

	  if(time_step==GLOBAL.TIME_OMICRON_BA4){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_OMICRON_BA4);
                   if(is_new_strain){
                    nodes[j].new_strain = 5;
		    //nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_OMICRON_NEW;

		    nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA4;
                                  }
                              }
                    }
	  if(time_step==GLOBAL.TIME_OMICRON_BA5){
               if((nodes[j].infection_status == Progression::exposed)||(nodes[j].infection_status == Progression::infective)||(nodes[j].infection_status == Progression::symptomatic)||(nodes[j].infection_status == Progression::hospitalised)||(nodes[j].infection_status == Progression::critical)){
                 bool is_new_strain = bernoulli(GLOBAL.FRACTION_NEW_OMICRON_BA5);
                   if(is_new_strain){
                    nodes[j].new_strain = 6;
		    //nodes[j].infectiousness *= GLOBAL.INFECTIOUSNESS_OMICRON_NEW;

		    nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA5;
                                  }
                              }
                    }					
//---New suceptible individuals by Shakir--second wave
	  if(time_step==GLOBAL.TIME_ALPHA-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_ALPHA);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }
	  //---New exposed individuals by Shakir--second (alpha) wave
	  if(time_step==GLOBAL.TIME_ALPHA){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_ALPHA);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
                  nodes[j].new_strain = 1;		//new strain 1 is alpha
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_ALPHA;

				  nodes[j].time_of_infection=time_step;
			  }
		  }
	  }

//---New susceptible individuals by Shakir--- third (delta) wave
	  if(time_step==GLOBAL.TIME_DELTA-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_DELTA);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }
	  //---New exposed individuals by Shakir--third (delta) wave
	  if(time_step==GLOBAL.TIME_DELTA){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_DELTA);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
				  nodes[j].time_of_infection=time_step;
                  nodes[j].new_strain = 2;//new_strain=2 is delta		
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_DELTA;

			  }
		  }
	  }


	  //---New susceptible individuals by Shakir--fourth (omicron) wave
	  if(time_step==GLOBAL.TIME_OMICRON-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }	  
	  if(time_step==GLOBAL.TIME_OMICRON){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_OMICRON);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
				  nodes[j].time_of_infection=time_step;
                  nodes[j].new_strain = 3;//new_strain 3 is omicron		
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON;

			  }
		  }
	  }

		
	  //---New susceptible individuals by Shakir--fifth (Omicron_new) wave
	  if(time_step==GLOBAL.TIME_OMICRON_NEW-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_NEW);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }	  
	  if(time_step==GLOBAL.TIME_OMICRON_NEW){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_OMICRON_NEW);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
				  nodes[j].time_of_infection=time_step;
                  nodes[j].new_strain = 4;//new_strain=4 is omicron new.		
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_NEW;

			  }
		  }
	  }

	  //---New susceptible individuals by Shakir--fifth (Omicron_new) wave
	  if(time_step==GLOBAL.TIME_OMICRON_BA4-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_BA4);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }	  
	  if(time_step==GLOBAL.TIME_OMICRON_BA4){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_OMICRON_BA4);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
				  nodes[j].time_of_infection=time_step;
                  nodes[j].new_strain = 5;//new_strain=4 is omicron new.		
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA4;

			  }
		  }
	  }

	  //---New susceptible individuals by Shakir--fifth (Omicron_new) wave
	  if(time_step==GLOBAL.TIME_OMICRON_BA5-60){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.FRACTION_SUSCEPTIBLE_OMICRON_BA5);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::susceptible;
			  }
		  }
	  }	  
	  if(time_step==GLOBAL.TIME_OMICRON_BA5){
		  if(nodes[j].infection_status==Progression::recovered || nodes[j].infection_status==Progression::susceptible){
			  bool is_reinfection_candidate=bernoulli(GLOBAL.REINFECTION_OMICRON_BA5);
			  if(is_reinfection_candidate){
				  nodes[j].infection_status=Progression::exposed;
				  nodes[j].time_of_infection=time_step;
                  nodes[j].new_strain = 6;//new_strain=4 is omicron new.		
		          nodes[j].infectiousness = nodes[j].infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA5;

			  }
		  }
	  }


	  if(node_update_status.new_infection){
		++num_new_infections;
		++num_total_infections;
		if (nodes[j].new_strain==0){
		++num_new_strain0_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==1){
		++num_new_strain1_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==2){
		++num_new_strain2_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==3){
		++num_new_strain3_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==4){
		++num_new_strain4_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==5){
		++num_new_strain5_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		if (nodes[j].new_strain==6){
		++num_new_strain6_infections;//this needs to be taken care of to count number of original and alpha deta omicron cases
		}
		  
		// update the mean fractions with the new data point normalized_lambda
		  
		{
		  total_lambda_fraction_data += nodes[j].lambda_incoming;
		  auto normalized_lambda = (nodes[j].lambda_incoming / nodes[j].lambda);
		  mean_lambda_fraction_data.mean_update(normalized_lambda, num_new_infections);
		  cumulative_mean_lambda_fraction_data.mean_update(normalized_lambda, num_total_infections);
		}
	  }
	  if(node_update_status.new_symptomatic){
		++num_cases;
	  }
	  if(node_update_status.new_symptomatic && nodes[j].quarantined){
		++quarantined_num_cases;
	  }
	  if(node_update_status.new_hospitalization){
		++num_cumulative_hospitalizations;
	  }
	  if(node_update_status.new_infective){
		++num_cumulative_infective;
	  }
	}
	
	//----------Selecting individuals to be vaccinated based on number of doseses delivered, This is is given by the variables 'vaccFn',"vaccFn2","vaccFn3","vaccFn4"----------------------------//


	       if ((time_step >= Time_VaccStart && new_day==0)){

				vaccinate_firstdose(nodes, new_vaccinated1_candidates_LE20, vaccFn*0.2, time_step);
				vaccinate_firstdose(nodes, new_vaccinated1_candidates_30_to_60, vaccFn*0.3, time_step);
				vaccinate_firstdose(nodes, new_vaccinated1_candidates_GT60, vaccFn*0.5, time_step);

				vaccinate_second_dose(nodes, new_vaccinated2_candidates_LE20, vaccFn2*0.2, time_step);
				vaccinate_second_dose(nodes, new_vaccinated2_candidates_30_to_60, vaccFn2*0.3, time_step);
				vaccinate_second_dose(nodes, new_vaccinated2_candidates_GT60, vaccFn2*0.5, time_step);

				vaccinate_waning_candidates(nodes, new_waning_candidates_LE20, vaccFn2*0.2, time_step);

				vaccinate_booster_dose(nodes, new_boosted_candidates_30_to_60, vaccFn3*0.3, time_step);
				vaccinate_booster_dose(nodes, new_boosted_candidates_30_to_60, vaccFn3*0.3, time_step);
				vaccinate_booster_dose(nodes, new_boosted_candidates_GT60, vaccFn3*0.5, time_step);

				vaccinate_waning2_candidates(nodes, new_waning2_candidates_LE20, vaccFn2*0.2, time_step);

				vaccinate_booster2_dose(nodes, new_boosted2_candidates_30_to_60, vaccFn4*0.3, time_step);
				vaccinate_booster2_dose(nodes, new_boosted2_candidates_30_to_60, vaccFn4*0.3, time_step);
				vaccinate_booster2_dose(nodes, new_boosted2_candidates_GT60, vaccFn4*0.5, time_step);

		   }
	
	

	update_all_kappa(nodes, homes, workplaces, communities, nbr_cells, intv_params, time_step);
	if(GLOBAL.ENABLE_TESTING){
		update_test_status(nodes, time_step);
		update_infection_testing(nodes, homes, time_step);
	    update_test_request(nodes, homes, workplaces, communities, nbr_cells, time_step,testing_protocol_file_read);
	}
	if(GLOBAL.USE_AGE_DEPENDENT_MIXING){
	  for (count_type h = 0; h < GLOBAL.num_homes; ++h){
		updated_lambda_h_age_dependent(nodes, homes[h],
									   home_age_matrix.u,
									   home_age_matrix.sigma,
									   home_age_matrix.vT);
	  }
	  for (count_type w = 0; w < GLOBAL.num_schools + GLOBAL.num_workplaces; ++w){
		if(workplaces[w].workplace_type == WorkplaceType::school){
		  updated_lambda_w_age_dependent(nodes, workplaces[w],
										 school_age_matrix.u,
										 school_age_matrix.sigma,
										 school_age_matrix.vT);
		}
		else{
		  updated_lambda_w_age_dependent(nodes, workplaces[w],
										 workplace_age_matrix.u,
										 workplace_age_matrix.sigma,
										 workplace_age_matrix.vT);
		}
		updated_lambda_project(nodes, workplaces[w]);
	  }
	}
    else{
	  for (count_type h = 0; h < GLOBAL.num_homes; ++h){
		updated_lambda_h_age_independent(nodes, homes[h]);
		//FEATURE_PROPOSAL: make the mixing dependent on node.age_group;
	  }
	  
	  for (count_type w = 0; w < GLOBAL.num_schools + GLOBAL.num_workplaces; ++w){
		updated_lambda_w_age_independent(nodes, workplaces[w]);
		updated_lambda_project(nodes, workplaces[w]);
		//FEATURE_PROPOSAL: make the mixing dependent on node.age_group;
	  }
	}

	if(GLOBAL.ENABLE_NEIGHBORHOOD_SOFT_CONTAINMENT){
	  update_grid_cell_statistics(nbr_cells, homes, nodes,
								  GLOBAL.LOCKED_NEIGHBORHOOD_LEAKAGE,
								  GLOBAL.NEIGHBORHOOD_LOCK_THRESHOLD);
	}

	for (count_type c = 0; c < GLOBAL.num_communities; ++c){
	  auto temp_stats = get_infected_community(nodes, communities[c]);
	  //let row = [time_step/SIM_STEPS_PER_DAY,c,temp_stats[0],temp_stats[1],temp_stats[2],temp_stats[3],temp_stats[4]].join(",");
	  plot_data.nums["csvContent"].push_back({time_step, {
		  c,
		  temp_stats.affected,
		  temp_stats.susceptible,
		  temp_stats.unvaccinated,
		  temp_stats.vaccinated1,
		  temp_stats.vaccinated2,
		  temp_stats.waning,
		  temp_stats.boosted,
		  temp_stats.boosted2,
		  temp_stats.exposed,
		  temp_stats.infective,
		  temp_stats.symptomatic,
		  temp_stats.hospitalised,
		  temp_stats.critical,
		  temp_stats.dead,
		  temp_stats.recovered,
		  temp_stats.infected_age_group_1,temp_stats.infected_age_group_2,temp_stats.infected_age_group_3,temp_stats.infected_age_group_4,temp_stats.infected_age_group_5,temp_stats.infected_age_group_6,temp_stats.infected_age_group_7,temp_stats.infected_age_group_8,temp_stats.infected_age_group_9,temp_stats.infected_age_group_10,\
		  temp_stats.hospitalised_age_group_1,temp_stats.hospitalised_age_group_2,temp_stats.hospitalised_age_group_3,temp_stats.hospitalised_age_group_4,temp_stats.hospitalised_age_group_5,temp_stats.hospitalised_age_group_6,temp_stats.hospitalised_age_group_7,temp_stats.hospitalised_age_group_8,temp_stats.hospitalised_age_group_9,temp_stats.hospitalised_age_group_10,\
          temp_stats.dead_age_group_1,temp_stats.dead_age_group_2,temp_stats.dead_age_group_3,temp_stats.dead_age_group_4,temp_stats.dead_age_group_5,temp_stats.dead_age_group_6,temp_stats.dead_age_group_7,temp_stats.dead_age_group_8,temp_stats.dead_age_group_9,temp_stats.dead_age_group_10,\
		  temp_stats.infected_race_group_1,temp_stats.infected_race_group_2,temp_stats.infected_race_group_3,temp_stats.infected_race_group_4,temp_stats.infected_race_group_5,temp_stats.infected_race_group_6,temp_stats.infected_race_group_7,\
		  temp_stats.hospitalised_race_group_1,temp_stats.hospitalised_race_group_2,temp_stats.hospitalised_race_group_3,temp_stats.hospitalised_race_group_4,temp_stats.hospitalised_race_group_5,temp_stats.hospitalised_race_group_6,temp_stats.hospitalised_race_group_7,\
          temp_stats.dead_race_group_1,temp_stats.dead_race_group_2,temp_stats.dead_race_group_3,temp_stats.dead_race_group_4,temp_stats.dead_race_group_5,temp_stats.dead_race_group_6,temp_stats.dead_race_group_7,\
		  temp_stats.infected_income_group_1,temp_stats.infected_income_group_2,temp_stats.infected_income_group_3,temp_stats.infected_income_group_4,temp_stats.infected_income_group_5,temp_stats.infected_income_group_6,temp_stats.infected_income_group_7,temp_stats.infected_income_group_8,\
		  temp_stats.hospitalised_income_group_1,temp_stats.hospitalised_income_group_2,temp_stats.hospitalised_income_group_3,temp_stats.hospitalised_income_group_4,temp_stats.hospitalised_income_group_5,temp_stats.hospitalised_income_group_6,temp_stats.hospitalised_income_group_7,temp_stats.hospitalised_income_group_8,\
          temp_stats.dead_income_group_1,temp_stats.dead_income_group_2,temp_stats.dead_income_group_3,temp_stats.dead_income_group_4,temp_stats.dead_income_group_5,temp_stats.dead_income_group_6,temp_stats.dead_income_group_7,temp_stats.dead_income_group_8,\
		  temp_stats.infected_ethnicity_group_1,temp_stats.infected_ethnicity_group_2,\
		  temp_stats.hospitalised_ethnicity_group_1,temp_stats.hospitalised_ethnicity_group_2,\
          temp_stats.dead_ethnicity_group_1,temp_stats.dead_ethnicity_group_2,\
		  temp_stats.infected_gender_group_1,temp_stats.infected_gender_group_2,\
		  temp_stats.hospitalised_gender_group_1,temp_stats.hospitalised_gender_group_2,\
          temp_stats.dead_gender_group_1,temp_stats.dead_gender_group_2,\		  		  
		  temp_stats.recovered_from_infective,
		  temp_stats.recovered_from_symptomatic,
		  temp_stats.recovered_from_hospitalised,
		  temp_stats.recovered_from_critical,
		  temp_stats.hd_area_affected,
		  temp_stats.hd_area_susceptible,
		  temp_stats.hd_area_exposed,
          temp_stats.hd_area_infective,
		  temp_stats.hd_area_symptomatic,
          temp_stats.hd_area_hospitalised,
          temp_stats.hd_area_critical,
          temp_stats.hd_area_dead,
		  temp_stats.hd_area_recovered,
		  temp_stats.hd_area_recovered_from_infective,
		  temp_stats.hd_area_recovered_from_symptomatic,
		  temp_stats.hd_area_recovered_from_hospitalised,
		  temp_stats.hd_area_recovered_from_critical
		  }});

	  //Update w_c value for this community, followed by update of lambdas
	  if(communities[c].individuals.size()>0){
		communities[c].w_c = interpolate(1.0, GLOBAL.LOCKED_COMMUNITY_LEAKAGE,
										 double(temp_stats.hospitalised)/double(communities[c].individuals.size()),
										 GLOBAL.COMMUNITY_LOCK_THRESHOLD);
	  } else{
		communities[c].w_c = 1;
	  }

	  updated_lambda_c_local(nodes, communities[c]);
	}

	updated_lambda_c_local_random_community(nodes, communities, homes);
	update_lambda_c_global(communities, community_fk_matrix);
	update_lambda_nbr_cells(nodes, nbr_cells, homes, communities);

	travel_fraction = updated_travel_fraction(nodes,time_step);
	travel_fraction_higher = updated_travel_fraction_higher(nodes,time_step);

	// Update lambdas for the next step
#pragma omp parallel for default(none)									\
  shared(travel_fraction, travel_fraction_higher, time_step, homes, workplaces, communities, nbr_cells, nodes, NUM_PEOPLE)
	for (count_type j = 0; j < NUM_PEOPLE; ++j){
	  update_lambdas(nodes[j], homes, workplaces, communities, nbr_cells, travel_fraction, travel_fraction_higher, time_step);
	}
	// sk ///////////////////////////
		if (GLOBAL.ENABLE_COHORTS && GLOBAL.TRAINS_RUNNING)
		{//Cohort lambda for each node is updated only here.
			update_individual_lambda_cohort(nodes, time_step, cohorts);
		}
	/////////////////////////////////
	//Get data for this simulation step
	count_type n_infected = 0,
	  n_exposed = 0,
	  n_hospitalised = 0,
	  n_symptomatic = 0,
	  n_critical = 0,
	  n_fatalities = 0,
	  n_recovered = 0,
	////////////// sk:  begin: newly added data
	n_unvaccinated=0,
	n_vaccinated1=0,
	n_vaccinated2=0,
	n_waning=0,
	n_boosted=0,
	n_boosted2=0,
	n_full_vaccinated=0,
	n_infected_age_group_1=0,
	n_infected_age_group_2=0,
	n_infected_age_group_3=0,
	n_infected_age_group_4=0,                   
	n_infected_age_group_5=0,
	n_infected_age_group_6=0,
	n_infected_age_group_7=0,
	n_infected_age_group_8=0,
	n_infected_age_group_9=0,
	n_infected_age_group_10=0,

	n_hospitalised_age_group_1=0,
	n_hospitalised_age_group_2=0,
	n_hospitalised_age_group_3=0,
	n_hospitalised_age_group_4=0,                   
	n_hospitalised_age_group_5=0,
	n_hospitalised_age_group_6=0,
	n_hospitalised_age_group_7=0,
	n_hospitalised_age_group_8=0,
	n_hospitalised_age_group_9=0,
	n_hospitalised_age_group_10=0,

	n_dead_age_group_1=0,
	n_dead_age_group_2=0,
	n_dead_age_group_3=0,
	n_dead_age_group_4=0,                   
	n_dead_age_group_5=0,
	n_dead_age_group_6=0,
	n_dead_age_group_7=0,
	n_dead_age_group_8=0,
	n_dead_age_group_9=0,
	n_dead_age_group_10=0,                                      				                                         

//----------Diversity related infections and hospitalizations: Races--------------------//
          // Races: code are 1 to 7.
	n_infected_race_group_1=0,
	n_infected_race_group_2=0,
	n_infected_race_group_3=0,
	n_infected_race_group_4=0,                   
	n_infected_race_group_5=0,
	n_infected_race_group_6=0,
	n_infected_race_group_7=0,


	n_hospitalised_race_group_1=0,
	n_hospitalised_race_group_2=0,
	n_hospitalised_race_group_3=0,
	n_hospitalised_race_group_4=0,                   
	n_hospitalised_race_group_5=0,
	n_hospitalised_race_group_6=0,
	n_hospitalised_race_group_7=0,

	n_dead_race_group_1=0,
	n_dead_race_group_2=0,
	n_dead_race_group_3=0,
	n_dead_race_group_4=0,                   
	n_dead_race_group_5=0,
	n_dead_race_group_6=0,
	n_dead_race_group_7=0,

	//----Income---
	n_infected_income_group_1=0,
	n_infected_income_group_2=0,
	n_infected_income_group_3=0,
	n_infected_income_group_4=0,                   
	n_infected_income_group_5=0,
	n_infected_income_group_6=0,
	n_infected_income_group_7=0,
	n_infected_income_group_8=0,


	n_hospitalised_income_group_1=0,
	n_hospitalised_income_group_2=0,
	n_hospitalised_income_group_3=0,
	n_hospitalised_income_group_4=0,                   
	n_hospitalised_income_group_5=0,
	n_hospitalised_income_group_6=0,
	n_hospitalised_income_group_7=0,
	n_hospitalised_income_group_8=0,

	n_dead_income_group_1=0,
	n_dead_income_group_2=0,
	n_dead_income_group_3=0,
	n_dead_income_group_4=0,                   
	n_dead_income_group_5=0,
	n_dead_income_group_6=0,
	n_dead_income_group_7=0,
	n_dead_income_group_8=0,
	
	//-----Ethnicity----
	n_infected_ethnicity_group_1=0,
	n_infected_ethnicity_group_2=0,

	n_hospitalised_ethnicity_group_1=0,
	n_hospitalised_ethnicity_group_2=0,

	n_dead_ethnicity_group_1=0,
	n_dead_ethnicity_group_2=0,

	//---------Gender
	n_infected_gender_group_1=0,
	n_infected_gender_group_2=0,

	n_hospitalised_gender_group_1=0,
	n_hospitalised_gender_group_2=0,

	n_dead_gender_group_1=0,
	n_dead_gender_group_2=0,
	////////////////////// sk: end of new... data 
//-------------------------------------------------------------------------------//
	  n_affected = 0,
	  n_infective = 0,
	  quarantined_individuals = 0,
	  quarantined_infectious = 0,
	  quarantined_individuals_cohorts=0, // sk
	  quarantined_infectious_cohorts=0; // sk

	count_type num_active_infections = 0;
	count_type num_active_infections_new_strain0 = 0,num_active_infections_new_strain1 = 0,num_active_infections_new_strain2 = 0,num_active_infections_new_strain3 = 0,num_active_infections_new_strain4 = 0,num_active_infections_new_strain5 = 0,num_active_infections_new_strain6 = 0;

	count_type n_primary_contact = 0, 
	  n_mild_symptomatic_tested = 0, //CCC2
   	  n_moderate_symptomatic_tested = 0, //DCHC
      n_severe_symptomatic_tested = 0, //DCH
      n_icu = 0,
	  n_requested_tests = 0,
	  n_tested_positive = 0;

	double susceptible_lambda = 0,
	  susceptible_lambda_H = 0,
	  susceptible_lambda_W = 0,
	  susceptible_lambda_C = 0,
	  susceptible_lambda_T = 0,
	  susceptible_lambda_PROJECT = 0,
	  susceptible_lambda_NBR_CELL = 0,
	  susceptible_lambda_RANDOM_COMMUNITY = 0;
	double curtailed_interaction = 0, normal_interaction = 0;

#pragma omp parallel for default(none) shared(nodes, GLOBAL, NUM_PEOPLE)			\
  reduction(+: n_infected, n_exposed,									\
			n_hospitalised, n_symptomatic,								\
			n_critical, n_fatalities,									\
			n_recovered, n_affected, n_infective,						\
				n_unvaccinated,n_vaccinated1,n_vaccinated2,n_waning,n_boosted,n_boosted2,           \
				n_infected_age_group_1,n_infected_age_group_2,   								\
                n_infected_age_group_3,n_infected_age_group_4,   								\
				n_infected_age_group_5,n_infected_age_group_6,   								\
				n_infected_age_group_7,n_infected_age_group_8,   								\
				n_infected_age_group_9,n_infected_age_group_10,   								\
				n_hospitalised_age_group_1,n_hospitalised_age_group_2,   								\
                n_hospitalised_age_group_3,n_hospitalised_age_group_4,   								\
				n_hospitalised_age_group_5,n_hospitalised_age_group_6,   								\
				n_hospitalised_age_group_7,n_hospitalised_age_group_8,   								\
				n_hospitalised_age_group_9,n_hospitalised_age_group_10,   								\
				n_dead_age_group_1,n_dead_age_group_2,   								\
                n_dead_age_group_3,n_dead_age_group_4,   								\
				n_dead_age_group_5,n_dead_age_group_6,   								\
				n_dead_age_group_7,n_dead_age_group_8,   								\
				n_dead_age_group_9,n_dead_age_group_10,   								\
				n_infected_race_group_1,		\
				n_infected_race_group_2,		\
				n_infected_race_group_3,		\
				n_infected_race_group_4,      \            
				n_infected_race_group_5,		\
				n_infected_race_group_6,		\
				n_infected_race_group_7,		\
				n_hospitalised_race_group_1,	\
				n_hospitalised_race_group_2,	\
				n_hospitalised_race_group_3,	\
				n_hospitalised_race_group_4,  \                
				n_hospitalised_race_group_5,	\
				n_hospitalised_race_group_6,	\
				n_hospitalised_race_group_7,	\
				n_dead_race_group_1,		\
				n_dead_race_group_2,		\
				n_dead_race_group_3,		\
				n_dead_race_group_4,      \            
				n_dead_race_group_5,		\
				n_dead_race_group_6,		\
				n_dead_race_group_7,		\
				n_infected_income_group_1,	\
				n_infected_income_group_2,	\
				n_infected_income_group_3,	\
				n_infected_income_group_4,    \              
				n_infected_income_group_5,	\
				n_infected_income_group_6,	\
				n_infected_income_group_7,	\
				n_infected_income_group_8,	\
				n_hospitalised_income_group_1,	\
				n_hospitalised_income_group_2,	\
				n_hospitalised_income_group_3,	\
				n_hospitalised_income_group_4,    \              
				n_hospitalised_income_group_5,	\
				n_hospitalised_income_group_6,	\
				n_hospitalised_income_group_7,	\
				n_hospitalised_income_group_8,	\
				n_dead_income_group_1,		\
				n_dead_income_group_2,		\
				n_dead_income_group_3,		\
				n_dead_income_group_4,        \          
				n_dead_income_group_5,		\
				n_dead_income_group_6,		\
				n_dead_income_group_7,		\
				n_dead_income_group_8,		\
				n_infected_ethnicity_group_1,	\
				n_infected_ethnicity_group_2,	\
				n_hospitalised_ethnicity_group_1,	\
				n_hospitalised_ethnicity_group_2,	\
				n_dead_ethnicity_group_1,		\
				n_dead_ethnicity_group_2,		\
				n_infected_gender_group_1,	\
				n_infected_gender_group_2,	\
				n_hospitalised_gender_group_1,	\
				n_hospitalised_gender_group_2,	\
				n_dead_gender_group_1,	\
				n_dead_gender_group_2,	\
			susceptible_lambda, susceptible_lambda_H,					\
			susceptible_lambda_W, susceptible_lambda_C,					\
			susceptible_lambda_T, susceptible_lambda_PROJECT,				\
			susceptible_lambda_NBR_CELL, susceptible_lambda_RANDOM_COMMUNITY,	 \
			quarantined_infectious, quarantined_individuals,			\
			curtailed_interaction, normal_interaction, 					\
	   		n_primary_contact, \
	  		n_mild_symptomatic_tested,  \
   	  		n_moderate_symptomatic_tested,  \
      		n_severe_symptomatic_tested, \
      		n_icu, n_requested_tests,n_tested_positive)
	for(count_type j = 0; j < NUM_PEOPLE; ++j){
	  auto infection_status = nodes[j].infection_status;
	
	  int age_index = nodes[j].age_index;			
	  int gender =nodes[j].gender;
	  int race = nodes[j].race;
	  int income =nodes[j].income;
	  int ethnicity=nodes[j].ethnicity;
	  
	  if(infection_status == Progression::susceptible){
		susceptible_lambda += nodes[j].lambda;
		susceptible_lambda_H += nodes[j].lambda_incoming.home;
		susceptible_lambda_W += nodes[j].lambda_incoming.work;
		susceptible_lambda_C += nodes[j].lambda_incoming.community;
		susceptible_lambda_T += nodes[j].lambda_incoming.travel;
		susceptible_lambda_PROJECT += nodes[j].lambda_incoming.project;
		susceptible_lambda_NBR_CELL += nodes[j].lambda_incoming.nbr_cell;
		susceptible_lambda_RANDOM_COMMUNITY += nodes[j].lambda_incoming.random_community;
	  }
	  if(infection_status == Progression::infective
		 || infection_status == Progression::symptomatic
		 || infection_status == Progression::hospitalised
		 || infection_status == Progression::critical){
		n_infected += 1;
		home_ward_infected[nodes[j].home_ward] += 1; // sk
				if(age_index==0){                   //age group infected Shakir
					n_infected_age_group_1+=1;
				}
				if(age_index==1){
					n_infected_age_group_2+=1;
				}
				if(age_index==2){
					n_infected_age_group_3+=1;
				}
				if(age_index==3){
					n_infected_age_group_4+=1;
				}
				if(age_index==4){
					n_infected_age_group_5+=1;
				}
				if(age_index==5){
					n_infected_age_group_6+=1;
				}
				if(age_index==6){
					n_infected_age_group_7+=1;
				}
				if(age_index==7){
					n_infected_age_group_8+=1;
				}
				if(age_index==8){
					n_infected_age_group_9+=1;
				}
				if(age_index==9){
					n_infected_age_group_10+=1;
				}

//------------------Race groups--------------------//
				if(race==1){                   //race group infected Shakir
					n_infected_race_group_1+=1;
				}
				if(race==2){
					n_infected_race_group_2+=1;
				}
				if(race==3){
					n_infected_race_group_3+=1;
				}
				if(race==4){
					n_infected_race_group_4+=1;
				}
				if(race==5){
					n_infected_race_group_5+=1;
				}
				if(race==6){
					n_infected_race_group_6+=1;
				}
				if(race==7){
					n_infected_race_group_7+=1;
				}

//----------------Income groups------------------//
				if(income==1){                   //income group infected Shakir
					n_infected_income_group_1+=1;
				}
				if(income==2){
					n_infected_income_group_2+=1;
				}
				if(income==3){
					n_infected_income_group_3+=1;
				}
				if(income==4){
					n_infected_income_group_4+=1;
				}
				if(income==5){
					n_infected_income_group_5+=1;
				}
				if(income==6){
					n_infected_income_group_6+=1;
				}
				if(income==7){
					n_infected_income_group_7+=1;
                }
				if(income==8){
					n_infected_income_group_8+=1;
				}
//------------------Ethnicity groups-------------------------//
				if(ethnicity==1){                   //ethnicity group infected Shakir
					n_infected_ethnicity_group_1+=1;
				}
				if(ethnicity==2){
					n_infected_ethnicity_group_2+=1;
				}
//---------------Gender groups------------------------------
				if(gender==1){                   //gender group infected Shakir
					n_infected_gender_group_1+=1;
				}
				if(gender==2){
					n_infected_gender_group_2+=1;
				}



	  }else if(infection_status != Progression::dead){
		  curtailed_interaction+=(nodes[j].kappa_H_incoming * GLOBAL.BETA_H
		  		+ nodes[j].kappa_C_incoming * GLOBAL.BETA_C
				+ ((nodes[j].workplace_type == WorkplaceType::office)?GLOBAL.BETA_W:0)*nodes[j].kappa_W_incoming
				+ ((nodes[j].workplace_type == WorkplaceType::school)?GLOBAL.BETA_S:0)*nodes[j].kappa_W_incoming
				+ ((nodes[j].workplace_type == WorkplaceType::office)?GLOBAL.BETA_PROJECT:0)*nodes[j].kappa_W_incoming
				+ ((nodes[j].workplace_type == WorkplaceType::school)?GLOBAL.BETA_CLASS:0)*nodes[j].kappa_W_incoming
				+ nodes[j].kappa_C_incoming*GLOBAL.BETA_NBR_CELLS
				+ nodes[j].kappa_C_incoming*GLOBAL.BETA_RANDOM_COMMUNITY
				+ ((nodes[j].has_to_travel)?GLOBAL.BETA_TRAVEL:0)*nodes[j].travels());
		  normal_interaction+=(GLOBAL.BETA_H
		  		+ GLOBAL.BETA_C
				+ ((nodes[j].workplace_type == WorkplaceType::office)?GLOBAL.BETA_W:0)
				+ ((nodes[j].workplace_type == WorkplaceType::school)?GLOBAL.BETA_S:0)
				+ ((nodes[j].workplace_type == WorkplaceType::office)?GLOBAL.BETA_PROJECT:0)
				+ ((nodes[j].workplace_type == WorkplaceType::school)?GLOBAL.BETA_CLASS:0)
				+ GLOBAL.BETA_NBR_CELLS
				+ GLOBAL.BETA_RANDOM_COMMUNITY
				+ ((nodes[j].has_to_travel)?GLOBAL.BETA_TRAVEL:0));
	  }
	  if(infection_status == Progression::exposed){
		n_exposed += 1;
		num_active_infections += 1;
		if (nodes[j].new_strain==0){///checking for original strain for now
		num_active_infections_new_strain0 += 1;
		}
		if (nodes[j].new_strain==1){///checking for alpha strain for now
		num_active_infections_new_strain1 += 1;
		}
		if (nodes[j].new_strain==2){///checking for delta strain for now
		num_active_infections_new_strain2 += 1;
		}
		if (nodes[j].new_strain==3){///checking for Omicron strain for now
		num_active_infections_new_strain3 += 1;
		}
		if (nodes[j].new_strain==4){///checking for OmicronNew strain for now
		num_active_infections_new_strain4 += 1;
		}
		if (nodes[j].new_strain==5){///checking for Omicron BA4 strain for now
		num_active_infections_new_strain5 += 1;
		}
		if (nodes[j].new_strain==6){///checking for omicron BA5 strain for now
		num_active_infections_new_strain6 += 1;
		}

	  }
	  if(infection_status == Progression::hospitalised){
		n_hospitalised += 1;
		num_active_infections += 1;
		if (nodes[j].new_strain==0){///checking for original strain for now
		num_active_infections_new_strain0 += 1;
		}
		if (nodes[j].new_strain==1){///checking for alpha strain for now
		num_active_infections_new_strain1 += 1;
		}
		if (nodes[j].new_strain==2){///checking for delta strain for now
		num_active_infections_new_strain2 += 1;
		}
		if (nodes[j].new_strain==3){///checking for Omicron strain for now
		num_active_infections_new_strain3 += 1;
		}
		if (nodes[j].new_strain==4){///checking for OmicronNew strain for now
		num_active_infections_new_strain4 += 1;
		}
		if (nodes[j].new_strain==5){///checking for Omicron BA4 strain for now
		num_active_infections_new_strain5 += 1;
		}
		if (nodes[j].new_strain==6){///checking for omicron BA5 strain for now
		num_active_infections_new_strain6 += 1;
		}

				if(age_index==0){                    //age group hospitalised Shakir
					n_hospitalised_age_group_1+=1;
				}
				if(age_index==1){
					n_hospitalised_age_group_2+=1;
				}
				if(age_index==2){
					n_hospitalised_age_group_3+=1;
				}
				if(age_index==3){
					n_hospitalised_age_group_4+=1;
				}
				if(age_index==4){
					n_hospitalised_age_group_5+=1;
				}
				if(age_index==5){
					n_hospitalised_age_group_6+=1;
				}
				if(age_index==6){
					n_hospitalised_age_group_7+=1;
				}
				if(age_index==7){
					n_hospitalised_age_group_8+=1;
				}
				if(age_index==8){
					n_hospitalised_age_group_9+=1;
				}
				if(age_index==9){
					n_hospitalised_age_group_10+=1;
				}

//------------------Race groups--------------------//
				if(race==1){                   //race group hospitalised Shakir
					n_hospitalised_race_group_1+=1;
				}
				if(race==2){
					n_hospitalised_race_group_2+=1;
				}
				if(race==3){
					n_hospitalised_race_group_3+=1;
				}
				if(race==4){
					n_hospitalised_race_group_4+=1;
				}
				if(race==5){
					n_hospitalised_race_group_5+=1;
				}
				if(race==6){
					n_hospitalised_race_group_6+=1;
				}
				if(race==7){
					n_hospitalised_race_group_7+=1;
				}

//----------------Income groups------------------//
				if(income==1){                   //income group hospitalised Shakir
					n_hospitalised_income_group_1+=1;
				}
				if(income==2){
					n_hospitalised_income_group_2+=1;
				}
				if(income==3){
					n_hospitalised_income_group_3+=1;
				}
				if(income==4){
					n_hospitalised_income_group_4+=1;
				}
				if(income==5){
					n_hospitalised_income_group_5+=1;
				}
				if(income==6){
					n_hospitalised_income_group_6+=1;
				}
				if(income==7){
					n_hospitalised_income_group_7+=1;
                }
				if(income==8){
					n_hospitalised_income_group_8+=1;
				}
//------------------Ethnicity groups-------------------------//
				if(ethnicity==1){                   //ethnicity group hospitalised Shakir
					n_hospitalised_ethnicity_group_1+=1;
				}
				if(ethnicity==2){
					n_hospitalised_ethnicity_group_2+=1;
				}
//---------------Gender groups------------------------------
				if(gender==1){                   //gender group hospitalised Shakir
					n_hospitalised_gender_group_1+=1;
				}
				if(gender==2){
					n_hospitalised_gender_group_2+=1;
				}


	  }
	  if(infection_status == Progression::symptomatic){
		n_symptomatic += 1;
		num_active_infections += 1;
		if (nodes[j].new_strain==0){///checking for original strain for now
		num_active_infections_new_strain0 += 1;
		}
		if (nodes[j].new_strain==1){///checking for alpha strain for now
		num_active_infections_new_strain1 += 1;
		}
		if (nodes[j].new_strain==2){///checking for delta strain for now
		num_active_infections_new_strain2 += 1;
		}
		if (nodes[j].new_strain==3){///checking for Omicron strain for now
		num_active_infections_new_strain3 += 1;
		}
		if (nodes[j].new_strain==4){///checking for OmicronNew strain for now
		num_active_infections_new_strain4 += 1;
		}
		if (nodes[j].new_strain==5){///checking for Omicron BA4 strain for now
		num_active_infections_new_strain5 += 1;
		}
		if (nodes[j].new_strain==6){///checking for omicron BA5 strain for now
		num_active_infections_new_strain6 += 1;
		}
	  }
	  if(infection_status == Progression::critical){
		n_critical += 1;
		num_active_infections += 1;

		if (nodes[j].new_strain==0){///checking for original strain for now
		num_active_infections_new_strain0 += 1;
		}
		if (nodes[j].new_strain==1){///checking for alpha strain for now
		num_active_infections_new_strain1 += 1;
		}
		if (nodes[j].new_strain==2){///checking for delta strain for now
		num_active_infections_new_strain2 += 1;
		}
		if (nodes[j].new_strain==3){///checking for Omicron strain for now
		num_active_infections_new_strain3 += 1;
		}
		if (nodes[j].new_strain==4){///checking for OmicronNew strain for now
		num_active_infections_new_strain4 += 1;
		}
		if (nodes[j].new_strain==5){///checking for Omicron BA4 strain for now
		num_active_infections_new_strain5 += 1;
		}
		if (nodes[j].new_strain==6){///checking for omicron BA5 strain for now
		num_active_infections_new_strain6 += 1;
		}
	  }
	  if(infection_status == Progression::infective){
		num_active_infections += 1;
		if (nodes[j].new_strain==0){///checking for original strain for now
		num_active_infections_new_strain0 += 1;
		}		
		if (nodes[j].new_strain==1){///checking for alpha strain for now
		num_active_infections_new_strain1 += 1;
		}
		if (nodes[j].new_strain==2){///checking for delta strain for now
		num_active_infections_new_strain2 += 1;
		}
		if (nodes[j].new_strain==3){///checking for Omicron strain for now
		num_active_infections_new_strain3 += 1;
		}
		if (nodes[j].new_strain==4){///checking for OmicronNew strain for now
		num_active_infections_new_strain4 += 1;
		}
		if (nodes[j].new_strain==5){///checking for Omicron BA4 strain for now
		num_active_infections_new_strain5 += 1;
		}
		if (nodes[j].new_strain==6){///checking for omicron BA5 strain for now
		num_active_infections_new_strain6 += 1;
		}
	  }
	  if(infection_status == Progression::dead){
		n_fatalities += 1;

				if(age_index==0){            //age group dead Shakir
					n_dead_age_group_1+=1;
				}
				if(age_index==1){
					n_dead_age_group_2+=1;
				}
				if(age_index==2){
					n_dead_age_group_3+=1;
				}
				if(age_index==3){
					n_dead_age_group_4+=1;
				}
				if(age_index==4){
					n_dead_age_group_5+=1;
				}
				if(age_index==5){
					n_dead_age_group_6+=1;
				}
				if(age_index==6){
					n_dead_age_group_7+=1;
				}
				if(age_index==7){
					n_dead_age_group_8+=1;
				}
				if(age_index==8){
					n_dead_age_group_9+=1;
				}
				if(age_index==9){
					n_dead_age_group_10+=1;
				}
//------------------Race groups--------------------//
				if(race==1){                   //race group dead Shakir
					n_dead_race_group_1+=1;
				}
				if(race==2){
					n_dead_race_group_2+=1;
				}
				if(race==3){
					n_dead_race_group_3+=1;
				}
				if(race==4){
					n_dead_race_group_4+=1;
				}
				if(race==5){
					n_dead_race_group_5+=1;
				}
				if(race==6){
					n_dead_race_group_6+=1;
				}
				if(race==7){
					n_dead_race_group_7+=1;
				}

//----------------Income groups------------------//
				if(income==1){                   //income group dead Shakir
					n_dead_income_group_1+=1;
				}
				if(income==2){
					n_dead_income_group_2+=1;
				}
				if(income==3){
					n_dead_income_group_3+=1;
				}
				if(income==4){
					n_dead_income_group_4+=1;
				}
				if(income==5){
					n_dead_income_group_5+=1;
				}
				if(income==6){
					n_dead_income_group_6+=1;
				}
				if(income==7){
					n_dead_income_group_7+=1;
                }
				if(income==8){
					n_dead_income_group_8+=1;
				}
//------------------Ethnicity groups-------------------------//
				if(ethnicity==1){                   //ethnicity group dead Shakir
					n_dead_ethnicity_group_1+=1;
				}
				if(ethnicity==2){
					n_dead_ethnicity_group_2+=1;
				}
//---------------Gender groups------------------------------
				if(gender==1){                   //gender group dead Shakir
					n_dead_gender_group_1+=1;
				}
				if(gender==2){
					n_dead_gender_group_2+=1;
				}								
	  }
	  if(infection_status == Progression::recovered){
		n_recovered += 1;
	  }
	  if((infection_status != Progression::susceptible)&&(nodes[j].state_before_recovery != Progression::vaccinated1)){
		n_affected += 1;
	  }
	  if(nodes[j].infective){
		n_infective += 1;
	  }
	  if(nodes[j].quarantined){
		quarantined_individuals += 1;
	  }
	  if(nodes[j].quarantined && (infection_status == Progression::infective
								  || infection_status == Progression::symptomatic
								  || infection_status == Progression::hospitalised
								  || infection_status == Progression::critical)){
		quarantined_infectious += 1;
	  }

	  if(nodes[j].disease_label == DiseaseLabel::primary_contact){
		  n_primary_contact +=1;
	  }
	  if(nodes[j].disease_label == DiseaseLabel::mild_symptomatic_tested){
		  n_mild_symptomatic_tested +=1;
	  }
	  if(nodes[j].disease_label == DiseaseLabel::moderate_symptomatic_tested){
		  n_moderate_symptomatic_tested +=1;
	  }
      if(nodes[j].disease_label == DiseaseLabel::severe_symptomatic_tested){
		  n_severe_symptomatic_tested +=1;
	  }
	  if(nodes[j].disease_label == DiseaseLabel::icu){
		  n_icu +=1;
	  }
	  if(nodes[j].test_status.test_requested){
		  n_requested_tests +=1;
	  }
	  if(nodes[j].test_status.tested_positive){
		  n_tested_positive +=1;
	  }
	  
//---counting vaccinated classes-----	  
	  if((nodes[j].state_before_recovery == Progression::vaccinated1 || nodes[j].vaccinated1==true) && (new_day==0) && (nodes[j].new_vaccinated1==true)){
		n_vaccinated1 += 1;
	  }

	  if((nodes[j].state_before_recovery == Progression::vaccinated2 || nodes[j].vaccinated2==true) && (new_day==0) && (nodes[j].new_vaccinated2==true)){
		n_vaccinated2 += 1;
	  }
	  if((nodes[j].state_before_recovery == Progression::boosted|| nodes[j].boosted==true) && (new_day==0) && (nodes[j].new_boosted==true)){
		n_boosted += 1;
	  }
	  if((nodes[j].state_before_recovery == Progression::waning  || nodes[j].waning==true) && (new_day==0) && (nodes[j].new_waning==true)){
		n_waning += 1;
	  }
	  if((nodes[j].state_before_recovery == Progression::boosted2 ||  nodes[j].boosted2==true) && (new_day==0) && (nodes[j].new_boosted2==true)){
		n_boosted2 += 1;
	  }  	  	  	  
	}

	//Apportion new expected infections (in next time step) to currently
	//infective nodes
	if(n_infective){
	  long double expected_infections_per_infective_node
		= (long double)(susceptible_lambda)/n_infective;
	  for(const auto& node: nodes){
		if(node.infective){
		  infections_by_new_infectives[node.time_became_infective]
			+= expected_infections_per_infective_node;
		}
	  }
	}

	plot_data.nums["num_infected"].push_back({time_step, {n_infected}});
	plot_data.nums["num_exposed"].push_back({time_step, {n_exposed}});
	plot_data.nums["num_hospitalised"].push_back({time_step, {n_hospitalised}});

	plot_data.nums["num_unvaccinated"].push_back({time_step, {n_unvaccinated}});
	plot_data.nums["num_vaccinated1"].push_back({time_step, {n_vaccinated1}});
	plot_data.nums["num_vaccinated2"].push_back({time_step, {n_vaccinated2}});
	plot_data.nums["num_waning"].push_back({time_step, {n_waning}});
	plot_data.nums["num_boosted"].push_back({time_step, {n_boosted}});
	plot_data.nums["num_boosted2"].push_back({time_step, {n_boosted2}});

	plot_data.nums["num_infected_age_group_1"].push_back({time_step,{n_infected_age_group_1}});
	plot_data.nums["num_infected_age_group_2"].push_back({time_step,{n_infected_age_group_2}});
	plot_data.nums["num_infected_age_group_3"].push_back({time_step,{n_infected_age_group_3}});
	plot_data.nums["num_infected_age_group_4"].push_back({time_step,{n_infected_age_group_4}});
	plot_data.nums["num_infected_age_group_5"].push_back({time_step,{n_infected_age_group_5}});
	plot_data.nums["num_infected_age_group_6"].push_back({time_step,{n_infected_age_group_6}});
	plot_data.nums["num_infected_age_group_7"].push_back({time_step,{n_infected_age_group_7}});
	plot_data.nums["num_infected_age_group_8"].push_back({time_step,{n_infected_age_group_8}});                        
	plot_data.nums["num_infected_age_group_9"].push_back({time_step,{n_infected_age_group_9}});
	plot_data.nums["num_infected_age_group_10"].push_back({time_step,{n_infected_age_group_10}});                                    

	plot_data.nums["num_hospitalised_age_group_1"].push_back({time_step,{n_hospitalised_age_group_1}});
	plot_data.nums["num_hospitalised_age_group_2"].push_back({time_step,{n_hospitalised_age_group_2}});
	plot_data.nums["num_hospitalised_age_group_3"].push_back({time_step,{n_hospitalised_age_group_3}});
	plot_data.nums["num_hospitalised_age_group_4"].push_back({time_step,{n_hospitalised_age_group_4}});
	plot_data.nums["num_hospitalised_age_group_5"].push_back({time_step,{n_hospitalised_age_group_5}});
	plot_data.nums["num_hospitalised_age_group_6"].push_back({time_step,{n_hospitalised_age_group_6}});
	plot_data.nums["num_hospitalised_age_group_7"].push_back({time_step,{n_hospitalised_age_group_7}});
	plot_data.nums["num_hospitalised_age_group_8"].push_back({time_step,{n_hospitalised_age_group_8}});                        
	plot_data.nums["num_hospitalised_age_group_9"].push_back({time_step,{n_hospitalised_age_group_9}});
	plot_data.nums["num_hospitalised_age_group_10"].push_back({time_step,{n_hospitalised_age_group_10}});

	plot_data.nums["num_dead_age_group_1"].push_back({time_step,{n_dead_age_group_1}});
	plot_data.nums["num_dead_age_group_2"].push_back({time_step,{n_dead_age_group_2}});
	plot_data.nums["num_dead_age_group_3"].push_back({time_step,{n_dead_age_group_3}});
	plot_data.nums["num_dead_age_group_4"].push_back({time_step,{n_dead_age_group_4}});
	plot_data.nums["num_dead_age_group_5"].push_back({time_step,{n_dead_age_group_5}});
	plot_data.nums["num_dead_age_group_6"].push_back({time_step,{n_dead_age_group_6}});
	plot_data.nums["num_dead_age_group_7"].push_back({time_step,{n_dead_age_group_7}});
	plot_data.nums["num_dead_age_group_8"].push_back({time_step,{n_dead_age_group_8}});                        
	plot_data.nums["num_dead_age_group_9"].push_back({time_step,{n_dead_age_group_9}});
	plot_data.nums["num_dead_age_group_10"].push_back({time_step,{n_dead_age_group_10}});                                    

	plot_data.nums["num_infected_race_group_1"].push_back({time_step,{n_infected_race_group_1}});
	plot_data.nums["num_infected_race_group_2"].push_back({time_step,{n_infected_race_group_2}});
	plot_data.nums["num_infected_race_group_3"].push_back({time_step,{n_infected_race_group_3}});
	plot_data.nums["num_infected_race_group_4"].push_back({time_step,{n_infected_race_group_4}});
	plot_data.nums["num_infected_race_group_5"].push_back({time_step,{n_infected_race_group_5}});
	plot_data.nums["num_infected_race_group_6"].push_back({time_step,{n_infected_race_group_6}});
	plot_data.nums["num_infected_race_group_7"].push_back({time_step,{n_infected_race_group_7}});
                             

	plot_data.nums["num_hospitalised_race_group_1"].push_back({time_step,{n_hospitalised_race_group_1}});
	plot_data.nums["num_hospitalised_race_group_2"].push_back({time_step,{n_hospitalised_race_group_2}});
	plot_data.nums["num_hospitalised_race_group_3"].push_back({time_step,{n_hospitalised_race_group_3}});
	plot_data.nums["num_hospitalised_race_group_4"].push_back({time_step,{n_hospitalised_race_group_4}});
	plot_data.nums["num_hospitalised_race_group_5"].push_back({time_step,{n_hospitalised_race_group_5}});
	plot_data.nums["num_hospitalised_race_group_6"].push_back({time_step,{n_hospitalised_race_group_6}});
	plot_data.nums["num_hospitalised_race_group_7"].push_back({time_step,{n_hospitalised_race_group_7}});


	plot_data.nums["num_dead_race_group_1"].push_back({time_step,{n_dead_race_group_1}});
	plot_data.nums["num_dead_race_group_2"].push_back({time_step,{n_dead_race_group_2}});
	plot_data.nums["num_dead_race_group_3"].push_back({time_step,{n_dead_race_group_3}});
	plot_data.nums["num_dead_race_group_4"].push_back({time_step,{n_dead_race_group_4}});
	plot_data.nums["num_dead_race_group_5"].push_back({time_step,{n_dead_race_group_5}});
	plot_data.nums["num_dead_race_group_6"].push_back({time_step,{n_dead_race_group_6}});
	plot_data.nums["num_dead_race_group_7"].push_back({time_step,{n_dead_race_group_7}});

	plot_data.nums["num_infected_income_group_1"].push_back({time_step,{n_infected_income_group_1}});
	plot_data.nums["num_infected_income_group_2"].push_back({time_step,{n_infected_income_group_2}});
	plot_data.nums["num_infected_income_group_3"].push_back({time_step,{n_infected_income_group_3}});
	plot_data.nums["num_infected_income_group_4"].push_back({time_step,{n_infected_income_group_4}});
	plot_data.nums["num_infected_income_group_5"].push_back({time_step,{n_infected_income_group_5}});
	plot_data.nums["num_infected_income_group_6"].push_back({time_step,{n_infected_income_group_6}});
	plot_data.nums["num_infected_income_group_7"].push_back({time_step,{n_infected_income_group_7}});
	plot_data.nums["num_infected_income_group_8"].push_back({time_step,{n_infected_income_group_8}});                        
                             

	plot_data.nums["num_hospitalised_income_group_1"].push_back({time_step,{n_hospitalised_income_group_1}});
	plot_data.nums["num_hospitalised_income_group_2"].push_back({time_step,{n_hospitalised_income_group_2}});
	plot_data.nums["num_hospitalised_income_group_3"].push_back({time_step,{n_hospitalised_income_group_3}});
	plot_data.nums["num_hospitalised_income_group_4"].push_back({time_step,{n_hospitalised_income_group_4}});
	plot_data.nums["num_hospitalised_income_group_5"].push_back({time_step,{n_hospitalised_income_group_5}});
	plot_data.nums["num_hospitalised_income_group_6"].push_back({time_step,{n_hospitalised_income_group_6}});
	plot_data.nums["num_hospitalised_income_group_7"].push_back({time_step,{n_hospitalised_income_group_7}});
	plot_data.nums["num_hospitalised_income_group_8"].push_back({time_step,{n_hospitalised_income_group_8}});                        


	plot_data.nums["num_dead_income_group_1"].push_back({time_step,{n_dead_income_group_1}});
	plot_data.nums["num_dead_income_group_2"].push_back({time_step,{n_dead_income_group_2}});
	plot_data.nums["num_dead_income_group_3"].push_back({time_step,{n_dead_income_group_3}});
	plot_data.nums["num_dead_income_group_4"].push_back({time_step,{n_dead_income_group_4}});
	plot_data.nums["num_dead_income_group_5"].push_back({time_step,{n_dead_income_group_5}});
	plot_data.nums["num_dead_income_group_6"].push_back({time_step,{n_dead_income_group_6}});
	plot_data.nums["num_dead_income_group_7"].push_back({time_step,{n_dead_income_group_7}});
	plot_data.nums["num_dead_income_group_8"].push_back({time_step,{n_dead_income_group_8}});                        

	plot_data.nums["num_infected_ethnicity_group_1"].push_back({time_step,{n_infected_ethnicity_group_1}});
	plot_data.nums["num_infected_ethnicity_group_2"].push_back({time_step,{n_infected_ethnicity_group_2}});
                         
	plot_data.nums["num_hospitalised_ethnicity_group_1"].push_back({time_step,{n_hospitalised_ethnicity_group_1}});
	plot_data.nums["num_hospitalised_ethnicity_group_2"].push_back({time_step,{n_hospitalised_ethnicity_group_2}});

	plot_data.nums["num_dead_ethnicity_group_1"].push_back({time_step,{n_dead_ethnicity_group_1}});
	plot_data.nums["num_dead_ethnicity_group_2"].push_back({time_step,{n_dead_ethnicity_group_2}});

	plot_data.nums["num_infected_gender_group_1"].push_back({time_step,{n_infected_gender_group_1}});
	plot_data.nums["num_infected_gender_group_2"].push_back({time_step,{n_infected_gender_group_2}});
                         
	plot_data.nums["num_hospitalised_gender_group_1"].push_back({time_step,{n_hospitalised_gender_group_1}});
	plot_data.nums["num_hospitalised_gender_group_2"].push_back({time_step,{n_hospitalised_gender_group_2}});

	plot_data.nums["num_dead_gender_group_1"].push_back({time_step,{n_dead_gender_group_1}});
	plot_data.nums["num_dead_gender_group_2"].push_back({time_step,{n_dead_gender_group_2}});


	plot_data.nums["num_symptomatic"].push_back({time_step, {n_symptomatic}});
	plot_data.nums["num_critical"].push_back({time_step, {n_critical}});
	plot_data.nums["num_fatalities"].push_back({time_step, {n_fatalities}});
	plot_data.nums["num_recovered"].push_back({time_step, {n_recovered}});
	plot_data.nums["num_affected"].push_back({time_step, {n_affected}});
	plot_data.nums["num_cases"].push_back({time_step, {num_cases}});
	plot_data.nums["num_cumulative_hospitalizations"].push_back({time_step, {num_cumulative_hospitalizations}});
	plot_data.nums["num_cumulative_infective"].push_back({time_step, {num_cumulative_infective}});

	plot_data.susceptible_lambdas["susceptible_lambda"].push_back({time_step, {susceptible_lambda}});
	plot_data.susceptible_lambdas["susceptible_lambda_H"].push_back({time_step, {susceptible_lambda_H}});
	plot_data.susceptible_lambdas["susceptible_lambda_W"].push_back({time_step, {susceptible_lambda_W}});
	plot_data.susceptible_lambdas["susceptible_lambda_C"].push_back({time_step, {susceptible_lambda_C}});
	plot_data.susceptible_lambdas["susceptible_lambda_T"].push_back({time_step, {susceptible_lambda_T}});
	plot_data.susceptible_lambdas["susceptible_lambda_PROJECT"].push_back({time_step, {susceptible_lambda_PROJECT}});
	plot_data.susceptible_lambdas["susceptible_lambda_NBR_CELL"].push_back({time_step, {susceptible_lambda_NBR_CELL}});
	plot_data.susceptible_lambdas["susceptible_lambda_RANDOM_COMMUNITY"].push_back({time_step, {susceptible_lambda_RANDOM_COMMUNITY}});

	// disease label stats
	plot_data.disease_label_stats["disease_label_stats"].push_back({time_step, {n_primary_contact,
							n_mild_symptomatic_tested, n_moderate_symptomatic_tested, 
							n_severe_symptomatic_tested, n_icu, n_requested_tests, n_tested_positive}});

	//Convert to fraction
	auto total_lambda_fraction_data_sum = total_lambda_fraction_data.sum();
	total_lambda_fraction_data /= total_lambda_fraction_data_sum;


	plot_data.total_lambda_fractions["total_fraction_lambda_H"].push_back({time_step, {total_lambda_fraction_data.home}});
	plot_data.total_lambda_fractions["total_fraction_lambda_W"].push_back({time_step, {total_lambda_fraction_data.work}});
	plot_data.total_lambda_fractions["total_fraction_lambda_C"].push_back({time_step, {total_lambda_fraction_data.community}});
	plot_data.total_lambda_fractions["total_fraction_lambda_T"].push_back({time_step, {total_lambda_fraction_data.travel}});
	plot_data.total_lambda_fractions["total_fraction_lambda_PROJECT"].push_back({time_step, {total_lambda_fraction_data.project}});
	plot_data.total_lambda_fractions["total_fraction_lambda_NBR_CELL"].push_back({time_step, {total_lambda_fraction_data.nbr_cell}});
	plot_data.total_lambda_fractions["total_fraction_lambda_RANDOM_COMMUNITY"].push_back({time_step, {total_lambda_fraction_data.random_community}});

	plot_data.mean_lambda_fractions["mean_fraction_lambda_H"].push_back({time_step, {mean_lambda_fraction_data.home}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_W"].push_back({time_step, {mean_lambda_fraction_data.work}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_C"].push_back({time_step, {mean_lambda_fraction_data.community}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_T"].push_back({time_step, {mean_lambda_fraction_data.travel}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_PROJECT"].push_back({time_step, {mean_lambda_fraction_data.project}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_NBR_CELL"].push_back({time_step, {mean_lambda_fraction_data.nbr_cell}});
	plot_data.mean_lambda_fractions["mean_fraction_lambda_RANDOM_COMMUNITY"].push_back({time_step, {mean_lambda_fraction_data.random_community}});

	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_H"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.home}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_W"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.work}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_C"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.community}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_T"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.travel}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_PROJECT"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.project}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_NBR_CELL"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.nbr_cell}});
	plot_data.cumulative_mean_lambda_fractions["cumulative_mean_fraction_lambda_RANDOM_COMMUNITY"].push_back({time_step,
																							   {cumulative_mean_lambda_fraction_data.random_community}});
	plot_data.quarantined_stats["quarantined_stats"].push_back({time_step, {
                 quarantined_individuals,
				 quarantined_infectious,
				 quarantined_num_cases
                  }});
	plot_data.curtailment_stats["curtailment_stats"].push_back({time_step, {
				 normal_interaction,
				 curtailed_interaction
                  }});
				  
	logging1 = std::to_string(time_step)+","+std::to_string(num_new_infections)+","+std::to_string(num_new_strain0_infections)+","+std::to_string(num_new_strain1_infections)+","+std::to_string(num_new_strain2_infections)+","+std::to_string(num_new_strain3_infections)+","+std::to_string(num_new_strain4_infections)+","+std::to_string(num_new_strain5_infections)+","+std::to_string(num_new_strain6_infections)+","+std::to_string(num_active_infections)+","+std::to_string(num_active_infections_new_strain0)+","+std::to_string(num_active_infections_new_strain1)+","+","+std::to_string(num_active_infections_new_strain2)+","+std::to_string(num_active_infections_new_strain3)+","+std::to_string(num_active_infections_new_strain4)+","+std::to_string(num_active_infections_new_strain5)+","+std::to_string(num_active_infections_new_strain6);
    logger1.push_back(logging1);

#ifdef DEBUG
	cerr<<std::endl<<"time_step: "<<time_step;
	auto end_time_timestep = std::chrono::high_resolution_clock::now();
  	cerr << "Time step: simulation time (ms): " << duration(start_time_timestep, end_time_timestep) << "\n";
#endif
  }
	// sk ////////////////////
	for(count_type nwards = 0; nwards < GLOBAL.num_wards; nwards++){
		//std::cout << "Ward: " << nwards <<" Infected: "<< home_ward_infected[nwards] << "\n";
		plot_data.ward_wise_stats["ward_infected"].push_back({nwards,{home_ward_infected[nwards]}});
	}
	/////////////////////////

  //Create CSV data out of the date for infections per new infective node
  plot_data.infections_by_new_infectives = {
	{"infections_by_new_infectives", {}}
  };
  for(count_type time_step = GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY; time_step < GLOBAL.START_DAY*GLOBAL.SIM_STEPS_PER_DAY+GLOBAL.NUM_TIMESTEPS; ++time_step){
	plot_data.infections_by_new_infectives["infections_by_new_infectives"].push_back({time_step,
																					  {infections_by_new_infectives[time_step]}});
  }
 plot_data.logger1=logger1;
#ifdef TIMING
  end_time = std::chrono::high_resolution_clock::now();
  cerr << "simulator: simulation time (ms): " << duration(start_time, end_time) << "\n";
#endif
  

  return plot_data;
}

