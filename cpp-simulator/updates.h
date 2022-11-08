//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#ifndef UPDATES_H_
#define UPDATES_H_

#include "models.h"
#include <vector>
#include <unordered_map> 
#include "train_loader.h"

double update_individual_lambda_h(const agent& node, int cur_time);

double update_individual_lambda_w(const agent& node, int cur_time,bool mask_wearing, double mask_scaling);

double update_individual_lambda_c(const agent& node, int cur_time,bool mask_wearing,double mask_scaling);

double update_individual_lambda_nbr_cell(const agent& node, int cur_time, bool mask_wearing,double mask_scaling);

struct node_update_status{
  bool new_infection = false;
  bool new_symptomatic = false;
  bool new_hospitalization = false;
  bool new_infective = false;
};

void update_scales_Intv(std::vector<house>& homes, std::vector<workplace>& workplaces,std::vector<community>& communities,std::vector<house>& houses, matrix<nbr_cell>& nbr_cells);

void vaccinate_firstdose(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);
void vaccinate_second_dose(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);
void vaccinate_waning_candidates(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);
void vaccinate_booster_dose(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);
void vaccinate_booster2_dose(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);
void vaccinate_waning2_candidates(std::vector<agent>& nodes, std::vector<count_type> new_vaccinated1_candidates, count_type vaccFn, count_type time_step);

//void update_vaccination_step(std::vector<house>& homes, std::vector<workplace>& workplaces,std::vector<community>& communities,std::vector<house>& houses, matrix<nbr_cell>& nbr_cells);
void new_strain_initiate(std::vector<agent>& nodes, std::vector<count_type> new_strain_candidates,int strain, count_type num_new_infections,count_type time_step);

//Returns whether the node was infected or turned symptomatic in this time step
node_update_status update_infection(agent& node, int cur_time,bool mask_wearing,double mask_scaling);

void update_all_kappa(std::vector<agent>& nodes, std::vector<house>& homes, std::vector<workplace>& workplaces, std::vector<community>& communities, matrix<nbr_cell>& nbr_cells, std::vector<intervention_params>& intv_params, int cur_time, std::vector<mask>& mask);

void updated_lambda_w_age_independent(const std::vector<agent>& nodes, workplace& workplace);

void updated_lambda_h_age_independent(const std::vector<agent>& nodes, house& home);

double updated_travel_fraction(const std::vector<agent>& nodes, int cur_time,const std::vector<mask>& mask);


std::vector<double> updated_travel_fraction_higher(const std::vector<agent>& nodes, int cur_time,const std::vector<mask>& mask);

void update_lambdas(agent&node, const std::vector<house>& homes, const std::vector<workplace>& workplaces, const std::vector<community>& communities, const std::vector<std::vector<nbr_cell>>& nbr_cells, double travel_fraction, std::vector<double> travel_fraction_higher, int cur_time,bool mask_wearing,double mask_scaling);

void updated_lambda_c_local(const std::vector<agent>& nodes, community& community);
void updated_lambda_c_local_random_community(const std::vector<agent>& nodes, const std::vector<community>& communities, std::vector<house>& houses);
void update_lambda_nbr_cells(const std::vector<agent>& nodes, std::vector<std::vector<nbr_cell>>& nbr_cells, const std::vector<house>& houses, const std::vector<community>& communities);

// sk ///////////////////////////////
//Updating the cohort lambdas
void update_cohort_edge_weights(
      std::unordered_map<count_type, std::vector<cohort_space>>& cohorts,
      const std::vector<agent>& nodes);
void update_lambda_intra_cohort(
      std::unordered_map<count_type, std::vector<cohort_space>>& cohorts, 
      std::vector<agent>& nodes, int cur_time);
void update_lambda_inter_cohort(
    const std::unordered_map<count_type, std::vector<train_coach>>& am_coachs,
    const std::unordered_map<count_type, std::vector<train_coach>>& pm_coachs,
    std::unordered_map<count_type, std::vector<cohort_space>>& cohorts,
    const TrainLoader& trains, int cur_time);
    /////////////////////////////
//Update test request and test status
void update_test_request(std::vector<agent>& nodes, const std::vector<house>& homes, const std::vector<workplace>& workplaces, const std::vector<community>& communities, std::vector<std::vector<nbr_cell>>& nbr_cells, const count_type current_time, const std::vector<testing_probability>& testing_protocol);
void update_test_status(std::vector<agent>& nodes, count_type current_time);


// Age stratification update functions.
void updated_lambda_w_age_dependent(const std::vector<agent>& nodes, workplace& workplace, const matrix<double>& workplace_tx_u, const std::vector<double>& workplace_tx_sigma, const matrix<double>& workplace_tx_vT);
void updated_lambda_project(const std::vector<agent>& nodes, workplace& workplace);

void updated_lambda_h_age_dependent(const std::vector<agent>& nodes, house& home, const matrix<double>& home_tx_u, const std::vector<double>& home_tx_sigma, const matrix<double>& home_tx_vT);

std::vector<double> updated_lambda_c_local_age_dependent(const std::vector<agent>& nodes, const community& community, const matrix<double>& community_tx_u, const std::vector<double>& community_tx_sigma, const matrix<double>& community_tx_vT);

void update_lambda_c_global(std::vector<community>& communities, const matrix<double>& community_distance_matrix);

struct casualty_stats{
  count_type affected = 0;
  count_type hd_area_affected = 0;

  count_type susceptible = 0;
  // sk /////////////////
  count_type unvaccinated = 0;//Unvaccinated counting Shakir May 20
  count_type vaccinated1 = 0;
  count_type vaccinated2 =0;
  count_type waning=0;
  count_type boosted=0;
  count_type waning2=0;
  count_type boosted2=0;

  count_type hd_area_susceptible = 0;
  count_type exposed = 0;
  count_type hd_area_exposed = 0;
  count_type infective = 0;
  count_type hd_area_infective = 0;
  count_type symptomatic = 0;
  count_type hd_area_symptomatic = 0;
  count_type hospitalised = 0;
  count_type hd_area_hospitalised = 0;
  count_type critical = 0;
  count_type hd_area_critical = 0;
  count_type dead = 0;
  count_type hd_area_dead = 0;
  count_type recovered = 0;
  count_type hd_area_recovered = 0;
  // sk //////////////////////
  count_type infected_age_group_1 = 0;
  count_type infected_age_group_2 = 0;//Age group output added by shakir
  count_type infected_age_group_3 = 0;
  count_type infected_age_group_4 = 0;
  count_type infected_age_group_5 = 0;
  count_type infected_age_group_6 = 0;
  count_type infected_age_group_7 = 0;
  count_type infected_age_group_8 = 0;
  count_type infected_age_group_9 = 0;
  count_type infected_age_group_10 = 0;


  count_type hospitalised_age_group_1 = 0;
  count_type hospitalised_age_group_2 = 0;
  count_type hospitalised_age_group_3 = 0;
  count_type hospitalised_age_group_4 = 0;
  count_type hospitalised_age_group_5 = 0;
  count_type hospitalised_age_group_6 = 0;//Age group output added by shakir
  count_type hospitalised_age_group_7 = 0;
  count_type hospitalised_age_group_8 = 0;
  count_type hospitalised_age_group_9 = 0;
  count_type hospitalised_age_group_10 = 0;

  count_type dead_age_group_1 = 0;
  count_type dead_age_group_2 = 0;
  count_type dead_age_group_3 = 0;
  count_type dead_age_group_4 = 0;//Age group output added by shakir
  count_type dead_age_group_5 = 0;
  count_type dead_age_group_6 = 0;
  count_type dead_age_group_7 = 0;
  count_type dead_age_group_8 = 0;
  count_type dead_age_group_9 = 0;
  count_type dead_age_group_10 = 0;

  count_type infected_race_group_1 = 0;
  count_type infected_race_group_2 = 0;//race group output added by shakir
  count_type infected_race_group_3 = 0;
  count_type infected_race_group_4 = 0;
  count_type infected_race_group_5 = 0;
  count_type infected_race_group_6 = 0;
  count_type infected_race_group_7 = 0;




  count_type hospitalised_race_group_1 = 0;
  count_type hospitalised_race_group_2 = 0;
  count_type hospitalised_race_group_3 = 0;
  count_type hospitalised_race_group_4 = 0;
  count_type hospitalised_race_group_5 = 0;
  count_type hospitalised_race_group_6 = 0;//race group output added by shakir
  count_type hospitalised_race_group_7 = 0;



  count_type dead_race_group_1 = 0;
  count_type dead_race_group_2 = 0;
  count_type dead_race_group_3 = 0;
  count_type dead_race_group_4 = 0;//race group output added by shakir
  count_type dead_race_group_5 = 0;
  count_type dead_race_group_6 = 0;
  count_type dead_race_group_7 = 0;


  count_type infected_income_group_1 = 0;
  count_type infected_income_group_2 = 0;//income group output added by shakir
  count_type infected_income_group_3 = 0;
  count_type infected_income_group_4 = 0;
  count_type infected_income_group_5 = 0;
  count_type infected_income_group_6 = 0;
  count_type infected_income_group_7 = 0;
  count_type infected_income_group_8 = 0;



  count_type hospitalised_income_group_1 = 0;
  count_type hospitalised_income_group_2 = 0;
  count_type hospitalised_income_group_3 = 0;
  count_type hospitalised_income_group_4 = 0;
  count_type hospitalised_income_group_5 = 0;
  count_type hospitalised_income_group_6 = 0;//income group output added by shakir
  count_type hospitalised_income_group_7 = 0;
  count_type hospitalised_income_group_8 = 0;


  count_type dead_income_group_1 = 0;
  count_type dead_income_group_2 = 0;
  count_type dead_income_group_3 = 0;
  count_type dead_income_group_4 = 0;//income group output added by shakir
  count_type dead_income_group_5 = 0;
  count_type dead_income_group_6 = 0;
  count_type dead_income_group_7 = 0;
  count_type dead_income_group_8 = 0;

  count_type infected_ethnicity_group_1 = 0;
  count_type infected_ethnicity_group_2 = 0;//ethnicity group output added by shakir

  count_type hospitalised_ethnicity_group_1 = 0;
  count_type hospitalised_ethnicity_group_2 = 0;

  count_type dead_ethnicity_group_1 = 0;
  count_type dead_ethnicity_group_2 = 0;

  count_type infected_gender_group_1 = 0;
  count_type infected_gender_group_2 = 0;//gender group output added by shakir

  count_type hospitalised_gender_group_1 = 0;
  count_type hospitalised_gender_group_2 = 0;

  count_type dead_gender_group_1 = 0;
  count_type dead_gender_group_2 = 0;
///////////////////////////////////////



  count_type recovered_from_infective = 0;
  count_type recovered_from_symptomatic = 0;
  count_type recovered_from_hospitalised = 0;
  count_type recovered_from_critical = 0;
  count_type hd_area_recovered_from_infective = 0;
  count_type hd_area_recovered_from_symptomatic = 0;
  count_type hd_area_recovered_from_hospitalised = 0;
  count_type hd_area_recovered_from_critical = 0;

};

casualty_stats get_infected_community(const std::vector<agent>& nodes, const community& community);

void update_grid_cell_statistics(matrix<nbr_cell>& nbr_cells,
								 std::vector<house>& homes,
								 std::vector<agent>& nodes,
								 double locked_neighborhood_leakage,
								 double locked_neighborhood_threshold);
void update_individual_lambda_cohort(std::vector<agent>&nodes, 
                      const int cur_time, 
                      std::unordered_map<count_type, std::vector<cohort_space>>& cohorts);

#endif
