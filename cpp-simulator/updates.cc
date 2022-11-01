//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <unordered_map>

#include "updates.h"
#include "interventions.h"
#include "testing.h"
#include "train_loader.h"

using std::cerr;
using std::vector;

namespace {
  void reset_cohort_lambdas(
    std::unordered_map<count_type, vector<cohort_space>>& cohorts) {
    for (auto& iter1 : cohorts) {
      for (auto& iter2 : iter1.second) {
        iter2.lambda_interaction_external = 0.0;
      }
    }
  }
}

bool mask_active(int cur_time){
	int mask_start_date = GLOBAL.MASK_START_DATE;
	int MASK_ON_TIME = mask_start_date * GLOBAL.SIM_STEPS_PER_DAY;
	return (cur_time >= MASK_ON_TIME && GLOBAL.MASK_ACTIVE);
}

double update_individual_lambda_h(const agent& node,int cur_time){
  return (node.infective?1.0:0.0)
	* node.kappa_T
	* node.infectiousness
	* (1 + node.severity)
	* node.kappa_H;
}

double update_individual_lambda_w(const agent& node, int cur_time){
  double mask_factor = 1.0;
  if(mask_active(cur_time) && node.compliant){
	  mask_factor = GLOBAL.MASK_FACTOR;
  }
  return (node.infective?1.0:0.0)
    * (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)
	* node.kappa_T
	* node.infectiousness
	* (1 + node.severity*(2*node.psi_T-1))
	* node.kappa_W
	* mask_factor;
}

double update_individual_lambda_c(const agent& node, int cur_time){
  double mask_factor = 1.0;
  if(mask_active(cur_time) && node.compliant){
	  mask_factor = GLOBAL.MASK_FACTOR;
  }
  return (node.infective?1.0:0.0)
	* node.kappa_T
	* node.infectiousness
	* node.funct_d_ck
	* (1 + node.severity)
	* node.kappa_C
	* mask_factor
  * node.zeta_a;
	// optimised version: return node.lambda_h * node.funct_d_ck;
}

double update_individual_lambda_nbr_cell(const agent& node, int cur_time){
  double mask_factor = 1.0;
  if(mask_active(cur_time) && node.compliant){
	  mask_factor = GLOBAL.MASK_FACTOR;
  }
  return (node.infective?1.0:0.0)
	* node.kappa_T
	* node.infectiousness
	* (1 + node.severity)
	* node.kappa_C
	* mask_factor
  * node.zeta_a;
}

//Returns whether the node was infected or turned symptomatic in this time step
node_update_status update_infection(agent& node, int cur_time){
  int age_index = node.age_index;
  bool transition = false;
  node_update_status update_status;
  //console.log(1-Math.exp(-node['lambda']/SIM_STEPS_PER_DAY))
  ///TODO: Parametrise transition times
  if (node.infection_status==Progression::susceptible){
	//#pragma omp critical
	{
	  transition = bernoulli(1-exp(-node.lambda/GLOBAL.SIM_STEPS_PER_DAY));
	}
	if(transition){
	  node.infection_status = Progression::exposed; //move to exposed state
	  node.time_of_infection = cur_time;
	  node.infective = false;
	  update_status.new_infection = true;
	  bool infected_with_new_strain = false;
	  double mPs,P_or,P_s[6]={node.lambda_higher1/node.lambda,node.lambda_higher2/node.lambda,node.lambda_higher3/node.lambda,node.lambda_higher4/node.lambda,node.lambda_higher5/node.lambda,node.lambda_higher6/node.lambda};

	  P_or=std::accumulate(P_s, P_s+6, P_or);//calculates the sum of probabilities from all variants except the original.
	  P_or=1-P_or;
	  int ds;

	  ds=std::distance(P_s, std::max_element(P_s, P_s + 6));///identifies the variant
	  ds+=1;
	  mPs=*std::max_element(P_s, P_s + 6);//---calculates the probability of dominant variant
	//   infected_with_new_strain[1-1] = bernoulli(node.lambda_higher1/node.lambda);
	//   infected_with_new_strain[2-1] = bernoulli(node.lambda_higher2/node.lambda);
	//   infected_with_new_strain[3-1] = bernoulli(node.lambda_higher3/node.lambda);
	//   infected_with_new_strain[4-1] = bernoulli(node.lambda_higher4/node.lambda);
	//   infected_with_new_strain[5-1] = bernoulli(node.lambda_higher5/node.lambda);
	//   infected_with_new_strain[6-1] = bernoulli(node.lambda_higher6/node.lambda);
	  infected_with_new_strain=bernoulli(mPs);

	  if (infected_with_new_strain && P_or<1 && ds==1){//--- This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 1;//alpha strain is 1
		  node.infectiousness = node.infectiousness_original* GLOBAL.INFECTIOUSNESS_ALPHA;
	  }
	  else if (infected_with_new_strain && P_or<1 && ds==2){//---This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 2;//delta strain is 2
		  node.infectiousness = node.infectiousness_original* GLOBAL.INFECTIOUSNESS_DELTA;
	  }
	  else if (infected_with_new_strain && P_or<1 && ds==3){//---Shakir: This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 3;//omicron is 3
		  node.infectiousness = node.infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON;
	  }
	  else if (infected_with_new_strain && P_or<1 && ds==4){//---Shakir: This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 4;//omicron_new is 4
		  node.infectiousness = node.infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_NEW;
	  }
	  else if (infected_with_new_strain && P_or<1 && ds==5){//---Shakir: This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 5;//omicron_BA4 is 5
		  node.infectiousness = node.infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA4;
	  }
	  else if (infected_with_new_strain && P_or<1 && ds==6 ){//---Shakir:This could be improved to precisely get the neighbourhood based caclucations.
		  node.new_strain = 6;//omicron_BA5 is 6
		  node.infectiousness = node.infectiousness_original*GLOBAL.INFECTIOUSNESS_OMICRON_BA5;
	  }	  
	}
  }
  else if(node.infection_status==Progression::exposed
		  && (double(cur_time) - node.time_of_infection
			  > node.incubation_period)){
	node.infection_status = Progression::infective; //move to infective state
	node.infective = true;
	node.time_became_infective = cur_time;
	update_status.new_infective = true;
  }
  else if(node.infection_status==Progression::infective
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period))){
	//#pragma omp critical
	{
	  transition = bernoulli(GLOBAL.SYMPTOMATIC_FRACTION);
	}
	if(transition){
	  node.infection_status = Progression::symptomatic; //move to symptomatic
	  node.infective = true;
	  update_status.new_symptomatic = true;
	  node.entered_symptomatic_state = true;
	}
	else {
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::recovered; //move to recovered
	  node.infective = false;
	}
  }
    else if(node.infection_status==Progression::recovered
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.recovered_to_sususceptible_period))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//The recovered fraction from the asymptomatic stage have lost imminity after 

  else if(node.infection_status==Progression::symptomatic
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period))){
	//#pragma omp critical
	{
	if(node.new_strain==1){//--->Changed fron new_strain to alpha delta etc by shakir
		double probability_alpha_strain=pow(GLOBAL.VIRULENT_NEW_ALPHA,0.33)*STATE_TRAN1[age_index][0]*STATE_TRAN_CoMorb1[node.comorbidity][0];
	  transition = bernoulli(probability_alpha_strain);
	}
	else if(node.new_strain==2){
		double probability_delta_strain=pow(GLOBAL.VIRULENT_NEW_DELTA,0.33)*STATE_TRAN2[age_index][0]*STATE_TRAN_CoMorb2[node.comorbidity][0];
	  transition = bernoulli(probability_delta_strain);
	}	
	else if(node.new_strain==3){
		double probability_omicron_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON,0.33)*STATE_TRAN3[age_index][0]*STATE_TRAN_CoMorb3[node.comorbidity][0];
	  transition = bernoulli(probability_omicron_strain);
	}		
	else if(node.new_strain==4){
		double probability_omicronNew_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_NEW,0.33)*STATE_TRAN4[age_index][0]*STATE_TRAN_CoMorb4[node.comorbidity][0];
	  transition = bernoulli(probability_omicronNew_strain);
	}			
	else if(node.new_strain==5){
		double probability_omicronBA4_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA4,0.33)*STATE_TRAN5[age_index][0]*STATE_TRAN_CoMorb5[node.comorbidity][0];
	  transition = bernoulli(probability_omicronBA4_strain);
	}
	else if(node.new_strain==6){
		double probability_omicronBA5_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA5,0.33)*STATE_TRAN6[age_index][0]*STATE_TRAN_CoMorb6[node.comorbidity][0];
	  transition = bernoulli(probability_omicronBA5_strain);
	}					
	else {
		transition = bernoulli(STATE_TRAN[age_index][0]*STATE_TRAN_CoMorb[node.comorbidity][0]);//----Shakir:Comorbidity enhances chances of hospitlaizations
	}
	}
	if(transition){
	  node.infection_status = Progression::hospitalised; //move to hospitalisation
	  node.infective = false;
	  update_status.new_hospitalization = true;
	  node.entered_hospitalised_state = true;
	}
	else {
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::recovered; //move to recovered
	  node.infective = false;
	}
  }

  else if(node.infection_status==Progression::recovered
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period
				 + node.recovered_to_sususceptible_period))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//the individulas that recover from symptomatic infection could again loose immunity and become susceptible again.

  else if(node.infection_status==Progression::hospitalised
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period
				 + node.hospital_regular_period))){
	//#pragma omp critical
	{
	if(node.new_strain==1){//--->Changed fron new_strain to alpha delta etc by shakir
		double probability_alpha_strain=pow(GLOBAL.VIRULENT_NEW_ALPHA,0.33)*STATE_TRAN1[age_index][1]*STATE_TRAN_CoMorb1[node.comorbidity][1];
	  transition = bernoulli(probability_alpha_strain);
	}
	else if(node.new_strain==2){
		double probability_delta_strain=pow(GLOBAL.VIRULENT_NEW_DELTA,0.33)*STATE_TRAN2[age_index][1]*STATE_TRAN_CoMorb2[node.comorbidity][1];
	  transition = bernoulli(probability_delta_strain);
	}	
	else if(node.new_strain==3){
		double probability_omicron_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON,0.33)*STATE_TRAN3[age_index][1]*STATE_TRAN_CoMorb3[node.comorbidity][1];
	  transition = bernoulli(probability_omicron_strain);
	}		
	else if(node.new_strain==4){
		double probability_omicronNew_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_NEW,0.33)*STATE_TRAN4[age_index][1]*STATE_TRAN_CoMorb4[node.comorbidity][1];
	  transition = bernoulli(probability_omicronNew_strain);
	}
	else if(node.new_strain==5){
		double probability_omicronBA4_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA4,0.33)*STATE_TRAN5[age_index][1]*STATE_TRAN_CoMorb5[node.comorbidity][1];
	  transition = bernoulli(probability_omicronBA4_strain);
	}
	else if(node.new_strain==6){
		double probability_omicronBA5_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA5,0.33)*STATE_TRAN6[age_index][1]*STATE_TRAN_CoMorb6[node.comorbidity][1];
	  transition = bernoulli(probability_omicronBA5_strain);
	}
	else{
		transition=bernoulli(STATE_TRAN[age_index][1]*STATE_TRAN_CoMorb[node.comorbidity][1]);//---Shakir:comorbidity enhances chances for critical care
	}
	}
	if(transition){
	  node.infection_status = Progression::critical; //move to critical care
	  node.infective = false;
	}
	else {
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::recovered; //move to recovered
	  node.infective = false;
	}
  }//the fraction of hospitalized individuals can either become critical or recover after hospitalization.

  else if(node.infection_status==Progression::recovered
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period
				 + node.hospital_regular_period
				 + node.recovered_to_sususceptible_period))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
  }
  }//recovered individulas after hospitalization can become susceptible again

  else if(node.infection_status==Progression::critical
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period
				 + node.hospital_regular_period
				 + node.hospital_critical_period))){
	//#pragma omp critical
	{
	  if(node.new_strain==1){//--->new_strain changed to alpha delta etc by shakir
		double probability_alpha_strain=pow(GLOBAL.VIRULENT_NEW_ALPHA,0.33)*STATE_TRAN1[age_index][2]*STATE_TRAN_CoMorb1[node.comorbidity][2];
	  transition = bernoulli(probability_alpha_strain);
	}
	  else if(node.new_strain==2){//--->new_strain changed to alpha delta etc by shakir
		double probability_delta_strain=pow(GLOBAL.VIRULENT_NEW_DELTA,0.33)*STATE_TRAN2[age_index][2]*STATE_TRAN_CoMorb2[node.comorbidity][2];
	  transition = bernoulli(probability_delta_strain);
	}
	  else if(node.new_strain==3){//--->new_strain changed to alpha delta etc by shakir
		double probability_omicron_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON,0.33)*STATE_TRAN3[age_index][2]*STATE_TRAN_CoMorb3[node.comorbidity][2];
	  transition = bernoulli(probability_omicron_strain);
	}	
	  else if(node.new_strain==4){//--->new_strain changed to alpha delta etc by shakir
		double probability_omicronNew_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_NEW,0.33)*STATE_TRAN4[age_index][2]*STATE_TRAN_CoMorb4[node.comorbidity][2];
	  transition = bernoulli(probability_omicronNew_strain);
	}		
	else if(node.new_strain==5){
		double probability_omicronBA4_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA4,0.33)*STATE_TRAN5[age_index][2]*STATE_TRAN_CoMorb5[node.comorbidity][2];
	  transition = bernoulli(probability_omicronBA4_strain);
	}
	else if(node.new_strain==6){
		double probability_omicronBA5_strain=pow(GLOBAL.VIRULENT_NEW_OMICRON_BA5,0.33)*STATE_TRAN6[age_index][2]*STATE_TRAN_CoMorb6[node.comorbidity][2];
	  transition = bernoulli(probability_omicronBA5_strain);
	}
	else {
		transition = bernoulli(STATE_TRAN[age_index][2]*STATE_TRAN_CoMorb[node.comorbidity][2]);//----Shakir:Comorbidity enhances chances of death
	}
	}
	if(transition){
	  node.infection_status = Progression::dead;//move to dead
	  node.infective = false;
	}
	else {
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::recovered;//move to recovered
	  node.infective = false;
	}
  }

  //# This part is added to make recovered individual susceptible again
//#Date of adding this step is April 29, 2022. Code edited by Shakir.
    else if(node.infection_status==Progression::recovered
		  && (double(cur_time) - node.time_of_infection
			  > (node.incubation_period
				 + node.asymptomatic_period
				 + node.symptomatic_period
				 + node.hospital_regular_period
				 + node.hospital_critical_period+node.recovered_to_sususceptible_period))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//recovered individuals from critical care can also loose immunity.
  
  //              nodes[new_vaccinated1_candidates_30_to_60[a]].infection_status = Progression::recovered;
//#Date of adding this step is April 29, 2022. Code edited by Shakir.
    else if((node.infection_status==Progression::recovered) && (node.state_before_recovery==Progression::vaccinated1)
		  && (double(cur_time)-double(node.time_at_vaccine1)
			  > (node.recovered_to_sususceptible_period_vaccinated1))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//immune individuals from vaccinated1 dose can also loose immunity.
    else if(node.infection_status==Progression::recovered && (node.state_before_recovery==Progression::vaccinated2)
		  && (double(cur_time)-double(node.time_at_vaccine2)
			  > (node.recovered_to_sususceptible_period_vaccinated2))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//immune individuals from vaccinated2 dose can also loose immunity.

     else if(node.infection_status==Progression::recovered && (node.state_before_recovery==Progression::waning)
		  && (double(cur_time)-double(node.time_at_vaccine2)
			  -5*30.0*double(GLOBAL.SIM_STEPS_PER_DAY)> (node.recovered_to_sususceptible_period_waning))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//immune individuals from waning dose can also loose immunity.

    else if(node.infection_status==Progression::recovered && (node.state_before_recovery==Progression::boosted)
		  && (double(cur_time)-double(node.time_at_boosted)
			  > (node.recovered_to_sususceptible_period_boosted))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//immune individuals from boosted1 dose can also loose immunity.


    else if(node.infection_status==Progression::recovered && (node.state_before_recovery==Progression::boosted2)
		  && (double(cur_time)-double(node.time_at_boosted2)
			  > (node.recovered_to_sususceptible_period_boosted2))){
	//#pragma omp critical
	//{
	  //transition = bernoulli(STATE_TRAN[age_index][2]);//The loss of immunity is age dependent, as given by STATE_TRAN matrix as defined in models.h. second column is for symptomatic to hostpilaised.
	//}
	//if(transition)
	{
	  node.state_before_recovery = node.infection_status;
	  node.infection_status = Progression::susceptible;//loose immunity and become susceptible again
	  node.infective = false;
	}
  }//immune individuals from boosted2 dose can also loose immunity.


  node.lambda_h = update_individual_lambda_h(node,cur_time);
  node.lambda_w = update_individual_lambda_w(node,cur_time);
  node.lambda_c = update_individual_lambda_c(node,cur_time);
  node.lambda_nbr_cell = update_individual_lambda_nbr_cell(node,cur_time);

  return update_status;
}

void update_all_kappa(vector<agent>& nodes, vector<house>& homes, vector<workplace>& workplaces, vector<community>& communities, matrix<nbr_cell>& nbr_cells, vector<intervention_params>& intv_params, int cur_time){
  intervention_params intv_params_local;
  if(cur_time < GLOBAL.NUM_DAYS_BEFORE_INTERVENTIONS*GLOBAL.SIM_STEPS_PER_DAY){
    //get_kappa_no_intervention(nodes, homes, workplaces, communities,cur_time);
    get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
  }
  else{
    switch(GLOBAL.INTERVENTION){
    case Intervention::no_intervention:
      //get_kappa_no_intervention(nodes, homes, workplaces, communities, cur_time);
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::case_isolation:
      intv_params_local.case_isolation = true;
      //get_kappa_case_isolation(nodes, homes, workplaces, communities, cur_time);
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::home_quarantine:
      //get_kappa_home_quarantine(nodes, homes, workplaces, communities, cur_time);
      intv_params_local.home_quarantine = true;
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::lockdown:
      //get_kappa_lockdown(nodes, homes, workplaces, communities, cur_time);
      intv_params_local.lockdown = true;
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::case_isolation_and_home_quarantine:
      //get_kappa_CI_HQ(nodes, homes, workplaces, communities, cur_time);
      intv_params_local.case_isolation = true;
      intv_params_local.home_quarantine = true;
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::case_isolation_and_home_quarantine_sd_65_plus:
      //get_kappa_CI_HQ_65P(nodes, homes, workplaces, communities, cur_time);
      intv_params_local.case_isolation = true;
      intv_params_local.home_quarantine = true;
      intv_params_local.social_dist_elderly = true;
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    case Intervention::lockdown_fper_ci_hq_sd_65_plus_sper_ci:
      get_kappa_LOCKDOWN_fper_CI_HQ_SD_65_PLUS_sper_CI(nodes, homes, workplaces, communities, cur_time,
                                                       GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD);
      break;
    case Intervention::lockdown_fper:
      get_kappa_LOCKDOWN_fper(nodes, homes, workplaces, communities, cur_time, GLOBAL.FIRST_PERIOD);
      break;
    case Intervention::ld_fper_ci_hq_sd65_sc_sper_sc_tper:
      get_kappa_LD_fper_CI_HQ_SD65_SC_sper_SC_tper(nodes, homes, workplaces, communities, cur_time,
                                                   GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD, GLOBAL.THIRD_PERIOD);
      break;
    case Intervention::ld_fper_ci_hq_sd65_sc_sper:
      get_kappa_LD_fper_CI_HQ_SD65_SC_sper(nodes, homes, workplaces, communities, cur_time,
                                           GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD);
      break;
    case Intervention::ld_fper_ci_hq_sd65_sc_oe_sper:
      get_kappa_LD_fper_CI_HQ_SD65_SC_OE_sper(nodes, homes, workplaces, communities, cur_time,
                                              GLOBAL.FIRST_PERIOD, GLOBAL.OE_SECOND_PERIOD);
	  break;
    case Intervention::intv_fper_intv_sper_intv_tper:
      get_kappa_intv_fper_intv_sper_intv_tper(nodes, homes, workplaces, communities, cur_time,
                                                   GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD, GLOBAL.THIRD_PERIOD);
      break;
    case Intervention::intv_NYC:
      get_kappa_NYC(nodes, homes, workplaces, communities, cur_time);
      break;
    case Intervention::intv_Mum:
	  get_kappa_Mumbai_alternative_version(nodes, homes, workplaces, communities, nbr_cells, cur_time,
                                                   GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD);
      break;
	case Intervention::intv_Mum_cyclic:
      get_kappa_Mumbai_cyclic(nodes, homes, workplaces, communities, nbr_cells, cur_time,
							  GLOBAL.FIRST_PERIOD, GLOBAL.SECOND_PERIOD);
	  break;
	case Intervention::intv_Hillsborough:
      get_kappa_Hillsborough(nodes, homes, workplaces, communities, cur_time);
	  break;
    case Intervention::intv_nbr_containment:
      get_kappa_containment(nodes, homes, workplaces, communities, nbr_cells, cur_time, GLOBAL.FIRST_PERIOD, Intervention::intv_nbr_containment);
      break;
    case Intervention::intv_ward_containment:
      get_kappa_containment(nodes, homes, workplaces, communities, nbr_cells, cur_time, GLOBAL.FIRST_PERIOD, Intervention::intv_ward_containment);
      break;
    case Intervention::intv_file_read:
      get_kappa_file_read(nodes, homes, workplaces, communities, nbr_cells, intv_params, cur_time);
      break;
    default:
      //get_kappa_no_intervention(nodes, homes, workplaces, communities, cur_time);
      get_kappa_custom_modular(nodes, homes, workplaces, communities, nbr_cells, cur_time, intv_params_local);
      break;
    }
  }
}


void updated_lambda_project(const vector<agent>& nodes, workplace& workplace){
  for(count_type i=0; i < workplace.projects.size(); ++i){
	  double sum_value_project = 0;
//	  double sum_value_project_higher=0,sum_alpha=0,sum_delta = 0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;
//-------lambda_incoming_higher for variants begins--shakir------//
	  double sum_value_project_higher[6]={0};

//------lambda_incoming_higher for variants ends--shakir------//  
	  
	  for(count_type j=0; j < workplace.projects[i].individuals.size(); ++j){
		  sum_value_project += nodes[workplace.projects[i].individuals[j]].lambda_w;
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==1){
			  sum_value_project_higher[1-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==2){
			  sum_value_project_higher[2-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }		  
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==3){
			  sum_value_project_higher[3-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }		  
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==4){
			  sum_value_project_higher[4-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==5){
			  sum_value_project_higher[5-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }
		  if (nodes[workplace.projects[i].individuals[j]].new_strain==6){
			  sum_value_project_higher[6-1] += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		  }
		//   if (nodes[workplace.projects[i].individuals[j]].new_strain){
		// 	  sum_value_project_higher += nodes[workplace.projects[i].individuals[j]].lambda_w;  
		//   }	//----originally oly this condition was there	  		  		  		  		  
	  }

	  workplace.projects[i].age_independent_mixing = workplace.projects[i].scale*sum_value_project;
//	  workplace.projects[i].age_independent_mixing_higher = workplace.projects[i].scale*sum_value_project_higher;//---delete later shakir

//-------lambda_incoming_higher for variants begins--shakir------//
	  workplace.projects[i].age_independent_mixing_higher1 = workplace.projects[i].scale*sum_value_project_higher[1-1];
	  workplace.projects[i].age_independent_mixing_higher2 = workplace.projects[i].scale*sum_value_project_higher[2-1];
	  workplace.projects[i].age_independent_mixing_higher3 = workplace.projects[i].scale*sum_value_project_higher[3-1];
	  workplace.projects[i].age_independent_mixing_higher4 = workplace.projects[i].scale*sum_value_project_higher[4-1];
	  workplace.projects[i].age_independent_mixing_higher5 = workplace.projects[i].scale*sum_value_project_higher[5-1];
	  workplace.projects[i].age_independent_mixing_higher6 = workplace.projects[i].scale*sum_value_project_higher[6-1];


//------lambda_incoming_higher for variants ends--shakir------//  

            }
          }
double get_individual_lambda_cohort(const agent& node, int cur_time){
  double mask_factor = 1.0;
  if(mask_active(cur_time) && node.compliant){
	  mask_factor = GLOBAL.MASK_FACTOR;
  }
  if (!GLOBAL.TRAINS_RUNNING) {
	  return 0;
  }
  return (node.infective?1.0:0.0)
    * node.my_cohort.edge_weight //TODO[NKV]: We would need to update this edge weight from cohorts.cc, I guess!
	* node.kappa_T
	* node.infectiousness
	* mask_factor
	* node.kappa_W //Nodes contribution is weighted by kappa_w
	* node.attending;
	// /	* (1 + node.severity)
}

void update_lambda_intra_cohort(std::unordered_map<count_type, vector<cohort_space>>& cohorts, vector<agent>& nodes, int cur_time){

	for (auto& it1: cohorts) { //cohort_hash (src*100) + dst
		for (auto& cohort_it: it1.second) {
			double sum_value = 0;
			if(cur_time % 2){
				for (auto& j: cohort_it.internal_nodes){
					sum_value += get_individual_lambda_cohort(nodes[j], cur_time);
				}
			}
			cohort_it.lambda_interaction_internal = cohort_it.scale * sum_value;
			// if (bernoulli(0.01)){std::cout<<cohort_it.lambda_interaction_internal<<cohort_it.internal_nodes.size()<<std::endl;}
		}
	}
}

void update_cohort_edge_weights(std::unordered_map<count_type, vector<cohort_space>>& cohorts, const vector<agent>& nodes){
	//TODO[v2]: update this logic for inter cohort interactions - overlap time
	return;
}

//TODO[v2]: Update the code to work on the unordered_map
void update_lambda_inter_cohort(
    const std::unordered_map<count_type, std::vector<train_coach>>& am_coachs,
    const std::unordered_map<count_type, std::vector<train_coach>>& pm_coachs,
    std::unordered_map<count_type, vector<cohort_space>>& cohorts,const TrainLoader& trains, int cur_time){
  reset_cohort_lambdas(cohorts);
  //morning trip
  if (cur_time % GLOBAL.SIM_STEPS_PER_DAY == 1){
    for (auto& coach_line : am_coachs) {
      for (auto& coach : coach_line.second) {
        for (auto& cohort_ptr: coach.cohorts) {
          // auto current_cohort = cohorts[unordered_index][vector_index];
          double sum_value = 0.0;

          for (auto& cohort_ptr2: coach.cohorts) {
            if (cohort_ptr == cohort_ptr2) {
              continue;
            }
            // auto other_cohort = cohorts[unordered_index2][vector_index2];
            double overlap_time = trains.GetOverlapMinutesAlongLine(
              coach.trainLine,
              cohort_ptr->source_station,
              cohort_ptr->destination_station,
              cohort_ptr2->source_station,
              cohort_ptr2->destination_station);
            sum_value += cohort_ptr2->lambda_interaction_internal * overlap_time;
          }
          cohort_ptr->lambda_interaction_external += sum_value;
        }
      }
    }
  }
  //evening trip
  else if (cur_time % GLOBAL.SIM_STEPS_PER_DAY == 3){
    for (auto& coach_line : pm_coachs) {
      for (auto& coach : coach_line.second) {
        for (auto& cohort_ptr: coach.cohorts) {
          double sum_value = 0.0;

          for (auto& cohort_ptr2: coach.cohorts) {
            if (cohort_ptr == cohort_ptr2) {
              continue;
            }
			// Reversing src dest for PM.
            double overlap_time = trains.GetOverlapMinutesAlongLine(
              coach.trainLine,
              cohort_ptr->destination_station,
              cohort_ptr->source_station,
              cohort_ptr2->destination_station,
              cohort_ptr2->source_station);
            sum_value += cohort_ptr2->lambda_interaction_internal * overlap_time;
          }
          cohort_ptr->lambda_interaction_external += sum_value;
        }
      }
    }
  }
}

void updated_lambda_w_age_independent(const vector<agent>& nodes, workplace& workplace){
  double sum_value = 0;
//  double sum_value_higher = 0,sum_alpha=0,sum_delta=0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;
//-------lambda_incoming_higher for variants begins--shakir------//
double sum_value_higher[6] = {0};  

//------lambda_incoming_higher for variants ends--shakir------//  

  vector<double> lambda_age_group(GLOBAL.NUM_AGE_GROUPS);
  for (count_type i=0; i < workplace.individuals.size(); ++i){
	sum_value += nodes[workplace.individuals[i]].lambda_w;
		  if (nodes[workplace.individuals[i]].new_strain==1){
			  sum_value_higher[1-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }
		  if (nodes[workplace.individuals[i]].new_strain==2){
			  sum_value_higher[2-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }		  
		  if (nodes[workplace.individuals[i]].new_strain==3){
			  sum_value_higher[3-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }		  
		  if (nodes[workplace.individuals[i]].new_strain==4){
			  sum_value_higher[4-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }	
		  if (nodes[workplace.individuals[i]].new_strain==5){
			  sum_value_higher[5-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }
		  if (nodes[workplace.individuals[i]].new_strain==6){
			  sum_value_higher[6-1] += nodes[workplace.individuals[i]].lambda_w;  
		  }	
	// if (nodes[workplace.individuals[i]].new_strain){
	// 	sum_value_higher += nodes[workplace.individuals[i]].lambda_w;//-----original new_strain blcoked by shakir
	// }
  }


  workplace.age_independent_mixing = workplace.scale*sum_value; 
//  workplace.age_independent_mixing_higher = workplace.scale*sum_value_higher; 
//-------lambda_incoming_higher for variants begins--shakir------//
  workplace.age_independent_mixing_higher1 = workplace.scale*sum_value_higher[1-1]; 
  workplace.age_independent_mixing_higher2 = workplace.scale*sum_value_higher[2-1]; 
  workplace.age_independent_mixing_higher3 = workplace.scale*sum_value_higher[3-1]; 
  workplace.age_independent_mixing_higher4 = workplace.scale*sum_value_higher[4-1]; 
  workplace.age_independent_mixing_higher5 = workplace.scale*sum_value_higher[5-1]; 
  workplace.age_independent_mixing_higher6 = workplace.scale*sum_value_higher[6-1]; 


//------lambda_incoming_higher for variants ends--shakir------//  

}

void updated_lambda_h_age_independent(const vector<agent>& nodes,  house& home){
  double sum_value = 0;

//  double sum_value_higher = 0,sum_alpha=0,sum_delta=0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;
//-------lambda_incoming_higher for variants begins--shakir------//

double sum_value_higher[6] = {0};
//------lambda_incoming_higher for variants ends--shakir------//  

  for (count_type i=0; i<home.individuals.size(); ++i){
	sum_value += nodes[home.individuals[i]].lambda_h;

		  if (nodes[home.individuals[i]].new_strain==1){
			  sum_value_higher[1-1] += nodes[home.individuals[i]].lambda_h;  
		  }
		  if (nodes[home.individuals[i]].new_strain==2){
			  sum_value_higher[2-1] += nodes[home.individuals[i]].lambda_h;  
		  }		  
		  if (nodes[home.individuals[i]].new_strain==3){
			  sum_value_higher[3-1] += nodes[home.individuals[i]].lambda_h;  
		  }		  
		  if (nodes[home.individuals[i]].new_strain==4){
			  sum_value_higher[4-1] += nodes[home.individuals[i]].lambda_h;  
		  }	
		  if (nodes[home.individuals[i]].new_strain==5){
			  sum_value_higher[5-1] += nodes[home.individuals[i]].lambda_h;  
		  }
		  if (nodes[home.individuals[i]].new_strain==6){
			  sum_value_higher[6-1] += nodes[home.individuals[i]].lambda_h;  
		  }	
	// if (nodes[home.individuals[i]].new_strain){
	// 	sum_value_higher += nodes[home.individuals[i]].lambda_h;
	// }///----original new_strain if condition blocked by shakir
  }




  home.age_independent_mixing =  home.scale*sum_value;
//  home.age_independent_mixing_higher =  home.scale*sum_value_higher;//----delete shakir later
 //-------lambda_incoming_higher for variants begins--shakir------//
  home.age_independent_mixing_higher1 =  home.scale*sum_value_higher[1-1];
  home.age_independent_mixing_higher2 =  home.scale*sum_value_higher[2-1];
  home.age_independent_mixing_higher3 =  home.scale*sum_value_higher[3-1];
  home.age_independent_mixing_higher4 =  home.scale*sum_value_higher[4-1];
  home.age_independent_mixing_higher5 =  home.scale*sum_value_higher[5-1];
  home.age_independent_mixing_higher6 =  home.scale*sum_value_higher[6-1];


//------lambda_incoming_higher for variants ends--shakir------//  
 
}

void updated_lambda_h_age_dependent(const vector<agent>& nodes,  house& home, const matrix<double>& home_tx_u, const vector<double>& home_tx_sigma, const matrix<double>& home_tx_vT){
  auto size = home_tx_u.size();

  vector<double> age_component(GLOBAL.NUM_AGE_GROUPS, 0.0);
  vector<double> lambda_age_group(GLOBAL.NUM_AGE_GROUPS, 0.0);
  vector<double> V_tx(GLOBAL.SIGNIFICANT_EIGEN_VALUES, 0.0);

  for (count_type i=0; i<home.individuals.size(); ++i){
      int ind_age_group = nodes[home.individuals[i]].age_group;
      age_component[ind_age_group] += nodes[home.individuals[i]].lambda_h;
  }

  for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
    for(count_type count=0; count<size; ++count){
      V_tx[eigen_count] += home_tx_vT[eigen_count][count]
                           * age_component[count];
    }
  }

  for (count_type count=0; count<GLOBAL.NUM_AGE_GROUPS; ++count){
    for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
      lambda_age_group[count] += home_tx_u[count][eigen_count]
                          * home_tx_sigma[eigen_count]
                          * V_tx[eigen_count];
    }
	lambda_age_group[count] *= home.scale;
  }
  home.age_dependent_mixing = std::move(lambda_age_group);
}

void updated_lambda_w_age_dependent(const vector<agent>& nodes, workplace& workplace, const matrix<double>& workplace_tx_u, const vector<double>& workplace_tx_sigma, const matrix<double>& workplace_tx_vT){

    auto size = workplace_tx_u.size();

    vector<double> age_component(GLOBAL.NUM_AGE_GROUPS, 0.0);
    vector<double> lambda_age_group(GLOBAL.NUM_AGE_GROUPS, 0.0);
    vector<double> V_tx(GLOBAL.SIGNIFICANT_EIGEN_VALUES, 0.0);
    for (count_type i=0; i<workplace.individuals.size(); ++i){
        int ind_age_group = nodes[workplace.individuals[i]].age_group;
        age_component[ind_age_group] += nodes[workplace.individuals[i]].lambda_w;
    }

    for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
      for(count_type count=0; count<size; ++count){
        V_tx[eigen_count] += workplace_tx_vT[eigen_count][count]
                             * age_component[count];
      }
    }

    for (count_type count=0; count<GLOBAL.NUM_AGE_GROUPS; ++count){
      for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
        lambda_age_group[count] += workplace_tx_u[count][eigen_count]
                            * workplace_tx_sigma[eigen_count]
                            * V_tx[eigen_count];
      }
	  lambda_age_group[count] *=  workplace.scale;
    }
    workplace.age_dependent_mixing = std::move(lambda_age_group);
}

vector<double> updated_lambda_c_local_age_dependent(const vector<agent>& nodes, const community& community, const matrix<double>& community_tx_u, const vector<double>& community_tx_sigma, const matrix<double>& community_tx_vT){

  auto size = community_tx_u.size();

  vector<double> age_component(GLOBAL.NUM_AGE_GROUPS, 0.0);
  vector<double> lambda_age_group(GLOBAL.NUM_AGE_GROUPS, 0.0);
  vector<double> V_tx(GLOBAL.SIGNIFICANT_EIGEN_VALUES, 0.0);

  for (count_type i=0; i<community.individuals.size(); ++i){
      int ind_age_group = nodes[community.individuals[i]].age_group;
      age_component[ind_age_group] += nodes[community.individuals[i]].lambda_h;
  }

  for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
    for(count_type count=0; count<size; ++count){
      V_tx[eigen_count] += community_tx_vT[eigen_count][count]
                           * age_component[count];
    }
  }

  for (count_type count=0; count<GLOBAL.NUM_AGE_GROUPS; ++count){
    for (count_type eigen_count=0; eigen_count<GLOBAL.SIGNIFICANT_EIGEN_VALUES; ++eigen_count){
      lambda_age_group[count] += community_tx_u[count][eigen_count]
                          * community_tx_sigma[eigen_count]
                          * V_tx[eigen_count];
    }
	lambda_age_group[count] *=  community.scale;
  }
 return lambda_age_group;
}

double updated_travel_fraction(const vector<agent>& nodes, const int cur_time){
  double infected_distance = 0, total_distance = 0;
  count_type actual_travellers = 0, usual_travellers = 0;

  const auto SIZE = nodes.size();
  const auto MASK_FACTOR = GLOBAL.MASK_FACTOR;




#pragma omp parallel for default(none) shared(nodes, SIZE, cur_time, MASK_FACTOR, GLOBAL) \
  reduction (+: usual_travellers, actual_travellers,  \
			 infected_distance, total_distance)
  for(count_type i = 0; i < SIZE; ++i){
	if(nodes[i].has_to_travel){
	  ++usual_travellers;
	}
	if(nodes[i].travels()){
	  double mask_factor = 1.0;
	  if(mask_active(cur_time) && nodes[i].compliant){
		mask_factor = MASK_FACTOR;
	  }
	  ++actual_travellers;
	  total_distance += nodes[i].commute_distance;
	  if(nodes[i].infective){
		  if(nodes[i].new_strain==1){//---new_Strain replace by shakir. More if conditions are added for other strains
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_ALPHA;
		  }
		  if(nodes[i].new_strain==2){
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_DELTA;
		  }

		  if(nodes[i].new_strain==3){
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_OMICRON;
		  }
		  if(nodes[i].new_strain==4){
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_OMICRON_NEW;
		  }		  		  
		  if(nodes[i].new_strain==5){
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_OMICRON_BA4;
		  }
		  if(nodes[i].new_strain==6){
		infected_distance += nodes[i].commute_distance * mask_factor * GLOBAL.INFECTIOUSNESS_OMICRON_BA5;
		  }		  		  		  		  
		  else{
			  infected_distance += nodes[i].commute_distance * mask_factor;
		  }
	  }
	}
  }
  if(total_distance == 0 || usual_travellers == 0){
	  return 0;
  } else{
	  return (infected_distance/total_distance)
	* double(actual_travellers)/double(usual_travellers);
  }
  
}

vector<double>  updated_travel_fraction_higher(const vector<agent>& nodes, const int cur_time){
  //double infected_distance_higher, total_distance = 0;//---Blocked by shakir for new variants

  //-------lambda_incoming_higher for variants begins--shakir------//
double infected_distance_higher[7] = {0},total_distance = 0;
std::vector<double> zeros;
std::vector<double> travel_fraction_higher;

travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);
travel_fraction_higher.push_back(0);

zeros.push_back(0);
zeros.push_back(0);
zeros.push_back(0);
zeros.push_back(0);
zeros.push_back(0);
zeros.push_back(0);
zeros.push_back(0);
//double travel_fraction_higher[7] ={0};
//------lambda_incoming_higher for variants ends--shakir------//  

  count_type actual_travellers = 0, usual_travellers = 0;

  const auto SIZE = nodes.size();
  const auto MASK_FACTOR = GLOBAL.MASK_FACTOR;




#pragma omp parallel for default(none) shared(nodes, cur_time, MASK_FACTOR, GLOBAL, SIZE) \
  reduction (+: usual_travellers, actual_travellers,  \
			 infected_distance_higher, total_distance)
  for(count_type i = 0; i < SIZE; ++i){
	if(nodes[i].has_to_travel){
	  ++usual_travellers;
	}
	if(nodes[i].travels()){
	  double mask_factor = 1.0;
	  if(mask_active(cur_time) && nodes[i].compliant){
		mask_factor = MASK_FACTOR;
	  }
	  ++actual_travellers;
	  total_distance += nodes[i].commute_distance;
	  if((nodes[i].infective)&&(nodes[i].new_strain==1)){//---new_strain replaced by shakir
		infected_distance_higher[1-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_ALPHA;
	  }
	  if((nodes[i].infective)&&(nodes[i].new_strain==2)){//---new_strain replaced by shakir
		infected_distance_higher[2-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_DELTA;
	  }
	  if((nodes[i].infective)&&(nodes[i].new_strain==3)){//---new_strain replaced by shakir
		infected_distance_higher[3-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_OMICRON;
	  }
	  if((nodes[i].infective)&&(nodes[i].new_strain==4)){//---new_strain replaced by shakir
		infected_distance_higher[4-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_OMICRON_NEW;
	  }
	  if((nodes[i].infective)&&(nodes[i].new_strain==5)){//---new_strain replaced by shakir
		infected_distance_higher[5-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_OMICRON_BA4;
	  }
	  if((nodes[i].infective)&&(nodes[i].new_strain==6)){//---new_strain replaced by shakir
		infected_distance_higher[6-1] += nodes[i].commute_distance * mask_factor  * GLOBAL.INFECTIOUSNESS_OMICRON_BA5;
	  }	  	  	  
	}
  }
  if(total_distance == 0 || usual_travellers == 0){
	  return zeros;//0;
  } else{
	for(int i=1;i<=7;i++){
		travel_fraction_higher[i-1]=infected_distance_higher[i-1]/total_distance 
		* double(actual_travellers)/double(usual_travellers);
	}
	  return (travel_fraction_higher);//(infected_distance_higher/total_distance)
	//* double(actual_travellers)/double(usual_travellers);
  }
  
}


void update_lambdas(agent&node, const vector<house>& homes, const vector<workplace>& workplaces, const vector<community>& communities, const vector<vector<nbr_cell>>& nbr_cells, const double travel_fraction, std::vector<double> travel_fraction_higher, const int cur_time){
  node.lambda_incoming.set_zero();
 // node.lambda_incoming_higher.set_zero();//---delete shakir later

//-------lambda_incoming_higher for variants--shakir------//
  node.lambda_incoming_higher1.set_zero();
  node.lambda_incoming_higher2.set_zero();
  node.lambda_incoming_higher3.set_zero();
  node.lambda_incoming_higher4.set_zero();
  node.lambda_incoming_higher5.set_zero();
  node.lambda_incoming_higher6.set_zero();

//-----------------------------------//

  //Contributions from home, workplace, community, and travel
  if (GLOBAL.USE_AGE_DEPENDENT_MIXING){
    node.lambda_incoming.home = node.kappa_H_incoming
	  * homes[node.home].age_dependent_mixing[node.age_group]
	  * node.hd_area_factor;

    if(node.workplace != WORKPLACE_HOME) {
	  node.lambda_incoming.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_dependent_mixing[node.age_group];
    }
    
  }
  else {
	//No null check for home as every agent has a home
	node.lambda_incoming.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing
	  * node.hd_area_factor;
	  

//-------lambda_incoming_higher for variants begins--shakir------//
	node.lambda_incoming_higher1.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher1
	  * node.hd_area_factor;

	node.lambda_incoming_higher2.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher2
	  * node.hd_area_factor;

	node.lambda_incoming_higher3.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher3
	  * node.hd_area_factor;

	node.lambda_incoming_higher4.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher4
	  * node.hd_area_factor;

	node.lambda_incoming_higher5.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher5
	  * node.hd_area_factor;

	node.lambda_incoming_higher6.home = node.kappa_H_incoming
	  * homes[node.home].age_independent_mixing_higher6
	  * node.hd_area_factor; 
	  
//------lambda_incoming_higher for variants ends--shakir------//

	//If the agent lives in a high population density area, eg, a slum

	//FEATURE_PROPOSAL: make the mixing dependent on node.age_group;
	if(node.workplace != WORKPLACE_HOME) {
	  node.lambda_incoming.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing;
	  //FEATURE_PROPOSAL: make the mixing dependent on node.age_group;
	  

//-------lambda_incoming_higher for variants begins--shakir------//
	  node.lambda_incoming_higher1.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher1;

	  node.lambda_incoming_higher2.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher2;

	  node.lambda_incoming_higher3.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher3;

	  node.lambda_incoming_higher4.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher4;

	  node.lambda_incoming_higher5.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher5;

	  node.lambda_incoming_higher6.work = (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].age_independent_mixing_higher6;

//------lambda_incoming_higher for variants ends--shakir------//

	}
  }
  if(node.workplace != WORKPLACE_HOME){
	  node.lambda_incoming.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing;
  

//-------lambda_incoming_higher for variants begins--shakir------//
  node.lambda_incoming_higher1.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher1;

  node.lambda_incoming_higher2.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher2;

  node.lambda_incoming_higher3.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher3;

  node.lambda_incoming_higher4.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher4;

  node.lambda_incoming_higher5.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher5;

  node.lambda_incoming_higher6.project =  (node.attending?1.0:GLOBAL.ATTENDANCE_LEAKAGE)*node.kappa_W_incoming
		* workplaces[node.workplace].projects[node.workplace_subnetwork].age_independent_mixing_higher6;

//------lambda_incoming_higher for variants ends--shakir------//  
  }
  // No null check for community as every node has a community.
  //
  // For all communities add the community lambda with a distance
  // related scaling factor
  node.lambda_incoming.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);
  //If the agent lives in a high population density area, eg, a slum
  
//   node.lambda_incoming_higher.community = node.kappa_C_incoming
// 	* node.zeta_a
// 	* node.funct_d_ck
// 	* communities[node.community].lambda_community_global_higher
// 	* node.hd_area_factor
// 	* pow(communities[node.community].individuals.size(),
// 		  node.hd_area_exponent);//----delete shakir later---//

//-------lambda_incoming_higher for variants begins--shakir------//
  node.lambda_incoming_higher1.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher1
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);
  
  node.lambda_incoming_higher2.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher2
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);		  

  node.lambda_incoming_higher3.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher3
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);

  node.lambda_incoming_higher4.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher4
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);

  node.lambda_incoming_higher5.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher5
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);

  node.lambda_incoming_higher6.community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* communities[node.community].lambda_community_global_higher6
	* node.hd_area_factor
	* pow(communities[node.community].individuals.size(),
		  node.hd_area_exponent);

//------lambda_incoming_higher for variants ends--shakir------//  
  

  node.lambda_incoming.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community
	* node.hd_area_factor;
	
   node.lambda_incoming_higher.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher
	* node.hd_area_factor;

//-------lambda_incoming_higher for variants begins--shakir------//
   node.lambda_incoming_higher1.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher1
	* node.hd_area_factor;

   node.lambda_incoming_higher2.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher2
	* node.hd_area_factor;

   node.lambda_incoming_higher3.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher3
	* node.hd_area_factor;

   node.lambda_incoming_higher4.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher4
	* node.hd_area_factor;

   node.lambda_incoming_higher5.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher5
	* node.hd_area_factor;

   node.lambda_incoming_higher6.random_community = node.kappa_C_incoming
	* node.zeta_a
	* node.funct_d_ck
	* homes[node.home].random_households.lambda_random_community_higher6
	* node.hd_area_factor;
										

//------lambda_incoming_higher for variants ends--shakir------//  
//std::cout<<"Neighbour size"<<"\t"<<nbr_cells.size()<<std::endl;

  if(nbr_cells.size()>0){
	node.lambda_incoming.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher
	  * node.hd_area_factor;
 
//-------lambda_incoming_higher for variants begins--shakir------//
	
	node.lambda_incoming_higher1.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher1
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher2.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher2
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher3.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher3
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher4.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher4
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher5.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher5
	  * node.hd_area_factor;
	
	node.lambda_incoming_higher6.nbr_cell = node.kappa_C_incoming
	  * node.zeta_a
	  * nbr_cells[homes[node.home].neighbourhood.cell_x][homes[node.home].neighbourhood.cell_y].lambda_nbr_higher6
	  * node.hd_area_factor;
      

//------lambda_incoming_higher for variants ends--shakir------//  

  }
  else{
	node.lambda_incoming.nbr_cell = 0;
	node.lambda_incoming_higher.nbr_cell = 0;

//-------lambda_incoming_higher for variants begins--shakir------//
	node.lambda_incoming_higher1.nbr_cell = 0;
	node.lambda_incoming_higher2.nbr_cell = 0;
	node.lambda_incoming_higher3.nbr_cell = 0;
	node.lambda_incoming_higher4.nbr_cell = 0;
	node.lambda_incoming_higher5.nbr_cell = 0;
	node.lambda_incoming_higher6.nbr_cell = 0;

//------lambda_incoming_higher for variants ends--shakir------//  

  }


  //Travel only happens at "odd" times, twice a day
//std::cout<<"curr time, cur_time%2, has to travel, attending, compliant, and quarantined "<<cur_time<<"\t"<<cur_time%2<<"\t"<<node.has_to_travel<<"\t"<<node.has_to_travel<<"\t"<<node.attending<<"\t"<<node.compliant<<"\t"<<node.quarantined<<"\t"<<(!(node.quarantined && node.compliant))<<"\t"<<node.travels()<<std::endl;

  if((cur_time % 2) && node.travels()){
	node.lambda_incoming.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction;
	
//std::cout<<"travel force"<<node.lambda_incoming.travel<<std::endl;


//-------lambda_incoming_higher for variants begins--shakir------//
 // std::cout<<"travel force ends"<<node.lambda_incoming.travel<<"\t travel fraction vector casing"<<"\t"<<travel_fraction_higher[1-1]<<std::endl;
	  
	node.lambda_incoming_higher1.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[1-1];
 // std::cout<<"travel force ends"<<node.lambda_incoming.travel<<std::endl;
	  
	node.lambda_incoming_higher2.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[2-1];
	  
	node.lambda_incoming_higher3.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[3-1];
	  
	node.lambda_incoming_higher4.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[4-1];
	  
	node.lambda_incoming_higher5.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[5-1];
	  
	node.lambda_incoming_higher6.travel = GLOBAL.BETA_TRAVEL
	  * node.commute_distance
	  * travel_fraction_higher[6-1];
 // std::cout<<"travel force ends"<<node.lambda_incoming.travel<<std::endl;


//------lambda_incoming_higher for variants ends--shakir------//  

  }

  if(mask_active(cur_time) && node.compliant){
	node.lambda_incoming.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming.nbr_cell *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher.nbr_cell *= GLOBAL.MASK_FACTOR;
//-------lambda_incoming_higher for variants begins--shakir------//
	node.lambda_incoming_higher1.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher1.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher1.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher1.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher1.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher1.nbr_cell *= GLOBAL.MASK_FACTOR;

	node.lambda_incoming_higher2.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher2.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher2.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher2.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher2.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher2.nbr_cell *= GLOBAL.MASK_FACTOR;

	node.lambda_incoming_higher3.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher3.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher3.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher3.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher3.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher3.nbr_cell *= GLOBAL.MASK_FACTOR;

	node.lambda_incoming_higher4.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher4.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher4.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher4.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher4.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher4.nbr_cell *= GLOBAL.MASK_FACTOR;

	node.lambda_incoming_higher5.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher5.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher5.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher5.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher5.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher5.nbr_cell *= GLOBAL.MASK_FACTOR;

	node.lambda_incoming_higher6.work *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher6.community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher6.travel *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher6.project *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher6.random_community *= GLOBAL.MASK_FACTOR;
	node.lambda_incoming_higher6.nbr_cell *= GLOBAL.MASK_FACTOR;

//------lambda_incoming_higher for variants ends--shakir------//  

  }

  node.lambda = node.lambda_incoming.sum();
  //node.lambda_higher = node.lambda_incoming_higher.sum();//---delete shakir later after verifying
//-------lambda_incoming_higher for variants begins--shakir------//
  node.lambda_higher1 = node.lambda_incoming_higher1.sum();
  node.lambda_higher2 = node.lambda_incoming_higher2.sum();
  node.lambda_higher3 = node.lambda_incoming_higher3.sum();
  node.lambda_higher4 = node.lambda_incoming_higher4.sum();
  node.lambda_higher5 = node.lambda_incoming_higher5.sum();
  node.lambda_higher6 = node.lambda_incoming_higher6.sum();

//------lambda_incoming_higher for variants ends--shakir------//  

}
// sk ///////////////////////
void update_individual_lambda_cohort(vector<agent>&nodes, const int cur_time, std::unordered_map<count_type, vector<cohort_space>>& cohorts){
	if (!GLOBAL.ENABLE_COHORTS || !GLOBAL.TRAINS_RUNNING || ((cur_time % 2 )==0) ){
		return;
	}
	for (auto& it1:cohorts){
		for (auto& it2:it1.second){
			for (auto indiv:it2.internal_nodes){				
				auto& node=nodes[indiv];
				if(node.attending){
					double sum_interactions = it2.lambda_interaction_internal * it2.commute_time + it2.lambda_interaction_external;
					//FIX: Possible problem related to GLOBAL.SIM_STEPS_PER_DAY. SHould cohorts lambda scale with GLOBAL.SIM_STEPS_PER_DAY?
					node.lambda_incoming.cohorts = sum_interactions * node.kappa_W_incoming; //
				
					if(mask_active(cur_time) && node.compliant){	
						node.lambda_incoming.cohorts *= GLOBAL.MASK_FACTOR;
  					}
				}				
				node.lambda = node.lambda_incoming.sum();
			}
		}
	}
	return;
}
/////////////////////////////
void updated_lambda_c_local(const vector<agent>& nodes, community& community){
  double sum_value = 0;
  //double sum_value_higher = 0,sum_alpha=0,sum_delta=0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;
  //-------lambda_incoming_higher for variants begins--shakir------//
double sum_value_higher[6]={0};

//------lambda_incoming_higher for variants ends--shakir------//  

  const auto SIZE = community.individuals.size();

#pragma omp parallel for default(none) shared(nodes, community) reduction (+: sum_value)
  for(count_type i = 0; i < SIZE; ++i){
	sum_value
	  += nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	if (nodes[community.individuals[i]].new_strain==1){			 
	//sum_alpha
	sum_value_higher[1-1]
	  += nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}

	if (nodes[community.individuals[i]].new_strain==2){			 
	//sum_delta
	  sum_value_higher[2-1]+= nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}
	if (nodes[community.individuals[i]].new_strain==3){			 
	//sum_omicron
	  sum_value_higher[3-1]+= nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}
	if (nodes[community.individuals[i]].new_strain==4){			 
	//sum_omicronnew
	  sum_value_higher[4-1]+= nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}
	if (nodes[community.individuals[i]].new_strain==5){			 
	//sum_omicronBA4
	  sum_value_higher[5-1]+= nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}

	if (nodes[community.individuals[i]].new_strain==6){			 
	//sum_omicronBA5
	  sum_value_higher[6-1]+= nodes[community.individuals[i]].lambda_c
	  * std::min(community.w_c,
				 nodes[community.individuals[i]].neighborhood_access_factor);
	}							 
	// if (nodes[community.individuals[i]].new_strain){---orignal new_strain if condition blocked by shakir			 
	// sum_value_higher
	//   += nodes[community.individuals[i]].lambda_c
	//   * std::min(community.w_c,
	// 			 nodes[community.individuals[i]].neighborhood_access_factor);
	// }				 
  }


	// if(sum_alpha>sum_delta && sum_alpha>sum_omicron && sum_alpha>sum_omicronnew
	// 	&& sum_alpha>sum_omicronBA4 && sum_alpha>sum_omicronBA5)//---these comparisons to find maximum introduced by shakir
	// {
	// sum_value_higher=sum_alpha;
	// dominant_var=1;
	// }

	// if(sum_delta>sum_alpha && sum_delta>sum_omicron && sum_delta>sum_omicronnew
	// && sum_delta>sum_omicronBA4 && sum_delta>sum_omicronBA5)
	// {
	// sum_value_higher=sum_delta;
	// dominant_var=2;
	// }
	// if(sum_omicron>sum_alpha && sum_omicron>sum_delta && sum_omicron>sum_omicronnew
	// && sum_omicron>sum_omicronBA4 && sum_omicron>sum_omicronBA5)
	// {
	// sum_value_higher=sum_omicron;
	// dominant_var=3;
	// }
	// if(sum_omicronnew>sum_alpha && sum_omicronnew>sum_delta && sum_omicronnew>sum_omicron
	// && sum_omicronnew>sum_omicronBA4 && sum_omicronnew>sum_omicronBA5)
	// {
	// sum_value_higher=sum_omicronnew;
	// dominant_var=4;
	// }
	// if(sum_omicronBA4>sum_alpha && sum_omicronBA4>sum_delta && sum_omicronBA4>sum_omicron
	// && sum_omicronBA4>sum_omicronnew && sum_omicronBA4>sum_omicronBA5)
	// {
	// sum_value_higher=sum_omicronBA4;
	// dominant_var=5;
	// }
	// if(sum_omicronBA5>sum_alpha && sum_omicronBA5>sum_delta && sum_omicronBA5>sum_omicron
	// && sum_omicronBA5>sum_omicronnew && sum_omicronBA5>sum_omicronBA4)
	// {
	// sum_value_higher=sum_omicronBA5;
	// dominant_var=6;
	// }

  community.lambda_community = community.scale*sum_value;
  //community.lambda_community_higher = community.scale*sum_value_higher;
//-------lambda_incoming_higher for variants begins--shakir------//
  community.lambda_community_higher1 = community.scale*sum_value_higher[1-1];
  community.lambda_community_higher2 = community.scale*sum_value_higher[2-1];
  community.lambda_community_higher3 = community.scale*sum_value_higher[3-1];
  community.lambda_community_higher4 = community.scale*sum_value_higher[4-1];
  community.lambda_community_higher5 = community.scale*sum_value_higher[5-1];
  community.lambda_community_higher6 = community.scale*sum_value_higher[6-1];  

//------lambda_incoming_higher for variants ends--shakir------//  


}

void updated_lambda_c_local_random_community(const vector<agent>& nodes, const vector<community>& communities, vector<house>& houses){
  const auto HOUSES_SIZE = houses.size();
#pragma omp parallel for default(none) shared(houses, nodes)
  for(count_type i = 0;  i < HOUSES_SIZE; ++i){
	double lambda_random_community_outgoing = 0;
//	double lambda_random_community_outgoing_higher = 0,sum_alpha=0,sum_delta=0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;

//-------lambda_incoming_higher for variants begins--shakir------//
	double lambda_random_community_outgoing_higher[6] = {0};

//------lambda_incoming_higher for variants ends--shakir------//  


	for(const auto& indiv: houses[i].individuals){
	  lambda_random_community_outgoing += nodes[indiv].lambda_c;
	//   if (nodes[indiv].new_strain){ original new_strain if condition blocked by shakir
	//      lambda_random_community_outgoing_higher += nodes[indiv].lambda_c;  
	//   }
	  if (nodes[indiv].new_strain==1){
	     lambda_random_community_outgoing_higher[1-1] += nodes[indiv].lambda_c;  
	  }
	  if (nodes[indiv].new_strain==2){
	     lambda_random_community_outgoing_higher[2-1] += nodes[indiv].lambda_c;  
	  }
	  if (nodes[indiv].new_strain==3){
	     lambda_random_community_outgoing_higher[3-1] += nodes[indiv].lambda_c;  
	  }
	  if (nodes[indiv].new_strain==4){
	     lambda_random_community_outgoing_higher[4-1] += nodes[indiv].lambda_c;  
	  }	  	  	  	  
	  if (nodes[indiv].new_strain==5){
	     lambda_random_community_outgoing_higher[5-1] += nodes[indiv].lambda_c;  
	}
	  if (nodes[indiv].new_strain==6){
	     lambda_random_community_outgoing_higher[6-1] += nodes[indiv].lambda_c;  
	  }		  		   	  	  	  
	}



    // if(sum_alpha>sum_delta && sum_alpha>sum_omicron && sum_alpha>sum_omicronnew
    //     && sum_alpha>sum_omicronBA4 && sum_alpha>sum_omicronBA5)//---these comparisons to find maximum introduced by shakir
	// {
    // lambda_random_community_outgoing_higher=sum_alpha;
	// dominant_var=1;
	// }
    // if(sum_delta>sum_alpha && sum_delta>sum_omicron && sum_delta>sum_omicronnew
    // && sum_delta>sum_omicronBA4 && sum_delta>sum_omicronBA5)
	// {
    // lambda_random_community_outgoing_higher=sum_delta;
	// dominant_var=2;
	// }
    // if(sum_omicron>sum_alpha && sum_omicron>sum_delta && sum_omicron>sum_omicronnew
    // && sum_omicron>sum_omicronBA4 && sum_omicron>sum_omicronBA5)
	// {
    // lambda_random_community_outgoing_higher=sum_omicron;
	// dominant_var=3;
	// }
    // if(sum_omicronnew>sum_alpha && sum_omicronnew>sum_delta && sum_omicronnew>sum_omicron
    // && sum_omicronnew>sum_omicronBA4 && sum_omicronnew>sum_omicronBA5)
	// {
    // lambda_random_community_outgoing_higher=sum_omicronnew;
	// dominant_var=4;
	// }
    // if(sum_omicronBA4>sum_alpha && sum_omicronBA4>sum_delta && sum_omicronBA4>sum_omicron
    // && sum_omicronBA4>sum_omicronnew && sum_omicronBA4>sum_omicronBA5)
	// {
    // lambda_random_community_outgoing_higher=sum_omicronBA4;
	// dominant_var=5;
	// }
    // if(sum_omicronBA5>sum_alpha && sum_omicronBA5>sum_delta && sum_omicronBA5>sum_omicron
    // && sum_omicronBA5>sum_omicronnew && sum_omicronBA5>sum_omicronBA4)
	// {
    // lambda_random_community_outgoing_higher=sum_omicronBA5;
	// dominant_var=6;
	// }
//------------------------------------------------------------------------//

	houses[i].lambda_random_community_outgoing = lambda_random_community_outgoing;
//	houses[i].lambda_random_community_outgoing_higher = lambda_random_community_outgoing_higher;---blocked by shakir--to be deleted later
//-------lambda_incoming_higher for variants begins--shakir------//

	houses[i].lambda_random_community_outgoing_higher1 = lambda_random_community_outgoing_higher[1-1];
	houses[i].lambda_random_community_outgoing_higher2 = lambda_random_community_outgoing_higher[2-1];
	houses[i].lambda_random_community_outgoing_higher3 = lambda_random_community_outgoing_higher[3-1];
	houses[i].lambda_random_community_outgoing_higher4 = lambda_random_community_outgoing_higher[4-1];
	houses[i].lambda_random_community_outgoing_higher5 = lambda_random_community_outgoing_higher[5-1];
	houses[i].lambda_random_community_outgoing_higher6 = lambda_random_community_outgoing_higher[6-1];

//------lambda_incoming_higher for variants ends--shakir------//  

  }
#pragma omp parallel for default(none) shared(houses, communities)
  for(count_type i = 0; i < HOUSES_SIZE; ++i){
	double sum_value_household = 0;
	double sum_value_household_higher[6] = {0};
	for(const auto& neighbouring_household: houses[i].random_households.households){
	  sum_value_household += houses[neighbouring_household].lambda_random_community_outgoing;
//	  sum_value_household_higher += houses[neighbouring_household].lambda_random_community_outgoing_higher;
	  
//-------lambda_incoming_higher for variants begins--shakir------//
	  sum_value_household_higher[1-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher1;
	  sum_value_household_higher[2-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher2;
	  sum_value_household_higher[3-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher3;
	  sum_value_household_higher[4-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher4;
	  sum_value_household_higher[5-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher5;
	  sum_value_household_higher[6-1] += houses[neighbouring_household].lambda_random_community_outgoing_higher6;


//------lambda_incoming_higher for variants ends--shakir------//  

	}
	houses[i].random_households.lambda_random_community = houses[i].random_households.scale
	  * sum_value_household
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);
	  
	// houses[i].random_households.lambda_random_community_higher = houses[i].random_households.scale
	//   * sum_value_household_higher
	//   * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);//-----delete shakir later
//-------lambda_incoming_higher for variants begins--shakir------//
	houses[i].random_households.lambda_random_community_higher1 = houses[i].random_households.scale
	  * sum_value_household_higher[1-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);
    
	houses[i].random_households.lambda_random_community_higher2 = houses[i].random_households.scale
	  * sum_value_household_higher[2-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);

	houses[i].random_households.lambda_random_community_higher3 = houses[i].random_households.scale
	  * sum_value_household_higher[3-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);

	houses[i].random_households.lambda_random_community_higher4 = houses[i].random_households.scale
	  * sum_value_household_higher[4-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);

	houses[i].random_households.lambda_random_community_higher5 = houses[i].random_households.scale
	  * sum_value_household_higher[5-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);

	houses[i].random_households.lambda_random_community_higher6 = houses[i].random_households.scale
	  * sum_value_household_higher[6-1]
	  * std::min(communities[houses[i].community].w_c, houses[i].neighborhood_access_factor);
//------lambda_incoming_higher for variants ends--shakir------//  
    
  }
}

void update_lambda_nbr_cells(const vector<agent>& nodes, vector<vector<nbr_cell>>& nbr_cells, const vector<house>& houses, const vector<community>& communities){
  for(count_type i=0; i<nbr_cells.size(); ++i){
	for(count_type j=0; j<nbr_cells[i].size(); ++j){
	  double sum_values = 0;
//	  double sum_values_higher = 0,sum_alpha=0,sum_delta=0,sum_omicron=0,sum_omicronnew=0,sum_omicronBA4=0,sum_omicronBA5=0;

//-------lambda_incoming_higher for variants begins--shakir------//
	  double sum_values_higher[6] = {0};
//------lambda_incoming_higher for variants ends--shakir------//  

#pragma omp parallel for default(none)					\
  shared(nbr_cells, communities, nodes,					\
		 houses, i, j)									\
  reduction (+: sum_values)
	  for(count_type h=0; h<nbr_cells[i][j].houses_list.size(); ++h){
		const auto house_index = nbr_cells[i][j].houses_list[h];
		for(count_type k=0; k<houses[house_index].individuals.size(); ++k){
		  sum_values += nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			* std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
//			if (nodes[houses[house_index].individuals[k]].new_strain)//-----Changed by shakir new_strain to alpha_strain
if (nodes[houses[house_index].individuals[k]].new_strain==1)
			{
				//sum_alpha
				 sum_values_higher[1-1]+= nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
			}
if (nodes[houses[house_index].individuals[k]].new_strain==2)
			{
				//sum_delta
				sum_values_higher[2-1] += nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
			}
if (nodes[houses[house_index].individuals[k]].new_strain==3)
			{
				//sum_omicron
				sum_values_higher[3-1] += nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
			}
if (nodes[houses[house_index].individuals[k]].new_strain==4)
			{
				//sum_omicronnew 
				sum_values_higher[4-1]+= nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
			}													
if (nodes[houses[house_index].individuals[k]].new_strain==5)
			{
				//sum_omicronBA4 
				sum_values_higher[5-1]+= nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
		}
if (nodes[houses[house_index].individuals[k]].new_strain==6)
			{
				//sum_omicronBA5
				sum_values_higher[6-1] += nodes[houses[house_index].individuals[k]].lambda_nbr_cell
			        * std::min(communities[houses[house_index].community].w_c,
					   houses[house_index].neighborhood_access_factor);
	  }
		}
	  }




	  nbr_cells[i][j].lambda_nbr = nbr_cells[i][j].scale*sum_values;
	  
	  //nbr_cells[i][j].lambda_nbr_higher = nbr_cells[i][j].scale*sum_values_higher;//----delete shakir later
//-------lambda_incoming_higher for variants begins--shakir------//

	  nbr_cells[i][j].lambda_nbr_higher1 = nbr_cells[i][j].scale*sum_values_higher[1-1];
	  nbr_cells[i][j].lambda_nbr_higher2 = nbr_cells[i][j].scale*sum_values_higher[2-1];
	  nbr_cells[i][j].lambda_nbr_higher3 = nbr_cells[i][j].scale*sum_values_higher[3-1];
	  nbr_cells[i][j].lambda_nbr_higher4 = nbr_cells[i][j].scale*sum_values_higher[4-1];
	  nbr_cells[i][j].lambda_nbr_higher5 = nbr_cells[i][j].scale*sum_values_higher[5-1];
	  nbr_cells[i][j].lambda_nbr_higher6 = nbr_cells[i][j].scale*sum_values_higher[6-1];

//------lambda_incoming_higher for variants ends--shakir------//  

	}
  }
}


void update_lambda_c_global(vector<community>& communities,
							const matrix<double>& community_distance_fk_matrix){
  const auto SIZE = communities.size();
  for (count_type c1 = 0; c1 < SIZE; ++c1){
	double num = 0;
	//double num_higher = 0;
//-------lambda_incoming_higher for variants begins--shakir------//
	double num_higher[6] = {0};


//------lambda_incoming_higher for variants ends--shakir------//  

	double denom = 0;

	for (count_type c2 = 0; c2 < SIZE; ++c2){
	  double fk_val = community_distance_fk_matrix[c1][c2];
	  num += fk_val * communities[c2].lambda_community;
	  //num_higher 
//-------lambda_incoming_higher for variants begins--shakir------//
	  num_higher[1-1]+= fk_val * communities[c2].lambda_community_higher1;
	  num_higher[2-1]+= fk_val * communities[c2].lambda_community_higher2;
	  num_higher[3-1]+= fk_val * communities[c2].lambda_community_higher3;
	  num_higher[4-1]+= fk_val * communities[c2].lambda_community_higher4;
	  num_higher[5-1]+= fk_val * communities[c2].lambda_community_higher5;
	  num_higher[6-1]+= fk_val * communities[c2].lambda_community_higher6;


//------lambda_incoming_higher for variants ends--shakir------//  
	  
	  denom += fk_val;
	}
	if(denom==0){		
		communities[c1].lambda_community_global = 0;
//		communities[c1].lambda_community_global_higher = 0;
//-------lambda_incoming_higher for variants begins--shakir------//
		communities[c1].lambda_community_global_higher1 = 0;
		communities[c1].lambda_community_global_higher2 = 0;
		communities[c1].lambda_community_global_higher3 = 0;
		communities[c1].lambda_community_global_higher4 = 0;
		communities[c1].lambda_community_global_higher5 = 0;
		communities[c1].lambda_community_global_higher6 = 0;

//------lambda_incoming_higher for variants ends--shakir------//  

	} else{		
		communities[c1].lambda_community_global = communities[c1].w_c*num/denom;
		//communities[c1].lambda_community_global_higher = communities[c1].w_c*num_higher/denom;
//-------lambda_incoming_higher for variants begins--shakir------//
		communities[c1].lambda_community_global_higher1 = communities[c1].w_c*num_higher[1-1]/denom;
		communities[c1].lambda_community_global_higher2 = communities[c1].w_c*num_higher[2-1]/denom;
		communities[c1].lambda_community_global_higher3 = communities[c1].w_c*num_higher[3-1]/denom;
		communities[c1].lambda_community_global_higher4 = communities[c1].w_c*num_higher[4-1]/denom;
		communities[c1].lambda_community_global_higher5 = communities[c1].w_c*num_higher[5-1]/denom;
		communities[c1].lambda_community_global_higher6 = communities[c1].w_c*num_higher[6-1]/denom;

//------lambda_incoming_higher for variants ends--shakir------//  

	}
	
  }
}

void update_test_request(vector<agent>& nodes, const vector<house>& homes,
						 const vector<workplace>& workplaces, const vector<community>& communities,
						 vector<vector<nbr_cell>>& nbr_cells, const count_type current_time, const vector<testing_probability>& testing_protocol){
  testing_probability probabilities;
  if(current_time >= GLOBAL.NUM_DAYS_BEFORE_INTERVENTIONS*GLOBAL.SIM_STEPS_PER_DAY){
	switch(GLOBAL.TESTING_PROTOCOL){
	case Testing_Protocol::no_testing:
		break;
	case Testing_Protocol::test_household:
		probabilities.prob_test_index_hospitalised = 1;
		probabilities.prob_test_household_symptomatic_symptomatic = 0;
		probabilities.prob_test_household_symptomatic_asymptomatic = 0;

		probabilities.prob_test_household_hospitalised_symptomatic = 1;
		probabilities.prob_test_household_hospitalised_asymptomatic = 0;
		probabilities.prob_test_household_positive_symptomatic = 1;
		probabilities.prob_test_household_positive_asymptomatic = 0;

		probabilities.prob_test_neighbourhood_hospitalised_symptomatic = 0;
		probabilities.prob_test_neighbourhood_hospitalised_asymptomatic = 0;
		probabilities.prob_test_neighbourhood_positive_symptomatic = 0;
		probabilities.prob_test_neighbourhood_positive_asymptomatic = 0;

		probabilities.prob_contact_trace_household_hospitalised = 1;
		probabilities.prob_contact_trace_household_positive = 1;

		probabilities.prob_retest_recovered = 1;
		set_test_request(nodes, homes, workplaces, nbr_cells, communities, probabilities, current_time);
		break;
	case Testing_Protocol::testing_protocol_file_read:
		set_test_request_fileread(nodes, homes, workplaces, nbr_cells, communities, testing_protocol, current_time);
		break;
	default:
		break;
	}
  }
}

void update_test_status(vector<agent>& nodes, count_type current_time){
  for(auto& node: nodes){
    if(node.test_status.test_requested){
	  if(node.infection_status == Progression::infective
		 || node.infection_status == Progression::symptomatic
		 || node.infection_status == Progression::hospitalised
		 || node.infection_status == Progression::critical){
		node.test_status.state = bernoulli(GLOBAL.TEST_FALSE_NEGATIVE)?test_result::negative:test_result::positive;
		node.test_status.tested_positive = node.test_status.tested_positive || (node.test_status.state == test_result::positive);
		node.test_status.tested_epoch = current_time;
	  }
	  else if(node.infection_status == Progression::exposed
			  && current_time-node.time_of_infection > GLOBAL.SIM_STEPS_PER_DAY*GLOBAL.TIME_TO_TEST_POSITIVE){
		node.test_status.state = bernoulli(GLOBAL.TEST_FALSE_NEGATIVE)?test_result::negative:test_result::positive;
		node.test_status.tested_positive = node.test_status.tested_positive || (node.test_status.state == test_result::positive);
		//We might want to have higher false negative rate here, depending upon updates in the data.
		node.test_status.tested_epoch = current_time;
	  }
	  else{
		// Test could come positive for a succeptible/recovered/dead person
		node.test_status.state = bernoulli(GLOBAL.TEST_FALSE_POSITIVE)?test_result::positive:test_result::negative;
		node.test_status.tested_positive = node.test_status.tested_positive || (node.test_status.state == test_result::positive);
		node.test_status.tested_epoch = current_time;
	  }
	  node.test_status.test_requested = false;
    }
  }
}
casualty_stats get_infected_community(const vector<agent>& nodes, const community& community){
  count_type affected = 0;
  count_type hd_area_affected = 0;
  count_type susceptible = 0;

  count_type unvaccinated = 0;
  count_type vaccinated1=0;
  count_type vaccinated2=0;
  count_type waning=0;
  count_type boosted=0;
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
  count_type recovered_from_infective = 0;
  count_type recovered_from_symptomatic = 0;
  count_type recovered_from_hospitalised = 0;
  count_type recovered_from_critical = 0;
  count_type recovered_from_vaccination=0;
  count_type hd_area_recovered_from_infective = 0;
  count_type hd_area_recovered_from_symptomatic = 0;
  count_type hd_area_recovered_from_hospitalised = 0;
  count_type hd_area_recovered_from_critical = 0;
  count_type hd_area_recovered_from_vaccination=0;

  count_type errors = 0;
  
  //---Age index counts----/Introduced by Shakir 0n May 4, 2022.
  count_type  infected_age_group_1 = 0;
  count_type  infected_age_group_2 = 0;
  count_type  infected_age_group_3 = 0;
  count_type  infected_age_group_4 = 0;//infected=exposed + infective + symptomatic+ hospitalised + critical+dead
  count_type  infected_age_group_5 = 0;
  count_type  infected_age_group_6 = 0;
  count_type  infected_age_group_7 = 0;
  count_type  infected_age_group_8 = 0;
  count_type  infected_age_group_9 = 0;
  count_type  infected_age_group_10 = 0;


  count_type  hospitalised_age_group_1 = 0;
  count_type  hospitalised_age_group_2 = 0;
  count_type  hospitalised_age_group_3 = 0;
  count_type  hospitalised_age_group_4 = 0;
  count_type  hospitalised_age_group_5 = 0;
  count_type  hospitalised_age_group_6 = 0;
  count_type  hospitalised_age_group_7 = 0;
  count_type  hospitalised_age_group_8 = 0;
  count_type  hospitalised_age_group_9 = 0;
  count_type  hospitalised_age_group_10 = 0;

  count_type  dead_age_group_1 = 0;
  count_type  dead_age_group_2 = 0;
  count_type  dead_age_group_3 = 0;
  count_type  dead_age_group_4 = 0;
  count_type  dead_age_group_5 = 0;
  count_type  dead_age_group_6 = 0;
  count_type  dead_age_group_7 = 0;
  count_type  dead_age_group_8 = 0;
  count_type  dead_age_group_9 = 0;
  count_type  dead_age_group_10 = 0;

  //---income index counts----/Introduced by Shakir 0n May 4, 2022.
  count_type  infected_income_group_1 = 0;
  count_type  infected_income_group_2 = 0;
  count_type  infected_income_group_3 = 0;
  count_type  infected_income_group_4 = 0;//infected=exposed + infective + symptomatic+ hospitalised + critical+dead
  count_type  infected_income_group_5 = 0;
  count_type  infected_income_group_6 = 0;
  count_type  infected_income_group_7 = 0;
  count_type  infected_income_group_8 = 0;


  count_type  hospitalised_income_group_1 = 0;
  count_type  hospitalised_income_group_2 = 0;
  count_type  hospitalised_income_group_3 = 0;
  count_type  hospitalised_income_group_4 = 0;
  count_type  hospitalised_income_group_5 = 0;
  count_type  hospitalised_income_group_6 = 0;
  count_type  hospitalised_income_group_7 = 0;
  count_type  hospitalised_income_group_8 = 0;


  count_type  dead_income_group_1 = 0;
  count_type  dead_income_group_2 = 0;
  count_type  dead_income_group_3 = 0;
  count_type  dead_income_group_4 = 0;
  count_type  dead_income_group_5 = 0;
  count_type  dead_income_group_6 = 0;
  count_type  dead_income_group_7 = 0;
  count_type  dead_income_group_8 = 0;

  //---race index counts----/Introduced by Shakir 0n May 4, 2022.
  count_type  infected_race_group_1 = 0;
  count_type  infected_race_group_2 = 0;
  count_type  infected_race_group_3 = 0;
  count_type  infected_race_group_4 = 0;//infected=exposed + infective + symptomatic+ hospitalised + critical+dead
  count_type  infected_race_group_5 = 0;
  count_type  infected_race_group_6 = 0;
  count_type  infected_race_group_7 = 0;

  count_type  hospitalised_race_group_1 = 0;
  count_type  hospitalised_race_group_2 = 0;
  count_type  hospitalised_race_group_3 = 0;
  count_type  hospitalised_race_group_4 = 0;
  count_type  hospitalised_race_group_5 = 0;
  count_type  hospitalised_race_group_6 = 0;
  count_type  hospitalised_race_group_7 = 0;

  count_type  dead_race_group_1 = 0;
  count_type  dead_race_group_2 = 0;
  count_type  dead_race_group_3 = 0;
  count_type  dead_race_group_4 = 0;
  count_type  dead_race_group_5 = 0;
  count_type  dead_race_group_6 = 0;
  count_type  dead_race_group_7 = 0;

  //---ethnicity index counts----/Introduced by Shakir 0n May 4, 2022.
  count_type  infected_ethnicity_group_1 = 0;
  count_type  infected_ethnicity_group_2 = 0;

  count_type  hospitalised_ethnicity_group_1 = 0;
  count_type  hospitalised_ethnicity_group_2 = 0;

  count_type  dead_ethnicity_group_1 = 0;
  count_type  dead_ethnicity_group_2 = 0;

    //---gender index counts----/Introduced by Shakir 0n May 4, 2022.
  count_type  infected_gender_group_1 = 0;
  count_type  infected_gender_group_2 = 0;

  count_type  hospitalised_gender_group_1 = 0;
  count_type  hospitalised_gender_group_2 = 0;

  count_type  dead_gender_group_1 = 0;
  count_type  dead_gender_group_2 = 0;

  const auto SIZE = community.individuals.size(); 

#pragma omp parallel for default(none) shared(nodes, community, SIZE)			\
  reduction(+: errors,													\
			susceptible,unvaccinated,vaccinated1,vaccinated2,waning,boosted,boosted2, hd_area_susceptible,							\
			exposed, hd_area_exposed,									\
			infective, hd_area_infective,								\
			symptomatic, hd_area_symptomatic,							\
			hospitalised, hd_area_hospitalised,							\
			critical, hd_area_critical,									\
			dead, hd_area_dead,											\
			recovered, hd_area_recovered,								\
			recovered_from_infective, recovered_from_symptomatic,		\
			recovered_from_hospitalised, recovered_from_critical,		\
			recovered_from_vaccination,                                  \
			hd_area_recovered_from_infective,							\
			hd_area_recovered_from_symptomatic,							\
			hd_area_recovered_from_hospitalised,						\
			hd_area_recovered_from_critical,                            \
			infected_age_group_1,infected_age_group_2,infected_age_group_3,infected_age_group_4,infected_age_group_5,		\   
			infected_age_group_6,infected_age_group_7,infected_age_group_8,infected_age_group_9,infected_age_group_10,      \
			hospitalised_age_group_1,hospitalised_age_group_2,hospitalised_age_group_3,hospitalised_age_group_4,hospitalised_age_group_5,	\   
			hospitalised_age_group_6,hospitalised_age_group_7,hospitalised_age_group_8,hospitalised_age_group_9,hospitalised_age_group_10,	\
			dead_age_group_1,dead_age_group_2,dead_age_group_3,dead_age_group_4,dead_age_group_5,		\   
			dead_age_group_6,dead_age_group_7,dead_age_group_8,dead_age_group_9,dead_age_group_10,      \	
			infected_race_group_1,infected_race_group_2,infected_race_group_3,infected_race_group_4,infected_race_group_5,		\   
			infected_race_group_6,infected_race_group_7,    \
			hospitalised_race_group_1,hospitalised_race_group_2,hospitalised_race_group_3,hospitalised_race_group_4,hospitalised_race_group_5,	\   
			hospitalised_race_group_6,hospitalised_race_group_7,	\
			dead_race_group_1,dead_race_group_2,dead_race_group_3,dead_race_group_4,dead_race_group_5,		\   
			dead_race_group_6,dead_race_group_7,    \				
			infected_income_group_1,infected_income_group_2,infected_income_group_3,infected_income_group_4,infected_income_group_5,		\   
			infected_income_group_6,infected_income_group_7,infected_income_group_8,     \
			hospitalised_income_group_1,hospitalised_income_group_2,hospitalised_income_group_3,hospitalised_income_group_4,hospitalised_income_group_5,	\   
			hospitalised_income_group_6,hospitalised_income_group_7,hospitalised_income_group_8,	\
			dead_income_group_1,dead_income_group_2,dead_income_group_3,dead_income_group_4,dead_income_group_5,		\   
			dead_income_group_6,dead_income_group_7,dead_income_group_8,    \
			infected_ethnicity_group_1,infected_ethnicity_group_2,  \
			hospitalised_ethnicity_group_1,hospitalised_ethnicity_group_2,	\
			dead_ethnicity_group_1,dead_ethnicity_group_2,\
			infected_gender_group_1,infected_gender_group_2,  \
			hospitalised_gender_group_1,hospitalised_gender_group_2,	\
			dead_gender_group_1,dead_gender_group_2,\												
			hd_area_recovered_from_vaccination)//added age_group_1-10, and also added unvaccinated group for outputs Shakir May 4
  for (count_type i=0; i<SIZE; ++i){
	bool hd_area_resident = nodes[community.individuals[i]].hd_area_resident;
	auto infection_status = nodes[community.individuals[i]].infection_status;
    int age_index = nodes[community.individuals[i]].age_index;
	int race =nodes[community.individuals[i]].race;
	int income =nodes[community.individuals[i]].income;
	int ethnicity =nodes[community.individuals[i]].ethnicity;
	int gender=	nodes[community.individuals[i]].gender;


	if (infection_status == Progression::susceptible){
	  susceptible += 1;
	  if(hd_area_resident){
		hd_area_susceptible += 1;
	  }
	}
	if (infection_status == Progression::exposed) {
	  exposed +=1;
	  if(hd_area_resident){
        hd_area_exposed += 1;
      }
	  	 	if(age_index==0){//----------------age group----
				  infected_age_group_1+=1;
	  			}
			if(age_index==1){
				  infected_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  infected_age_group_3+=1;
	  			}
			if(age_index==3){
				  infected_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  infected_age_group_5+=1;
	  			}
			if(age_index==5){
				  infected_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  infected_age_group_7+=1;
	  			}
			if(age_index==7){
				  infected_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  infected_age_group_9+=1;
	  			}
			if(age_index==9){
				  infected_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  infected_race_group_1+=1;
	  			}
	 		 if(race==2){
				  infected_race_group_2+=1;
	  			}
			if(race==3){
				  infected_race_group_3+=1;
	  			}
	 		 if(race==4){
				  infected_race_group_4+=1;
	  			}
			if(race==5){
				  infected_race_group_5+=1;
	  			}
	 		 if(race==6){
				  infected_race_group_6+=1;
	  			}
			if(race==7){
				  infected_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  infected_income_group_1+=1;
	  			}
	 		 if(income==2){
				  infected_income_group_2+=1;
	  			}
			if(income==3){
				  infected_income_group_3+=1;
	  			}
	 		 if(income==4){
				  infected_income_group_4+=1;
	  			}
			if(income==5){
				  infected_income_group_5+=1;
	  			}
	 		 if(income==6){
				  infected_income_group_6+=1;
	  			}
			if(income==7){
				  infected_income_group_7+=1;
	  			}
	 		 if(income==8){
				  infected_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  infected_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  infected_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  infected_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  infected_gender_group_2+=1;
             }
	}
	if (infection_status == Progression::infective) {
	  infective +=1;
	  if(hd_area_resident){
        hd_area_infective += 1;
      }

	  	 	if(age_index==0){
				  infected_age_group_1+=1;
	  			}
			if(age_index==1){
				  infected_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  infected_age_group_3+=1;
	  			}
			if(age_index==3){
				  infected_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  infected_age_group_5+=1;
	  			}
			if(age_index==5){
				  infected_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  infected_age_group_7+=1;
	  			}
			if(age_index==7){
				  infected_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  infected_age_group_9+=1;
	  			}
			if(age_index==9){
				  infected_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  infected_race_group_1+=1;
	  			}
	 		 if(race==2){
				  infected_race_group_2+=1;
	  			}
			if(race==3){
				  infected_race_group_3+=1;
	  			}
	 		 if(race==4){
				  infected_race_group_4+=1;
	  			}
			if(race==5){
				  infected_race_group_5+=1;
	  			}
	 		 if(race==6){
				  infected_race_group_6+=1;
	  			}
			if(race==7){
				  infected_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  infected_income_group_1+=1;
	  			}
	 		 if(income==2){
				  infected_income_group_2+=1;
	  			}
			if(income==3){
				  infected_income_group_3+=1;
	  			}
	 		 if(income==4){
				  infected_income_group_4+=1;
	  			}
			if(income==5){
				  infected_income_group_5+=1;
	  			}
	 		 if(income==6){
				  infected_income_group_6+=1;
	  			}
			if(income==7){
				  infected_income_group_7+=1;
	  			}
	 		 if(income==8){
				  infected_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  infected_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  infected_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  infected_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  infected_gender_group_2+=1;
             }
	}
	if (infection_status == Progression::symptomatic) {
	  symptomatic += 1;
      if(hd_area_resident){
        hd_area_symptomatic += 1;
      }

	  	 	if(age_index==0){
				  infected_age_group_1+=1;
	  			}
			if(age_index==1){
				  infected_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  infected_age_group_3+=1;
	  			}
			if(age_index==3){
				  infected_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  infected_age_group_5+=1;
	  			}
			if(age_index==5){
				  infected_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  infected_age_group_7+=1;
	  			}
			if(age_index==7){
				  infected_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  infected_age_group_9+=1;
	  			}
			if(age_index==9){
				  infected_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  infected_race_group_1+=1;
	  			}
	 		 if(race==2){
				  infected_race_group_2+=1;
	  			}
			if(race==3){
				  infected_race_group_3+=1;
	  			}
	 		 if(race==4){
				  infected_race_group_4+=1;
	  			}
			if(race==5){
				  infected_race_group_5+=1;
	  			}
	 		 if(race==6){
				  infected_race_group_6+=1;
	  			}
			if(race==7){
				  infected_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  infected_income_group_1+=1;
	  			}
	 		 if(income==2){
				  infected_income_group_2+=1;
	  			}
			if(income==3){
				  infected_income_group_3+=1;
	  			}
	 		 if(income==4){
				  infected_income_group_4+=1;
	  			}
			if(income==5){
				  infected_income_group_5+=1;
	  			}
	 		 if(income==6){
				  infected_income_group_6+=1;
	  			}
			if(income==7){
				  infected_income_group_7+=1;
	  			}
	 		 if(income==8){
				  infected_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  infected_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  infected_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  infected_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  infected_gender_group_2+=1;
             }
	}
	if (infection_status == Progression::recovered) {
	  recovered += 1;
	  if(hd_area_resident){
        hd_area_recovered += 1;
      }
	  auto state_before_recovery
		= nodes[community.individuals[i]].state_before_recovery;
	  switch(state_before_recovery){
	  case Progression::infective:
		recovered_from_infective += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_infective += 1;
		}
		break;
	  case Progression::symptomatic:
		recovered_from_symptomatic += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_symptomatic += 1;
		}
		break;
	  case Progression::hospitalised:
		recovered_from_hospitalised += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_hospitalised += 1;
		}
		break;
	  case Progression::critical:
		recovered_from_critical += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_critical += 1;
		}
		break;
	  case Progression::vaccinated1:
		recovered_from_vaccination += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_vaccination += 1;
		}
	  case Progression::vaccinated2:
		recovered_from_vaccination += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_vaccination += 1;
		}
	   case Progression::waning:
		recovered_from_vaccination += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_vaccination += 1;
		}
	   case Progression::boosted:
		recovered_from_vaccination += 1;
		if(hd_area_resident){
		  hd_area_recovered_from_vaccination += 1;
		}
		break;
		
	  default:
		errors += 1; //errors state_before_recovery
		break;
	  }
	}
	if (infection_status == Progression::hospitalised) {
	  hospitalised += 1;
      if(hd_area_resident){
        hd_area_hospitalised += 1;
      }
	 		 if(age_index==0){
				  hospitalised_age_group_1+=1;
	  			}
			if(age_index==1){
				  hospitalised_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  hospitalised_age_group_3+=1;
	  			}
			if(age_index==3){
				  hospitalised_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  hospitalised_age_group_5+=1;
	  			}
			if(age_index==5){
				  hospitalised_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  hospitalised_age_group_7+=1;
	  			}
			if(age_index==7){
				  hospitalised_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  hospitalised_age_group_9+=1;
	  			}
			if(age_index==9){
				  hospitalised_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  hospitalised_race_group_1+=1;
	  			}
	 		 if(race==2){
				  hospitalised_race_group_2+=1;
	  			}
			if(race==3){
				  hospitalised_race_group_3+=1;
	  			}
	 		 if(race==4){
				  hospitalised_race_group_4+=1;
	  			}
			if(race==5){
				  hospitalised_race_group_5+=1;
	  			}
	 		 if(race==6){
				  hospitalised_race_group_6+=1;
	  			}
			if(race==7){
				  hospitalised_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  hospitalised_income_group_1+=1;
	  			}
	 		 if(income==2){
				  hospitalised_income_group_2+=1;
	  			}
			if(income==3){
				  hospitalised_income_group_3+=1;
	  			}
	 		 if(income==4){
				  hospitalised_income_group_4+=1;
	  			}
			if(income==5){
				  hospitalised_income_group_5+=1;
	  			}
	 		 if(income==6){
				  hospitalised_income_group_6+=1;
	  			}
			if(income==7){
				  hospitalised_income_group_7+=1;
	  			}
	 		 if(income==8){
				  hospitalised_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  hospitalised_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  hospitalised_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  hospitalised_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  hospitalised_gender_group_2+=1;
             }
//-----------infected because hospitalized------
	  	 	if(age_index==0){
				  infected_age_group_1+=1;
	  			}
			if(age_index==1){
				  infected_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  infected_age_group_3+=1;
	  			}
			if(age_index==3){
				  infected_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  infected_age_group_5+=1;
	  			}
			if(age_index==5){
				  infected_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  infected_age_group_7+=1;
	  			}
			if(age_index==7){
				  infected_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  infected_age_group_9+=1;
	  			}
			if(age_index==9){
				  infected_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  infected_race_group_1+=1;
	  			}
	 		 if(race==2){
				  infected_race_group_2+=1;
	  			}
			if(race==3){
				  infected_race_group_3+=1;
	  			}
	 		 if(race==4){
				  infected_race_group_4+=1;
	  			}
			if(race==5){
				  infected_race_group_5+=1;
	  			}
	 		 if(race==6){
				  infected_race_group_6+=1;
	  			}
			if(race==7){
				  infected_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  infected_income_group_1+=1;
	  			}
	 		 if(income==2){
				  infected_income_group_2+=1;
	  			}
			if(income==3){
				  infected_income_group_3+=1;
	  			}
	 		 if(income==4){
				  infected_income_group_4+=1;
	  			}
			if(income==5){
				  infected_income_group_5+=1;
	  			}
	 		 if(income==6){
				  infected_income_group_6+=1;
	  			}
			if(income==7){
				  infected_income_group_7+=1;
	  			}
	 		 if(income==8){
				  infected_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  infected_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  infected_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  infected_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  infected_gender_group_2+=1;
             }
	}
	if (infection_status == Progression::critical) {
	  critical += 1;
      if(hd_area_resident){
        hd_area_critical += 1;
      }
	  	 	if(age_index==0){
				  infected_age_group_1+=1;
	  			}
			if(age_index==1){
				  infected_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  infected_age_group_3+=1;
	  			}
			if(age_index==3){
				  infected_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  infected_age_group_5+=1;
	  			}
			if(age_index==5){
				  infected_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  infected_age_group_7+=1;
	  			}
			if(age_index==7){
				  infected_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  infected_age_group_9+=1;
	  			}
			if(age_index==9){
				  infected_age_group_10+=1;
	  			}
		//----------------race group----

			if(race==1){
				  infected_race_group_1+=1;
	  			}
	 		 if(race==2){
				  infected_race_group_2+=1;
	  			}
			if(race==3){
				  infected_race_group_3+=1;
	  			}
	 		 if(race==4){
				  infected_race_group_4+=1;
	  			}
			if(race==5){
				  infected_race_group_5+=1;
	  			}
	 		 if(race==6){
				  infected_race_group_6+=1;
	  			}
			if(race==7){
				  infected_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  infected_income_group_1+=1;
	  			}
	 		 if(income==2){
				  infected_income_group_2+=1;
	  			}
			if(income==3){
				  infected_income_group_3+=1;
	  			}
	 		 if(income==4){
				  infected_income_group_4+=1;
	  			}
			if(income==5){
				  infected_income_group_5+=1;
	  			}
	 		 if(income==6){
				  infected_income_group_6+=1;
	  			}
			if(income==7){
				  infected_income_group_7+=1;
	  			}
	 		 if(income==8){
				  infected_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  infected_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  infected_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  infected_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  infected_gender_group_2+=1;
             }					  	  	  
	}
	if (infection_status == Progression::dead) {
	  dead += 1;
	  if(hd_area_resident){
        hd_area_dead += 1;
      }
	 		 if(age_index==0){
				  dead_age_group_1+=1;
	  			}
			if(age_index==1){
				  dead_age_group_2+=1;
	  			}
	 		 if(age_index==2){
				  dead_age_group_3+=1;
	  			}
			if(age_index==3){
				  dead_age_group_4+=1;
	  			}
	 		 if(age_index==4){
				  dead_age_group_5+=1;
	  			}
			if(age_index==5){
				  dead_age_group_6+=1;
	  			}
	 		 if(age_index==6){
				  dead_age_group_7+=1;
	  			}
			if(age_index==7){
				  dead_age_group_8+=1;
	  			}
	 		 if(age_index==8){
				  dead_age_group_9+=1;
	  			}
			if(age_index==9){
				  dead_age_group_10+=1;
	  			}

		//----------------race group----

			if(race==1){
				  dead_race_group_1+=1;
	  			}
	 		 if(race==2){
				  dead_race_group_2+=1;
	  			}
			if(race==3){
				  dead_race_group_3+=1;
	  			}
	 		 if(race==4){
				  dead_race_group_4+=1;
	  			}
			if(race==5){
				  dead_race_group_5+=1;
	  			}
	 		 if(race==6){
				  dead_race_group_6+=1;
	  			}
			if(race==7){
				  dead_race_group_7+=1;
	  			}

		//----------------income group----

			if(income==1){
				  dead_income_group_1+=1;
	  			}
	 		 if(income==2){
				  dead_income_group_2+=1;
	  			}
			if(income==3){
				  dead_income_group_3+=1;
	  			}
	 		 if(income==4){
				  dead_income_group_4+=1;
	  			}
			if(income==5){
				  dead_income_group_5+=1;
	  			}
	 		 if(income==6){
				  dead_income_group_6+=1;
	  			}
			if(income==7){
				  dead_income_group_7+=1;
	  			}
	 		 if(income==8){
				  dead_income_group_8+=1;
	  			}


		//----------------ethnicity group----

			if(ethnicity==1){
				  dead_ethnicity_group_1+=1;
	  			}
	 		 if(ethnicity==2){
				  dead_ethnicity_group_2+=1;
	  			}

		//----------------gender group----

			if(gender==1){
				  dead_gender_group_1+=1;
	  			}
	 		 if(gender==2){
				  dead_gender_group_2+=1;
             }				
//--------------counted in infected because dead------------------//

	  	//  	if(age_index==0){
		// 		  infected_age_group_1+=1;
	  	// 		}
		// 	if(age_index==1){
		// 		  infected_age_group_2+=1;
	  	// 		}
	 	// 	 if(age_index==2){
		// 		  infected_age_group_3+=1;
	  	// 		}
		// 	if(age_index==3){
		// 		  infected_age_group_4+=1;
	  	// 		}
	 	// 	 if(age_index==4){
		// 		  infected_age_group_5+=1;
	  	// 		}
		// 	if(age_index==5){
		// 		  infected_age_group_6+=1;
	  	// 		}
	 	// 	 if(age_index==6){
		// 		  infected_age_group_7+=1;
	  	// 		}
		// 	if(age_index==7){
		// 		  infected_age_group_8+=1;
	  	// 		}
	 	// 	 if(age_index==8){
		// 		  infected_age_group_9+=1;
	  	// 		}
		// 	if(age_index==9){
		// 		  infected_age_group_10+=1;
	  	// 		}
		// //----------------race group----

		// 	if(race==1){
		// 		  infected_race_group_1+=1;
	  	// 		}
	 	// 	 if(race==2){
		// 		  infected_race_group_2+=1;
	  	// 		}
		// 	if(race==3){
		// 		  infected_race_group_3+=1;
	  	// 		}
	 	// 	 if(race==4){
		// 		  infected_race_group_4+=1;
	  	// 		}
		// 	if(race==5){
		// 		  infected_race_group_5+=1;
	  	// 		}
	 	// 	 if(race==6){
		// 		  infected_race_group_6+=1;
	  	// 		}
		// 	if(race==7){
		// 		  infected_race_group_7+=1;
	  	// 		}

		//----------------income group----

		// 	if(income==1){
		// 		  infected_income_group_1+=1;
	  	// 		}
	 	// 	 if(income==2){
		// 		  infected_income_group_2+=1;
	  	// 		}
		// 	if(income==3){
		// 		  infected_income_group_3+=1;
	  	// 		}
	 	// 	 if(income==4){
		// 		  infected_income_group_4+=1;
	  	// 		}
		// 	if(income==5){
		// 		  infected_income_group_5+=1;
	  	// 		}
	 	// 	 if(income==6){
		// 		  infected_income_group_6+=1;
	  	// 		}
		// 	if(income==7){
		// 		  infected_income_group_7+=1;
	  	// 		}
	 	// 	 if(income==8){
		// 		  infected_income_group_8+=1;
	  	// 		}


		// //----------------ethnicity group----

		// 	if(ethnicity==1){
		// 		  infected_ethnicity_group_1+=1;
	  	// 		}
	 	// 	 if(ethnicity==2){
		// 		  infected_ethnicity_group_2+=1;
	  	// 		}

		// //----------------gender group----

		// 	if(gender==1){
		// 		  infected_gender_group_1+=1;
	  	// 		}
	 	// 	 if(gender==2){
		// 		  infected_gender_group_2+=1;
        //      }
	}
//----Counting vaccinated individuals in a community----------//
	auto state_before_recovery
		= nodes[community.individuals[i]].state_before_recovery;

	if((state_before_recovery == Progression::vaccinated1)){
		vaccinated1 += 1;
	  }

	if((state_before_recovery == Progression::vaccinated2)){
		vaccinated2 += 1;
	  }
	if((state_before_recovery == Progression::waning)){
		waning += 1;
	  }
	if((state_before_recovery == Progression::boosted)){
		boosted += 1;
	  }	  	  	  

	if((state_before_recovery == Progression::boosted2)){
		boosted2 += 1;
	  }	    	  	  
	
  }
  if(errors){
	cerr << "erroneous state_before_recovery found\n";
	assert(false);
  }
  
    affected = exposed + infective + symptomatic
	+ hospitalised + critical
	+ recovered + dead;

  hd_area_affected = hd_area_exposed + hd_area_infective + hd_area_symptomatic
	+ hd_area_hospitalised + hd_area_critical
	+ hd_area_recovered + hd_area_dead;

  
  casualty_stats stat;
  stat.affected = affected;
  stat.hd_area_affected = hd_area_affected;
  stat.susceptible = susceptible;
  stat.unvaccinated=unvaccinated;
  stat.vaccinated1=vaccinated1;
  stat.vaccinated2=vaccinated2;
  stat.waning=waning;
  stat.boosted=boosted;  
  stat.boosted2=boosted2;  
  stat.hd_area_susceptible = hd_area_susceptible;
  stat.exposed = exposed;
  stat.hd_area_exposed = hd_area_exposed;
  stat.infective = infective;
  stat.hd_area_infective = hd_area_infective;
  stat.symptomatic = symptomatic;
  stat.hd_area_symptomatic = hd_area_symptomatic;
  stat.hospitalised = hospitalised;
  stat.hd_area_hospitalised = hd_area_hospitalised;
  stat.critical = critical;
  stat.hd_area_critical = hd_area_critical;
  stat.dead = dead;
  stat.hd_area_dead = hd_area_dead;
  stat.recovered = recovered;

  stat.infected_age_group_1= infected_age_group_1;
  stat.infected_age_group_2=infected_age_group_2;  
  stat.infected_age_group_3= infected_age_group_3;
  stat.infected_age_group_4=infected_age_group_4;  
  stat.infected_age_group_5= infected_age_group_5;
  stat.infected_age_group_6=infected_age_group_6;//Age group shakir addition,exposed + infective + symptomatic+ hospitalised + critical
  stat.infected_age_group_7= infected_age_group_7;
  stat.infected_age_group_8=infected_age_group_8;
  stat.infected_age_group_9= infected_age_group_9;
  stat.infected_age_group_10=infected_age_group_10;

  stat.hospitalised_age_group_1= hospitalised_age_group_1;
  stat.hospitalised_age_group_2=hospitalised_age_group_2;  
  stat.hospitalised_age_group_3= hospitalised_age_group_3;
  stat.hospitalised_age_group_4=hospitalised_age_group_4;  
  stat.hospitalised_age_group_5= hospitalised_age_group_5;
  stat.hospitalised_age_group_6=hospitalised_age_group_6;//Age group shakir addition
  stat.hospitalised_age_group_7= hospitalised_age_group_7;
  stat.hospitalised_age_group_8=hospitalised_age_group_8;
  stat.hospitalised_age_group_9= hospitalised_age_group_9;
  stat.hospitalised_age_group_10=hospitalised_age_group_10;

  stat.dead_age_group_1= dead_age_group_1;
  stat.dead_age_group_2=dead_age_group_2;  
  stat.dead_age_group_3= dead_age_group_3;
  stat.dead_age_group_4=dead_age_group_4;  
  stat.dead_age_group_5= dead_age_group_5;
  stat.dead_age_group_6=dead_age_group_6;//Age group shakir addition
  stat.dead_age_group_7= dead_age_group_7;
  stat.dead_age_group_8=dead_age_group_8;  
  stat.dead_age_group_9= dead_age_group_9;
  stat.dead_age_group_10=dead_age_group_10;  

  stat.infected_race_group_1= infected_race_group_1;
  stat.infected_race_group_2=infected_race_group_2;  
  stat.infected_race_group_3= infected_race_group_3;
  stat.infected_race_group_4=infected_race_group_4;  
  stat.infected_race_group_5= infected_race_group_5;
  stat.infected_race_group_6=infected_race_group_6;//race group shakir addition,exposed + infective + symptomatic+ hospitalised + critical
  stat.infected_race_group_7= infected_race_group_7;

  stat.hospitalised_race_group_1= hospitalised_race_group_1;
  stat.hospitalised_race_group_2=hospitalised_race_group_2;  
  stat.hospitalised_race_group_3= hospitalised_race_group_3;
  stat.hospitalised_race_group_4=hospitalised_race_group_4;  
  stat.hospitalised_race_group_5= hospitalised_race_group_5;
  stat.hospitalised_race_group_6=hospitalised_race_group_6;//race group shakir addition
  stat.hospitalised_race_group_7= hospitalised_race_group_7;

  stat.dead_race_group_1= dead_race_group_1;
  stat.dead_race_group_2=dead_race_group_2;  
  stat.dead_race_group_3= dead_race_group_3;
  stat.dead_race_group_4=dead_race_group_4;  
  stat.dead_race_group_5= dead_race_group_5;
  stat.dead_race_group_6=dead_race_group_6;//race group shakir addition
  stat.dead_race_group_7= dead_race_group_7;
 
 
  stat.infected_income_group_1= infected_income_group_1;
  stat.infected_income_group_2=infected_income_group_2;  
  stat.infected_income_group_3= infected_income_group_3;
  stat.infected_income_group_4=infected_income_group_4;  
  stat.infected_income_group_5= infected_income_group_5;
  stat.infected_income_group_6=infected_income_group_6;//income group shakir addition,exposed + infective + symptomatic+ hospitalised + critical
  stat.infected_income_group_7= infected_income_group_7;
  stat.infected_income_group_8=infected_income_group_8;

  stat.hospitalised_income_group_1= hospitalised_income_group_1;
  stat.hospitalised_income_group_2=hospitalised_income_group_2;  
  stat.hospitalised_income_group_3= hospitalised_income_group_3;
  stat.hospitalised_income_group_4=hospitalised_income_group_4;  
  stat.hospitalised_income_group_5= hospitalised_income_group_5;
  stat.hospitalised_income_group_6=hospitalised_income_group_6;//income group shakir addition
  stat.hospitalised_income_group_7= hospitalised_income_group_7;
  stat.hospitalised_income_group_8=hospitalised_income_group_8;

  stat.dead_income_group_1= dead_income_group_1;
  stat.dead_income_group_2=dead_income_group_2;  
  stat.dead_income_group_3= dead_income_group_3;
  stat.dead_income_group_4=dead_income_group_4;  
  stat.dead_income_group_5= dead_income_group_5;
  stat.dead_income_group_6=dead_income_group_6;//income group shakir addition
  stat.dead_income_group_7= dead_income_group_7;
  stat.dead_income_group_8=dead_income_group_8;  
 
  stat.infected_ethnicity_group_1= infected_ethnicity_group_1;
  stat.infected_ethnicity_group_2=infected_ethnicity_group_2;  

  stat.hospitalised_ethnicity_group_1= hospitalised_ethnicity_group_1;
  stat.hospitalised_ethnicity_group_2=hospitalised_ethnicity_group_2;  

  stat.dead_ethnicity_group_1= dead_ethnicity_group_1;
  stat.dead_ethnicity_group_2=dead_ethnicity_group_2;  

  stat.infected_gender_group_1= infected_gender_group_1;
  stat.infected_gender_group_2=infected_gender_group_2;  

  stat.hospitalised_gender_group_1= hospitalised_gender_group_1;
  stat.hospitalised_gender_group_2=hospitalised_gender_group_2;  

  stat.dead_gender_group_1= dead_gender_group_1;
  stat.dead_gender_group_2=dead_gender_group_2;  

  stat.hd_area_recovered = hd_area_recovered;
  stat.recovered_from_infective = recovered_from_infective;
  stat.recovered_from_symptomatic = recovered_from_symptomatic;
  stat.recovered_from_hospitalised = recovered_from_hospitalised;
  stat.recovered_from_critical = recovered_from_critical;
  stat.hd_area_recovered_from_infective = hd_area_recovered_from_infective;
  stat.hd_area_recovered_from_symptomatic = hd_area_recovered_from_symptomatic;
  stat.hd_area_recovered_from_hospitalised = hd_area_recovered_from_hospitalised;
  stat.hd_area_recovered_from_critical = hd_area_recovered_from_critical;

  return stat;
  // Populate it afterwards...
}

void update_grid_cell_statistics(matrix<nbr_cell>& nbr_cells,
								 vector<house>& homes,
								 vector<agent>& nodes,
								 const double locked_neighborhood_leakage,
								 const double locked_neighborhood_threshold) {
  for(auto& nbr_cell_row: nbr_cells){
	for(auto& nbr_cell: nbr_cell_row){

	  const auto SIZE = nbr_cell.houses_list.size();
	  count_type num_active_hospitalisations = 0;

#pragma omp parallel for shared(homes, nodes, nbr_cell) \
  reduction(+: num_active_hospitalisations)
	  for(count_type i = 0; i < SIZE; ++i){
		for(const auto individual_index: homes[nbr_cell.houses_list[i]].individuals){
		  if(nodes[individual_index].infection_status
			 == Progression::hospitalised){
			++num_active_hospitalisations;
		  }
		}
	  }
	  nbr_cell.num_active_hospitalisations = num_active_hospitalisations;
	  nbr_cell.access_factor = interpolate(1.0, locked_neighborhood_leakage,
										   double(nbr_cell.num_active_hospitalisations)/double(nbr_cell.population),
										   locked_neighborhood_threshold);

#pragma omp parallel for shared(homes, nodes, nbr_cell)
	  for(count_type i = 0; i < SIZE; ++i){
		homes[nbr_cell.houses_list[i]].neighborhood_access_factor
		  = nbr_cell.access_factor;
		for(const auto individual_index: homes[nbr_cell.houses_list[i]].individuals){
		  nodes[individual_index].neighborhood_access_factor
			= nbr_cell.access_factor;
		}
	  }

	}
  }
}




//Addtional functions added by shakir. First is for reducing betas. The second is vaccine implementation via functions---

//updates the BETAS to reflect control measures--->social distancing and increased mask wearing--Shakir
void update_scales_Intv(vector<house>& homes, vector<workplace>& workplaces,vector<community>& communities, vector<house>& houses, vector<vector<nbr_cell>>& nbr_cells){
	if(GLOBAL.MEASURES==1)
	{
				for (count_type w = 0; w < homes.size(); ++w){
					homes[w].scale = homes[w].scale*(1-GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < workplaces.size(); ++w) {
					workplaces[w].scale=workplaces[w].scale*(1-GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < communities.size(); ++w) {
					communities[w].scale =communities[w].scale *(1-GLOBAL.REDUCTION_FACTOR);
				}

				for(count_type i=0; i<houses.size(); ++i){
				houses[i].random_households.scale=houses[i].random_households.scale*(1-GLOBAL.REDUCTION_FACTOR);
				}
				for(count_type i=0; i<nbr_cells.size(); ++i){
					for(count_type j=0; j<nbr_cells[i].size(); ++j){
					nbr_cells[i][j].scale=nbr_cells[i][j].scale*(1-GLOBAL.REDUCTION_FACTOR);
					}
				}
	}

	if(GLOBAL.MEASURES==2)
	{
				for (count_type w = 0; w < homes.size(); ++w){
					homes[w].scale = homes[w].scale*(1+GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < workplaces.size(); ++w) {
					workplaces[w].scale=workplaces[w].scale*(1+GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < communities.size(); ++w) {
					communities[w].scale =communities[w].scale *(1+GLOBAL.REDUCTION_FACTOR);
				}

				for(count_type i=0; i<houses.size(); ++i){
				houses[i].random_households.scale=houses[i].random_households.scale*(1+GLOBAL.REDUCTION_FACTOR);
				}
				for(count_type i=0; i<nbr_cells.size(); ++i){
					for(count_type j=0; j<nbr_cells[i].size(); ++j){
					nbr_cells[i][j].scale=nbr_cells[i][j].scale*(1+GLOBAL.REDUCTION_FACTOR);
					}
				}
	}
	if(GLOBAL.MEASURES==3)
	{
				for (count_type w=0; w < workplaces.size(); ++w) {
					workplaces[w].scale=workplaces[w].scale*(1-GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < communities.size(); ++w) {
					communities[w].scale =communities[w].scale *(1-GLOBAL.REDUCTION_FACTOR);
				}

				for(count_type i=0; i<houses.size(); ++i){
				houses[i].random_households.scale=houses[i].random_households.scale*(1-GLOBAL.REDUCTION_FACTOR);
				}
				for(count_type i=0; i<nbr_cells.size(); ++i){
					for(count_type j=0; j<nbr_cells[i].size(); ++j){
					nbr_cells[i][j].scale=nbr_cells[i][j].scale*(1-GLOBAL.REDUCTION_FACTOR);
					}
				}
	}	

	if(GLOBAL.MEASURES==4)
	{
				for (count_type w=0; w < workplaces.size(); ++w) {
					workplaces[w].scale=workplaces[w].scale*(1+GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < communities.size(); ++w) {
					communities[w].scale =communities[w].scale *(1+GLOBAL.REDUCTION_FACTOR);
				}

				for(count_type i=0; i<houses.size(); ++i){
				houses[i].random_households.scale=houses[i].random_households.scale*(1+GLOBAL.REDUCTION_FACTOR);
				}
				for(count_type i=0; i<nbr_cells.size(); ++i){
					for(count_type j=0; j<nbr_cells[i].size(); ++j){
					nbr_cells[i][j].scale=nbr_cells[i][j].scale*(1+GLOBAL.REDUCTION_FACTOR);
					}
				}
	}
		if(GLOBAL.MEASURES==5)
	{
				for (count_type w=0; w < workplaces.size(); ++w) {
					workplaces[w].scale=workplaces[w].scale*(1-2*GLOBAL.REDUCTION_FACTOR);
				}
				for (count_type w=0; w < communities.size(); ++w) {
					communities[w].scale =communities[w].scale *(1-2*GLOBAL.REDUCTION_FACTOR);
				}

				for(count_type i=0; i<houses.size(); ++i){
				houses[i].random_households.scale=houses[i].random_households.scale*(1-2*GLOBAL.REDUCTION_FACTOR);
				}
				for(count_type i=0; i<nbr_cells.size(); ++i){
					for(count_type j=0; j<nbr_cells[i].size(); ++j){
					nbr_cells[i][j].scale=nbr_cells[i][j].scale*(1-2*GLOBAL.REDUCTION_FACTOR);
					}
				}
	}			
}

//--------------------Interventions to reduce betas is over-----------//


void vaccinate_firstdose(vector<agent>& nodes, vector<count_type> new_vaccinated1_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated1_candidates_list_size = new_vaccinated1_candidates.size();
					  if (new_vaccinated1_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated1_candidates.begin(), new_vaccinated1_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated1 = vaccFn;
				  count_type num = std::min(new_vaccinated1_candidates_list_size, num_of_individuals_to_be_vaccinated1);
			       for(count_type a = 0; a < num; ++a){
				    nodes[new_vaccinated1_candidates[a]].vaccinated1 = true;
				    nodes[new_vaccinated1_candidates[a]].new_vaccinated1 = true;

				    nodes[new_vaccinated1_candidates[a]].time_at_vaccine1=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS1);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::vaccinated1;
				  }
				}
}
void vaccinate_second_dose(vector<agent>& nodes, vector<count_type> new_vaccinated2_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated2_candidates_list_size = new_vaccinated2_candidates.size();
					  if (new_vaccinated2_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated2_candidates.begin(), new_vaccinated2_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated2= vaccFn;
				  count_type num = std::min(new_vaccinated2_candidates_list_size, num_of_individuals_to_be_vaccinated2);
			       for(count_type a = 0; a < num; ++a){
				    nodes[new_vaccinated2_candidates[a]].vaccinated2 = true;
				    nodes[new_vaccinated2_candidates[a]].new_vaccinated2 = true;

				    nodes[new_vaccinated2_candidates[a]].time_at_vaccine2=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS2);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated2_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated2_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated2_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated2_candidates[a]].state_before_recovery = Progression::vaccinated2;
				  }
				}
}

void vaccinate_waning_candidates(vector<agent>& nodes, vector<count_type> new_vaccinated1_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated1_candidates_list_size = new_vaccinated1_candidates.size();
					  if (new_vaccinated1_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated1_candidates.begin(), new_vaccinated1_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated1 = vaccFn;
				  count_type num = std::min(new_vaccinated1_candidates_list_size, num_of_individuals_to_be_vaccinated1);
			       for(count_type a = 0; a < num; ++a){
					if(double(time_step)-nodes[new_vaccinated1_candidates[a]].time_at_vaccine2>=5.0*30*GLOBAL.SIM_STEPS_PER_DAY){
				    nodes[new_vaccinated1_candidates[a]].waning = true;
				    nodes[new_vaccinated1_candidates[a]].new_waning = true;

				    nodes[new_vaccinated1_candidates[a]].time_at_waning=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS_WANING);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::waning;
				  }
				  else if((!vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){
					 nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::vaccinated2;
					}
				}
				   }
}


void vaccinate_booster_dose(vector<agent>& nodes, vector<count_type> new_vaccinated1_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated1_candidates_list_size = new_vaccinated1_candidates.size();
					  if (new_vaccinated1_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated1_candidates.begin(), new_vaccinated1_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated1 = vaccFn;
				  count_type num = std::min(new_vaccinated1_candidates_list_size, num_of_individuals_to_be_vaccinated1);
			       for(count_type a = 0; a < num; ++a){
				   if(double(time_step)-nodes[new_vaccinated1_candidates[a]].time_at_vaccine2>=6.0*30*GLOBAL.SIM_STEPS_PER_DAY){

				    nodes[new_vaccinated1_candidates[a]].boosted = true;
				    nodes[new_vaccinated1_candidates[a]].new_boosted = true;

				    nodes[new_vaccinated1_candidates[a]].time_at_boosted=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS_BOOSTED);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::boosted;
				  }
				   }
				}
}

void vaccinate_booster2_dose(vector<agent>& nodes, vector<count_type> new_vaccinated1_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated1_candidates_list_size = new_vaccinated1_candidates.size();
					  if (new_vaccinated1_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated1_candidates.begin(), new_vaccinated1_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated1 = vaccFn;
				  count_type num = std::min(new_vaccinated1_candidates_list_size, num_of_individuals_to_be_vaccinated1);
			       for(count_type a = 0; a < num; ++a){
				   if(double(time_step)-nodes[new_vaccinated1_candidates[a]].boosted>=6.0*30*GLOBAL.SIM_STEPS_PER_DAY){

				    nodes[new_vaccinated1_candidates[a]].boosted2 = true;
				    nodes[new_vaccinated1_candidates[a]].new_boosted2 = true;

				    nodes[new_vaccinated1_candidates[a]].time_at_boosted2=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS_BOOSTED2);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::boosted2;
				  }
				   }
				}
}
void vaccinate_waning2_candidates(vector<agent>& nodes, vector<count_type> new_vaccinated1_candidates, count_type vaccFn,count_type time_step){
					  count_type new_vaccinated1_candidates_list_size = new_vaccinated1_candidates.size();
					  if (new_vaccinated1_candidates_list_size > vaccFn){
				     //Randomly permute the list of candidates
				     std::shuffle(new_vaccinated1_candidates.begin(), new_vaccinated1_candidates.end(), GENERATOR);	  
				   }
				   bool vaccination_works = false;
				  count_type num_of_individuals_to_be_vaccinated1 = vaccFn;
				  count_type num = std::min(new_vaccinated1_candidates_list_size, num_of_individuals_to_be_vaccinated1);
			       for(count_type a = 0; a < num; ++a){
					if(double(time_step)-nodes[new_vaccinated1_candidates[a]].boosted>=5.0*30*GLOBAL.SIM_STEPS_PER_DAY){

				    nodes[new_vaccinated1_candidates[a]].waning2 = true;
				    nodes[new_vaccinated1_candidates[a]].new_waning2 = true;

				    nodes[new_vaccinated1_candidates[a]].time_at_waning2=time_step;
					 vaccination_works = bernoulli(GLOBAL.VACCINATION_EFFECTIVENESS_WANING2);//1;//drand48();
				 if ((vaccination_works)&&(nodes[new_vaccinated1_candidates[a]].infection_status==Progression::susceptible ||nodes[new_vaccinated1_candidates[a]].infection_status==Progression::recovered)){ 
				     nodes[new_vaccinated1_candidates[a]].infection_status = Progression::recovered;
				     nodes[new_vaccinated1_candidates[a]].state_before_recovery = Progression::waning2;
				  }
				}
				   }
}