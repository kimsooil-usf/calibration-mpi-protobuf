//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#include "models.h"
#include "intervention_primitives.h"
#include <algorithm>
#include <vector>

using std::vector;
using std::min;

void set_kappa_base_node(agent& node, double community_factor, const int cur_time){
  //set the basic kappa values for this node according to current time
  node.quarantined=false;
  node.kappa_T = kappa_T(node, cur_time);
  node.kappa_H = 1.0;
  node.kappa_H_incoming = 1.0;
  node.kappa_W = 1.0;
  node.kappa_W_incoming = 1.0;
  if(node.compliant){
    node.kappa_C = community_factor;
    node.kappa_C_incoming = community_factor;
  }else{
    node.kappa_C = 1.0;
    node.kappa_C_incoming = 1.0;
  }
}

void set_kappa_lockdown_node(agent& node, const int cur_time, const intervention_params intv_params){
  node.kappa_T = kappa_T(node, cur_time);
  node.quarantined = true; //lockdown implies quarantined
  if(node.workplace_type==WorkplaceType::office){
    node.kappa_W = intv_params.lockdown_kappas_compliant.kappa_W;
    node.kappa_W_incoming = intv_params.lockdown_kappas_compliant.kappa_W_incoming;
  }else{
	//Schools and colleges are assumed closed in all lockdowns
    node.kappa_W = 0.0;
    node.kappa_W_incoming = 0.0;
  }

  if(node.compliant){
    node.kappa_H = intv_params.lockdown_kappas_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_compliant.kappa_C;
    node.kappa_C_incoming = intv_params.lockdown_kappas_compliant.kappa_C_incoming;
  }else{
    node.kappa_H = intv_params.lockdown_kappas_non_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_non_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_non_compliant.kappa_C;
    node.kappa_C_incoming = intv_params.lockdown_kappas_non_compliant.kappa_C_incoming;
  }
}


void set_kappa_lockdown_hillsborugh_node(agent& node, const int cur_time, const intervention_hillsborough_params intv_params){
  node.kappa_T = kappa_T(node, cur_time);
  node.quarantined = true; //lockdown implies quarantined
  if(cur_time<GLOBAL.FIRST_PERIOD*GLOBAL.SIM_STEPS_PER_DAY)
  {
  if(node.workplace_type==WorkplaceType::office){
    node.kappa_W = intv_params.lockdown_kappas_compliant.kappa_W/2.5;
    node.kappa_W_incoming = intv_params.lockdown_kappas_compliant.kappa_W_incoming/2.5;//10 percent of work places are banks, hospitals and grocery in the first lockdown 0.25/2.5=0.1 which is scaling the original parameter to 10 percent.
  }else{
	//Schools and colleges are assumed closed in all lockdowns
    node.kappa_W = 0.0;
    node.kappa_W_incoming = 0.0;
  }

  if(node.compliant){
    node.kappa_H = intv_params.lockdown_kappas_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_compliant.kappa_C;
    node.kappa_C_incoming = intv_params.lockdown_kappas_compliant.kappa_C_incoming;
  }else{
    node.kappa_H = intv_params.lockdown_kappas_non_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_non_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_non_compliant.kappa_C;
    node.kappa_C_incoming = intv_params.lockdown_kappas_non_compliant.kappa_C_incoming;
  }
  }//First lockdown ends--Shakir

   if(cur_time>GLOBAL.FIRST_PERIOD*GLOBAL.SIM_STEPS_PER_DAY && cur_time<(GLOBAL.FIRST_PERIOD+GLOBAL.SECOND_PERIOD)*GLOBAL.SIM_STEPS_PER_DAY)
  {
    if(node.workplace_type==WorkplaceType::office){
    node.kappa_W = intv_params.lockdown_kappas_compliant.kappa_W;
    node.kappa_W_incoming = intv_params.lockdown_kappas_compliant.kappa_W_incoming;//25 percent office are opened in the second lockdown
  }else{
	//Schools and colleges are assumed closed in all lockdowns
    node.kappa_W = 0.0;
    node.kappa_W_incoming = 0.0;
  }

  if(node.compliant){
    node.kappa_H = intv_params.lockdown_kappas_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_compliant.kappa_C;//25 percent of of the outdoor/restaurants are opepend in the second lockdown
    node.kappa_C_incoming = intv_params.lockdown_kappas_compliant.kappa_C_incoming;//25 percent of of the outdoor/restaurants are opepend in the second lockdown
  }else{
    node.kappa_H = intv_params.lockdown_kappas_non_compliant.kappa_H*2;
    node.kappa_H_incoming = intv_params.lockdown_kappas_non_compliant.kappa_H_incoming*2;
    node.kappa_C = intv_params.lockdown_kappas_non_compliant.kappa_C*2;
    node.kappa_C_incoming = intv_params.lockdown_kappas_non_compliant.kappa_C_incoming*2;
  }
  }//second lockdown ends--Shakir

   if(cur_time>(GLOBAL.FIRST_PERIOD+GLOBAL.SECOND_PERIOD)*GLOBAL.SIM_STEPS_PER_DAY && cur_time<(GLOBAL.FIRST_PERIOD+GLOBAL.SECOND_PERIOD+GLOBAL.THIRD_PERIOD)*GLOBAL.SIM_STEPS_PER_DAY)
  {
    if(node.workplace_type==WorkplaceType::office){
    node.kappa_W = intv_params.lockdown_kappas_compliant.kappa_W*2;
    node.kappa_W_incoming = intv_params.lockdown_kappas_compliant.kappa_W_incoming*2;//50 percent of (2*0.25) of the offices are opened in the third lockdown
  }else{
	//Schools and colleges are assumed closed in all lockdowns
    node.kappa_W = 0.0;
    node.kappa_W_incoming = 0.0;//schools are still closed in third lockdown
  }

  if(node.compliant){
    node.kappa_H = intv_params.lockdown_kappas_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_compliant.kappa_C;//50 percent of the outdoor/restaurants are opened in the third lockdown but complaint people do not go out
    node.kappa_C_incoming = intv_params.lockdown_kappas_compliant.kappa_C_incoming;//50 percent of of the outdoor/restaurants are opepend in the third lockdown
  }else{
    node.kappa_H = intv_params.lockdown_kappas_non_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_non_compliant.kappa_H_incoming;//50 percent of the outdoor/restaurants are opened in the third lockdown
    node.kappa_C = intv_params.lockdown_kappas_non_compliant.kappa_C*2;
    node.kappa_C_incoming = intv_params.lockdown_kappas_non_compliant.kappa_C_incoming*2;//50 percent of the outdoor/restaurants are opened in the third lockdown
  }
  }//Third lockdown ends--Shakir

     if(cur_time>(GLOBAL.FIRST_PERIOD+GLOBAL.SECOND_PERIOD+GLOBAL.THIRD_PERIOD)*GLOBAL.SIM_STEPS_PER_DAY && cur_time<(GLOBAL.FIRST_PERIOD+GLOBAL.SECOND_PERIOD+GLOBAL.THIRD_PERIOD+GLOBAL.FOURTH_PERIOD)*GLOBAL.SIM_STEPS_PER_DAY)
  {
    if(node.workplace_type==WorkplaceType::office){
    node.kappa_W = intv_params.lockdown_kappas_compliant.kappa_W*4;
    node.kappa_W_incoming = intv_params.lockdown_kappas_compliant.kappa_W_incoming*4;//100 percent of offices are opened (4*0.25) in the forth lockdown
  }else{
	//Schools and colleges are assumed closed in all lockdowns
    node.kappa_W = 0.5;//50 percent schools are opepend in fourth lockdown
    node.kappa_W_incoming = 0.5;//50 percent of schools are opened in fourth lockdown.
  }

  if(node.compliant){
    node.kappa_H = intv_params.lockdown_kappas_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_compliant.kappa_H_incoming;
    node.kappa_C = intv_params.lockdown_kappas_compliant.kappa_C;//75 percent of outdoor/restaurants are opepend in the fourth lockdown
    node.kappa_C_incoming = intv_params.lockdown_kappas_compliant.kappa_C_incoming;//75 percent of outdoor/restaurants are opepend in the fourth lockdown. But complaint people go out less.
  }else{
    node.kappa_H = intv_params.lockdown_kappas_non_compliant.kappa_H;
    node.kappa_H_incoming = intv_params.lockdown_kappas_non_compliant.kappa_H_incoming*3;
    node.kappa_C = intv_params.lockdown_kappas_non_compliant.kappa_C;
    node.kappa_C_incoming = intv_params.lockdown_kappas_non_compliant.kappa_C_incoming*3;//non compliant people go out in large numbers
  }
  }//Forth lockdown ends--Shakir
}

void modify_kappa_SDE_node(agent& node){
  if(node.age>= UPPER_AGE && node.compliant){
    node.kappa_W_incoming = min(0.25, node.kappa_W_incoming);
    node.kappa_C_incoming = min(0.25, node.kappa_C_incoming);
  }
}

void modify_kappa_SC_node(agent& node, double SC_factor){
  if (node.workplace_type==WorkplaceType::school){
    //school and colleges are closed
    node.kappa_W = min(SC_factor, node.kappa_W);
    node.kappa_W_incoming = min(SC_factor, node.kappa_W_incoming);
  }
}

void modify_kappa_OE_node(agent& node){
  if(node.workplace_type==WorkplaceType::office){
    //odd-even rule for workplaces. 50% interactions for workplaces.
    node.kappa_W = min(0.5, node.kappa_W);
    node.kappa_W_incoming = min(0.5,node.kappa_W_incoming);
  }
}

void reset_home_quarantines(vector<house>& homes){
#pragma omp parallel for default(none) shared(homes)
  for(count_type count = 0; count<homes.size(); ++count){
    homes[count].quarantined = false;
  }
}

void modify_kappa_case_isolate_node(agent& node){
  node.quarantined = true;
  node.kappa_H = min(0.75, node.kappa_H);
  node.kappa_W = min(0.0, node.kappa_W);
  node.kappa_C = min(0.1, node.kappa_C);
  node.kappa_H_incoming = min(0.75, node.kappa_H_incoming);
  node.kappa_W_incoming = min(0.0, node.kappa_W_incoming);
  node.kappa_C_incoming = min(0.1, node.kappa_C_incoming);
}

void modify_kappa_ward_containment(agent& node){
  node.quarantined = true;
  node.kappa_H = min(0.75, node.kappa_H);
  node.kappa_W = min(0.25, node.kappa_W);
  node.kappa_C = min(0.25, node.kappa_C);
  node.kappa_H_incoming = min(0.75, node.kappa_H_incoming);
  node.kappa_W_incoming = min(0.25, node.kappa_W_incoming);
  node.kappa_C_incoming = min(0.25, node.kappa_C_incoming);
}

bool should_be_isolated_node(const agent& node, const int cur_time, const int quarantine_days){
  double time_since_symptoms = cur_time
                              - (node.time_of_infection
                              + node.incubation_period
                              + node.asymptomatic_period);
  return (node.entered_symptomatic_state &&
   (time_since_symptoms > NUM_DAYS_TO_RECOG_SYMPTOMS*GLOBAL.SIM_STEPS_PER_DAY) &&
   (time_since_symptoms <= (NUM_DAYS_TO_RECOG_SYMPTOMS + quarantine_days)*GLOBAL.SIM_STEPS_PER_DAY));
}

void mark_homes_for_quarantine(const vector<agent>& nodes, vector<house>& homes, const int cur_time){
  //mark all homes for quarantine
  for (count_type count = 0; count < nodes.size(); ++count){
    if(nodes[count].compliant && should_be_isolated_node(nodes[count],cur_time, HOME_QUARANTINE_DAYS)){
       homes[nodes[count].home].quarantined = true;
     }
  }
}

void mark_neighbourhood_homes_for_quarantine(const vector<agent>& nodes, vector<house>& homes,
											 const matrix<nbr_cell>& nbr_cells, const int cur_time){
	for (count_type count = 0; count < nodes.size(); ++count){
		double time_since_hospitalised = cur_time
		- (nodes[count].time_of_infection
			+ nodes[count].incubation_period
			+ nodes[count].asymptomatic_period
			+ nodes[count].symptomatic_period);
		if(((nodes[count].entered_hospitalised_state) &&
		(time_since_hospitalised <= (HOME_QUARANTINE_DAYS)*GLOBAL.SIM_STEPS_PER_DAY)) ){
			homes[nodes[count].home].quarantined = true;		
			//TODO: Need to check if the nbr_cell's quarantined flag needs to be set.
			grid_cell my_nbr_grid_cell = homes[nodes[count].home].neighbourhood;
			nbr_cell my_nbr_cell = 	nbr_cells[my_nbr_grid_cell.cell_x][my_nbr_grid_cell.cell_y];	
			count_type num_homes_in_cell = my_nbr_cell.houses_list.size();
			for(count_type nbr_count = 0; nbr_count < num_homes_in_cell; ++nbr_count){
				count_type neighbour = my_nbr_cell.houses_list[nbr_count];
				homes[neighbour].quarantined = true;
			}		
		}
	}
}

void isolate_quarantined_residents(vector<agent>& nodes, const vector<house>& homes, const int cur_time){
  for (count_type count = 0; count < homes.size(); ++count){
    if(homes[count].quarantined){
      for(count_type resident = 0; resident < homes[count].individuals.size(); ++resident){
        modify_kappa_case_isolate_node(nodes[resident]);
      }
    }
  }
}
