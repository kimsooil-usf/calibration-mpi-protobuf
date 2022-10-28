//Copyright [2020] [Indian Institute of Science, Bangalore & Tata Institute of Fundamental Research, Mumbai]
//SPDX-License-Identifier: Apache-2.0
#ifndef INITIALIZERS_H_
#define INITIALIZERS_H_
#include "models.h"

#include <vector>


//Initialize the office attendance
void initialize_office_attendance();


std::vector<mask> init_mask();//---Mask wearing function by shakir. It reads data from masktrends.json file and uses it as input for mask wearing
std::vector<house> init_homes();
std::vector<workplace> init_workplaces();
std::vector<community> init_community();
std::vector<agent> init_nodes();
matrix<nbr_cell> init_nbr_cells();
std::vector<intervention_params> init_intervention_params();
std::vector<intervention_hillsborough_params> init_intervention_hillsborough_params();

std::vector<testing_probability> init_testing_protocol();

matrix<double> compute_community_distances(const std::vector<community>& communities);
matrix<double> compute_community_distances_fkernel(const matrix<double>& community_distances);

//Assign individuals to homes, workplace, community
void assign_individual_home_community(std::vector<agent>& nodes, std::vector<house>& homes, std::vector<workplace>& workplaces, std::vector<community>& communities);
void assign_homes_nbr_cell(const std::vector<house>& homes, matrix<nbr_cell>& nbr_cells);
void assign_individual_projects(std::vector<workplace>& workplaces, std::vector<agent>& nodes);
void assign_household_community(std::vector<community>& communities, const std::vector<agent>& nodes, std::vector<house>& homes);
void assign_household_random_community(std::vector<house>& homes, const std::vector<community>& communities);


// Compute scale factors for each home, workplace and community. Done once at the beginning.
void compute_scale_homes(std::vector<house>& homes);
void compute_scale_workplaces(std::vector<workplace>& workplaces);
void compute_scale_communities(const std::vector<agent>& nodes, std::vector<community>& communities);
void compute_scale_random_community(std::vector<house>& houses, std::vector<agent>& nodes);
void compute_scale_nbr_cells(std::vector<agent>& nodes, matrix<nbr_cell>& nbr_cells, const std::vector<house>& homes);


// Age stratification JSON read function.
svd init_home_age_interaction_matrix();
svd init_school_age_interaction_matrix();
svd init_workplace_age_interaction_matrix();
svd init_community_age_interaction_matrix();

void print_testing_protocol(const int index, const testing_probability probabilities);

#endif
