/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "./custom.h"


// macros for customizable cell parameters and scenarios
#define BIRTH_RATE_SENS 0.01
#define DEATH_RATE_SENS 0.002
#define BIRTH_RATE_RES 0.009
#define DEATH_RATE_RES 0.002

#define BIRTH_RATE_TRT_SENS 0.0005
#define DEATH_RATE_TRT_SENS 0.02
#define BIRTH_RATE_TRT_RES 0.005
#define DEATH_RATE_TRT_RES 0.002

// copy_number_alteration_model or cna_fitness_model
#define CNA_MODEL empty_cna_function;
#define CNA_PROB 0.00

Cell_Definition fixed_cell;
Cell_Definition stationary_cell;
Cell_Definition motile_cell;
Cycle_Model birth_death;

void create_cell_types( void )
{

	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs


	// housekeeping
	create_birthdeath_model();
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type

	cell_defaults.type = 0; // assign type and increment
	cell_defaults.name = "tumor cell";

	// End adding new live phase info for cell death/removal
	cell_defaults.functions.cycle_model = birth_death; // NOTE: Uses stochastic splitting according to birth death process
	cell_defaults.genotype.genotype_id="0";
	cell_defaults.genotype.alter_genotype = CNA_MODEL;
	cell_defaults.genotype.cna_prob = CNA_PROB;

	// set default_cell_functions;
	// cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.volume_update_function = empty_function;

	//* only needed for a 2-D simulation:
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;
	//*/

	// make sure the defaults are self-consistent.
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	// Set cell size to have radius of 10 - note radius is updated based on volume
	double radius = 10.0;
	cell_defaults.phenotype.geometry.radius = radius;
	static double four_thirds_pi = 4.188790204786391;
	cell_defaults.phenotype.volume.total = pow(radius, 3) * four_thirds_pi;

	// first find index for a few key variables.
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );

	// initially no necrosis
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	// set oxygen uptake / secretion parameters for the default cell type
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10;
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;


//	cell_definition_vector.push_back(cell_defaults);
	// Now, let's define another cell type.
	// It's best to just copy the default and modify it.

	// make this cell type randomly motile, less adhesive, greater survival,
	// and less proliferative
	motile_cell = cell_defaults;
	motile_cell.name = "motile tumor cell";
	motile_cell.type = 1;

	// make sure the new cell type has its own reference phenotyhpe

	motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype );

	// enable random motility
	motile_cell.phenotype.motility.is_motile = true;//true;
	motile_cell.phenotype.motility.persistence_time = 0; // 15 minutes
	motile_cell.phenotype.motility.migration_speed = 0.0; // 0.25 micron/minute
	motile_cell.phenotype.motility.migration_bias = 0.0;// completely random

	// Set cell-cell adhesion to 5% of other cells

	// NOTE: To alter the transition rate of a specific types
	//motile_cell.phenotype.cycle.data.transition_rate(0,0) = 1.0;

	// Set apoptosis to zero
	motile_cell.phenotype.death.rates[apoptosis_model_index] = 0.0;

	motile_cell.phenotype.cycle.data.transition_rate(0, 0) = 0.001;
	motile_cell.phenotype.cycle.data.transition_rate(0, 1) = 0;
	motile_cell.functions.custom_cell_rule = add_velocity;


	stationary_cell = motile_cell;
	stationary_cell.type = 2;
	stationary_cell.name = "stationary tumor cell";
	stationary_cell.phenotype.motility.is_motile = false;
	stationary_cell.phenotype.motility.persistence_time = 0; // 15 minutes
	stationary_cell.phenotype.motility.migration_speed = 0.0; // 0.25 micron/minute
	stationary_cell.phenotype.motility.migration_bias = 0.0;// completely random
	stationary_cell.functions.custom_cell_rule = NULL;
	stationary_cell.phenotype.cycle.data.transition_rate(0, 0) = 0;
	stationary_cell.phenotype.cycle.data.transition_rate(0, 1) = 0;
	stationary_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 10;

	fixed_cell = motile_cell;
	fixed_cell.type = 3;
	fixed_cell.name = "fixed tumor cell";
	fixed_cell.phenotype.motility.is_motile = false;
	fixed_cell.phenotype.motility.persistence_time = 0; // 15 minutes
	fixed_cell.phenotype.motility.migration_speed = 0.0; // 0.25 micron/minute
	fixed_cell.phenotype.motility.migration_bias = 0.0;// completely random
	fixed_cell.functions.custom_cell_rule = add_velocity;
	fixed_cell.phenotype.cycle.data.transition_rate(0, 0) = 0;
	fixed_cell.phenotype.cycle.data.transition_rate(0, 1) = 0;
	fixed_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 10;



	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

	default_microenvironment_options.X_range = {-200, 200};
	default_microenvironment_options.Y_range = {-200, 200};
	default_microenvironment_options.Z_range = {-10, 10};
	default_microenvironment_options.simulate_2D = true; // 3D!

	// no gradients need for this example

	default_microenvironment_options.calculate_gradients = false;

	// set Dirichlet conditions

	default_microenvironment_options.outer_Dirichlet_conditions = true;

	// if there are more substrates, resize accordingly
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{
	// create some cells near the origin

	Cell* pC;

	pC = create_cell( motile_cell );
	pC->assign_position( 0.0, 0.0, 0.0 );
	pC->genotype.genotype_id="0";

  /* Commented out to initialize with single cell
	pC = create_cell();
	pC->assign_position( -2.0, -6.0, 1.0 );

	pC = create_cell();
	pC->assign_position( -5.0, 8.0, -7.0 );

	// now create a motile cell

	pC = create_cell( motile_cell );
	pC->assign_position( 5.0, -8.0, 3.0 );
  */

	return;
}

// function to import from a file
void setup_tissue( std::string input_file )
{
	// create some cells near the origin
	std::ifstream infile(input_file);
	Cell* pC;
	int iter = 1;

	while (infile)
  {
    std::string s;
    if (!getline( infile, s )) break;

    std::istringstream ss( s );
		int x_coord;
		int y_coord;
		int z_coord;
		std::string cell_type;

		ss >> x_coord;
		ss >> y_coord;
		ss >> z_coord;
		ss >> cell_type;
		if(cell_type == "stationary")
		{
			pC = create_cell( stationary_cell ); // sensitive cell definition
		}
		else if(cell_type == "motile")
		{
			pC = create_cell( motile_cell ); // resistant cell definition
		}
		else if(cell_type == "fixed")
		{
			pC = create_cell(fixed_cell); // motile cell definition
		}
		else
		{
			pC = create_cell();
		}

		pC->assign_position( x_coord, y_coord, z_coord );
		pC->genotype.genotype_id="0"+ std::to_string(iter);
  }
	return;
}

void create_birthdeath_model( void )
{
	// Create a 2-phase model for live and "dead" cells
	// When live cells die, they transition to dead cells but are
	// removed in the process - trick to get around issues with
	// default exit upon division options for live cell types
	birth_death.code = PhysiCell_constants::live_cells_cycle_model;
	birth_death.name = "BirthDeath";

	birth_death.data.time_units = "min";

	birth_death.add_phase( PhysiCell_constants::live , "Live" );
	birth_death.add_phase(PhysiCell_constants::custom_phase , "Dead");

	birth_death.add_phase_link( 0 , 0 , NULL );
	birth_death.transition_rate(0, 0) = BIRTH_RATE_SENS;
	birth_death.add_phase_link( 0 , 1 , NULL );
	birth_death.transition_rate(0, 1) = DEATH_RATE_SENS;

	birth_death.phases[0].entry_function = standard_live_phase_entry_function;

	birth_death.phases[0].division_at_phase_exit = false;
	birth_death.phases[1].division_at_phase_exit = false;


	birth_death.phase_links[0][0].exit_function = phase_link_division;
	birth_death.phase_links[0][1].exit_function = phase_link_death;

	return;
}


void phase_link_death( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.flagged_for_removal = true;
	return;
}

void phase_link_division( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.flagged_for_division = true;
	return;
}


void empty_cna_function( Cell* pCell, Genotype& genotype)
{
	return;
}

void add_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
	if(pCell->type == 2)
	{
		return;
	}
	if(pCell->type == 3)
	{
		pCell->velocity = {0.0, 0.0, 0.0};
	}
	else
	{
		//std::vector<double> new_velocity = {-0.3, 0.0, 0.0};
		//pCell->velocity += new_velocity;
	}
}
