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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"
// put custom code modules here!

#include "./custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	bool XML_status = false;
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// PNRG setup
	SeedRandom( );

	// time setup
	std::string time_units = "min";

	/* XXX REMOVED
	double t = 0.0; // current simulation time
	double t_output_interval = 1.0; // Every hour
	double t_max = 150;//5 * 60 * 24; // 10 days
	double t_next_output_time = t;
	int output_index = 0; // used for creating unique output filenames
	*/

	/* Microenvironment setup */

	setup_microenvironment(); // modify this in the custom code

	/* PhysiCell setup */

	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 20;
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	/* Users typically start modifying here. START USERMODS */

	create_cell_types();

	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */

	// set MultiCellDS save options

	set_save_biofvm_mesh_as_matlab( false );
	set_save_biofvm_data_as_matlab( false );
	set_save_biofvm_cell_data( false );
	set_save_biofvm_cell_data_as_custom_matlab( false );

	// set the performance timers
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	std::ofstream report_file;
	char filename[1024];
	sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() );
	report_file.open(filename); 	// create the data log file
	report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;

	// main loop
	try
	{
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time.
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout );
				log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);

				char filename[1024];
				sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

				PhysiCell_globals.full_output_index++;
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}

			// run PhysiCell
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			PhysiCell_globals.current_time += diffusion_dt;
		}
		log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
		report_file.close();
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}

	// save a final simulation snapshot

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );
	std::cout << std::endl;

	return 0;
}
