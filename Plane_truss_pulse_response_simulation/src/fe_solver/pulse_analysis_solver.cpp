#include "pulse_analysis_solver.h"

pulse_analysis_solver::pulse_analysis_solver()
{
	// Empty constructor
}

pulse_analysis_solver::~pulse_analysis_solver()
{
	// Empty destructor
}

void pulse_analysis_solver::pulse_analysis_start(const nodes_list_store& model_nodes,
	const elementline_list_store& model_lineelements,
	const nodeconstraint_list_store& model_constarints,
	const nodeload_list_store& model_loads,
	const nodepointmass_list_store& model_ptmass,
	const nodeinlcond_list_store& model_inlcond,
	const std::unordered_map<int, material_data>& material_list,
	const bool& is_include_consistent_mass_matrix,
	const modal_analysis_solver& md_solver,
	const modal_analysis_result_store& modal_results,
	const double total_simulation_time,
	const double time_interval,
	const double damping_ratio,
	pulse_analysis_result_store& pulse_response_result,
	pulse_nodes_list_store& pulse_result_nodes,
	pulse_elementline_list_store& pulse_result_lineelements,
	bool& is_pulse_analysis_complete)
{
	// Main solver call
	is_pulse_analysis_complete = false;
	pulse_response_result.clear_results();

	// Check the model
	// Number of loads, initial condition (Exit if no load and no initial condition is present)
	if (model_loads.load_count == 0 && model_inlcond.inlcond_count == 0)
	{
		return;
	}

	// Assign the node id map
	this->nodeid_map = md_solver.nodeid_map;

	//--------------------------------------------------------------------------------------------------------------------
	// Create modal reduced intial condition matrices
	Eigen::MatrixXd modal_reducedInitialDisplacementMatrix(md_solver.reducedDOF, 1);
	Eigen::MatrixXd modal_reducedInitialVelocityMatrix(md_solver.reducedDOF, 1);

	create_initial_condition_matrices(modal_reducedInitialDisplacementMatrix,
		modal_reducedInitialVelocityMatrix,
		model_inlcond,
		model_nodes,
		md_solver.globalDOFMatrix,
		md_solver.reduced_eigenVectorsMatrix,
		md_solver.numDOF,
		md_solver.reducedDOF);
	//___________________________________________________________________________________
	// Create a file to keep track of frequency response matrices
	std::ofstream output_file;
	output_file.open("pulse_analysis_results.txt");

	if (print_matrix == true)
	{
		// Print the Modal Initial Displacement matrix
		output_file << "Modal Initial Displacement Matrix" << std::endl;
		output_file << modal_reducedInitialDisplacementMatrix << std::endl;
		output_file << std::endl;

		// Print the Modal Initial Velocity matrix
		output_file << "Modal Initial Velocity Matrix" << std::endl;
		output_file << modal_reducedInitialVelocityMatrix << std::endl;
		output_file << std::endl;
	}

	//____________________________________________________________________________________________________________________
	// Create the Pulse force data for all the individual 
	std::vector<pulse_load_data> pulse_loads(model_loads.load_count);
	int k = 0;

	for (auto& ld_m : model_loads.loadMap)
	{
		load_data ld = ld_m.second; // get the load data

		pulse_loads[k].load_id = k;
		create_pulse_load_matrices(pulse_loads[k],
			ld,
			model_nodes,
			md_solver.globalDOFMatrix,
			md_solver.reduced_eigenVectorsMatrix,
			md_solver.numDOF,
			md_solver.reducedDOF);

		k++; // iterate load id
	}

	//____________________________________________________________________________________________________________________
	// Pulse Response
	std::unordered_map<int, pulse_node_result> node_results;
	int r_id = 0;

	for (double time_t = 0.0; time_t <= total_simulation_time; time_t = time_t + time_interval)
	{

		Eigen::MatrixXd displ_ampl_RespMatrix_reduced_b4eig_trans(md_solver.reducedDOF, 1);
		displ_ampl_RespMatrix_reduced_b4eig_trans.setZero();

		for (int i = 0; i < md_solver.reducedDOF; i++)
		{
			double displ_resp_initial = 0.0; // Displacement response due to initial condition

			get_steady_state_initial_condition_soln(displ_resp_initial,
				time_t,
				md_solver.modalMass(i),
				md_solver.modalStiff(i),
				modal_reducedInitialDisplacementMatrix(i, 0),
				modal_reducedInitialVelocityMatrix(i, 0));

			//_______________________________________________________________________
			double displ_resp_force = 0.0; // Displacement response due to pulse force

			// get all the loads
			for (auto& pulse_load : pulse_loads)
			{
				// Go through all the force
				double at_force_displ_resp = 0.0;

				get_steady_state_pulse_soln(at_force_displ_resp,
					time_t,
					md_solver.modalMass(i),
					md_solver.modalStiff(i),
					pulse_load.modal_reducedLoadamplMatrix(i, 0),
					pulse_load.load_start_time,
					pulse_load.load_end_time);

				displ_resp_force = displ_resp_force + at_force_displ_resp;
			}

			// Add to the modal displ matrix
			displ_ampl_RespMatrix_reduced_b4eig_trans.coeffRef(i, 0) = displ_resp_initial + displ_resp_force ;
		}

		// Apply modal de-transformation
		Eigen::MatrixXd displ_ampl_RespMatrix_reduced(md_solver.reducedDOF, 1);
		displ_ampl_RespMatrix_reduced.setZero();

		displ_ampl_RespMatrix_reduced = md_solver.reduced_eigenVectorsMatrix * displ_ampl_RespMatrix_reduced_b4eig_trans;

		// Extend reduced modal displ matrix to global modal displ matrix
		Eigen::MatrixXd displ_ampl_RespMatrix_b4supp_trans(md_solver.numDOF, 1);
		displ_ampl_RespMatrix_b4supp_trans.setZero();

		get_global_resp_matrix(displ_ampl_RespMatrix_b4supp_trans,
			displ_ampl_RespMatrix_reduced,
			md_solver.globalDOFMatrix,
			md_solver.numDOF,
			md_solver.reducedDOF);

		// Apply support transformation
		Eigen::MatrixXd displ_ampl_RespMatrix(md_solver.numDOF, 1);
		displ_ampl_RespMatrix.setZero();

		displ_ampl_RespMatrix = md_solver.globalSupportInclinationMatrix * displ_ampl_RespMatrix_b4supp_trans;

		// Store the results to node results
		for (auto& nd_m : model_nodes.nodeMap)
		{
			// get the node id
			int nd_id = nd_m.second.node_id;
			int nd_index = nodeid_map[nd_id];

			// Node displacement response
			glm::vec2 node_displ = glm::vec2(displ_ampl_RespMatrix((nd_index * 2) + 0, 0),
				displ_ampl_RespMatrix((nd_index * 2) + 1, 0));

			// Add the index
			node_results[nd_id].index.push_back(r_id);
			// Add the time val
			node_results[nd_id].time_val.push_back(time_t);
			// Add the displacement
			node_results[nd_id].node_pulse_displ.push_back(node_displ);
		}

		r_id++;
	}

	// Map the results
	map_pulse_analysis_results(pulse_response_result,
		pulse_result_nodes,
		pulse_result_lineelements,
		r_id,
		model_nodes,
		model_lineelements,
		node_results);

	pulse_response_result.set_analysis_setting(r_id, time_interval, total_simulation_time);

	if (pulse_result_lineelements.max_line_displ == 0)
	{
		// Analysis failed 
		return;
	}

	// Analysis complete
	is_pulse_analysis_complete = true;
}

void pulse_analysis_solver::create_initial_condition_matrices(Eigen::MatrixXd& modal_reducedInitialDisplacementMatrix,
	Eigen::MatrixXd& modal_reducedInitialVelocityMatrix,
	const nodeinlcond_list_store& model_inlcond,
	const nodes_list_store& model_nodes,
	const Eigen::MatrixXd& globalDOFMatrix,
	const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
	const int& numDOF,
	const int& reducedDOF)
{
	// Create a global initial condition matrix
	Eigen::MatrixXd globalInitialDisplacementMatrix(numDOF, 1);
	Eigen::MatrixXd globalInitialVelocityMatrix(numDOF, 1);

	globalInitialDisplacementMatrix.setZero();
	globalInitialVelocityMatrix.setZero();

	for (auto& inlc_m : model_inlcond.inlcondMap)
	{
		nodeinl_condition_data inlc = inlc_m.second;

		// get the matrix id
		int n_id = nodeid_map[inlc.node_id]; // get the ordered map of the start node ID

		// Create a node initial displacement matrix
		Eigen::MatrixXd nodeinitialDisplacementMatrix(2, 1);
		nodeinitialDisplacementMatrix.setZero();

		nodeinitialDisplacementMatrix.coeffRef(0, 0) = inlc.inl_displacement_x;
		nodeinitialDisplacementMatrix.coeffRef(1, 0) = inlc.inl_displacement_y;

		// global initial displacement matrix
		globalInitialDisplacementMatrix.block<2, 1>(n_id * 2, 0) += nodeinitialDisplacementMatrix.block<2, 1>(0, 0);

		// Create a node initial velocity matrix
		Eigen::MatrixXd nodeinitialVelocityMatrix(2, 1);
		nodeinitialVelocityMatrix.setZero();

		nodeinitialVelocityMatrix.coeffRef(0, 0) = inlc.inl_velocity_x;
		nodeinitialVelocityMatrix.coeffRef(1, 0) = inlc.inl_velocity_y;

		// global initial velocity matrix
		globalInitialVelocityMatrix.block<2, 1>(n_id * 2, 0) += nodeinitialVelocityMatrix.block<2, 1>(0, 0);
	}

	// Reduce the intial condition matrix with the degree of freedom 
	Eigen::MatrixXd reducedInitialDisplacementMatrix(reducedDOF, 1);
	Eigen::MatrixXd reducedInitialVelocityMatrix(reducedDOF, 1);

	reducedInitialDisplacementMatrix.setZero();
	reducedInitialVelocityMatrix.setZero();

	// reduced initial displacement matrix
	get_reduced_global_matrix(reducedInitialDisplacementMatrix,
		globalInitialDisplacementMatrix,
		globalDOFMatrix,
		numDOF,
		reducedDOF);

	// reduced initial velocity matrix
	get_reduced_global_matrix(reducedInitialVelocityMatrix,
		globalInitialVelocityMatrix,
		globalDOFMatrix,
		numDOF,
		reducedDOF);

	// apply modal decomposition of the initial displacements
	modal_reducedInitialDisplacementMatrix = reduced_eigenVectorsMatrix * reducedInitialDisplacementMatrix;
	modal_reducedInitialVelocityMatrix = reduced_eigenVectorsMatrix * reducedInitialVelocityMatrix;
}

void pulse_analysis_solver::create_pulse_load_matrices(pulse_load_data& pulse_loads,
	const load_data& ld,
	const nodes_list_store& model_nodes,
	const Eigen::MatrixXd& globalDOFMatrix,
	const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
	const int& numDOF,
	const int& reducedDOF)
{
	// Create the global load amplitude matrix
	// Extract the line in which the load is applied
	node_store nd = model_nodes.nodeMap.at(ld.node_id);

	// Get the Matrix row ID
	int n_id = nodeid_map[nd.node_id]; // get the ordered map of the start node ID

	// Create a matrix for element load matrix
	Eigen::MatrixXd nodeLoadMatrix(2, 1);
	nodeLoadMatrix.setZero();

	double load_val = ld.load_value; // Load value
	double load_angle_rad = ld.load_angle * (m_pi / 180.0f);

	double f_x = load_val * std::cos(load_angle_rad);
	double f_y = load_val * std::sin(load_angle_rad);

	nodeLoadMatrix.coeffRef(0, 0) = f_x;
	nodeLoadMatrix.coeffRef(1, 0) = f_y;

	// global load matrix
	Eigen::MatrixXd globalLoadamplMatrix(numDOF, 1);
	globalLoadamplMatrix.setZero();

	globalLoadamplMatrix.block<2, 1>(n_id * 2, 0) += nodeLoadMatrix.block<2, 1>(0, 0);

	//______________________________________________________________________________________________
	// Reduce the global load matrix with DOF
	Eigen::MatrixXd reducedLoadamplMatrix(reducedDOF, 1);
	reducedLoadamplMatrix.setZero();

	int r = 0;
	// Loop throug the Degree of freedom of indices
	for (int i = 0; i < numDOF; i++)
	{
		if (globalDOFMatrix(i, 0) == 0)
		{
			// constrained row index, so skip
			continue;
		}
		else
		{
			reducedLoadamplMatrix.coeffRef(r, 0) = globalLoadamplMatrix(i, 0);
			r++;
		}
	}

	//______________________________________________________________________________________________
	// Apply modal transformation to the reduced load ampl matrix

	Eigen::MatrixXd modal_reducedLoadamplMatrix(reducedDOF, 1);
	modal_reducedLoadamplMatrix.setZero();


	modal_reducedLoadamplMatrix = reduced_eigenVectorsMatrix * reducedLoadamplMatrix;

	//______________________________________________________________________________________________
	// Copy the data to pulse load data variable
	pulse_loads.load_id = ld.load_id;
	pulse_loads.load_start_time = ld.load_start_time;
	pulse_loads.load_end_time = ld.load_end_time;
	pulse_loads.globalLoadamplMatrix = globalLoadamplMatrix; // Global load matrix for this load
	pulse_loads.reducedLoadamplMatrix = reducedLoadamplMatrix; // Reduced load matrix with constraint matrix
	pulse_loads.modal_reducedLoadamplMatrix = modal_reducedLoadamplMatrix; // Modal reduction applied to load matrix
}

void pulse_analysis_solver::get_reduced_global_matrix(Eigen::MatrixXd& reducedglobalMatrix,
	const Eigen::MatrixXd& globalMatrix,
	const Eigen::MatrixXd& globalDOFMatrix,
	const int& numDOF,
	const int& reducedDOF)
{
	// Get the reduced global matrix with the Degree of freedom
	int r = 0;

	// Loop throug the Degree of freedom of indices
	for (int i = 0; i < numDOF; i++)
	{
		if (globalDOFMatrix(i, 0) == 0)
		{
			// constrained row index, so skip
			continue;
		}
		else
		{
			// Get the reduced matrices
			reducedglobalMatrix.coeffRef(r, 0) = globalMatrix.coeffRef(i, 0);
			r++;
		}
	}
}

void pulse_analysis_solver::get_global_resp_matrix(Eigen::MatrixXd& displ_ampl_RespMatrix_b4supp_trans,
	const Eigen::MatrixXd& displ_ampl_RespMatrix_reduced,
	const Eigen::MatrixXd& globalDOFMatrix,
	const int& numDOF,
	const int& reducedDOF)
{
	// Get global response matrix from the reduced matrices
	// Loop throug the Degree of freedom of indices
	int r = 0;

	// Loop throug the Degree of freedom of indices
	for (int i = 0; i < numDOF; i++)
	{
		if (globalDOFMatrix(i, 0) == 0)
		{
			// constrained row index, so skip
			continue;
		}
		else
		{
			// Get the reduced matrices
			displ_ampl_RespMatrix_b4supp_trans.coeffRef(i, 0) = displ_ampl_RespMatrix_reduced(r, 0);
			r++;
		}
	}
}

void pulse_analysis_solver::get_steady_state_initial_condition_soln(double& steady_state_displ_resp,
	const double& time_t,
	const double& modal_mass,
	const double& modal_stiff,
	const double& modal_initial_displacement,
	const double& modal_initial_velocity)
{
	// Return the steady state solution for the intial displacment and velocity
	double modal_omega_n = std::sqrt(modal_stiff / modal_mass); // Modal omega n

	steady_state_displ_resp = (modal_initial_displacement * std::cos(modal_omega_n * time_t)) +
		((modal_initial_velocity / modal_omega_n) * std::sin(modal_omega_n * time_t));
}


void pulse_analysis_solver::get_steady_state_pulse_soln(double& steady_state_displ_resp,
	const double& time_t,
	const double& modal_mass,
	const double& modal_stiff,
	const double& modal_force_ampl,
	const double& modal_force_starttime,
	const double& modal_force_endtime)
{
	// Return the steady state solution for the half sine pulse
	double modal_omega_n = std::sqrt(modal_stiff / modal_mass); // Modal omega n
	double modal_omega_f = m_pi / (modal_force_endtime - modal_force_starttime);

	if (std::abs(modal_omega_f - modal_omega_n) < epsilon)
	{
		if (modal_omega_n < epsilon)
		{
			// For infinity solution state
			steady_state_displ_resp = ((modal_force_ampl / modal_omega_f) * (time_t - modal_force_starttime)) -
				((modal_force_ampl / std::pow(modal_omega_f, 2)) * std::sin(modal_omega_f * (time_t - modal_force_starttime)));
		}
		else
		{
			// omega_f / omega_n
			double omega_f_by_omega_n;
			omega_f_by_omega_n = modal_omega_f / modal_omega_n;

			steady_state_displ_resp = ((modal_force_ampl / modal_stiff) / (1 - std::pow(omega_f_by_omega_n, 2))) *
				(std::sin(modal_omega_f * (time_t - modal_force_starttime)) -
					omega_f_by_omega_n * std::sin(modal_omega_n * (time_t - modal_force_starttime)));
		}
	}
	else
	{
		// Resonsance state
		steady_state_displ_resp = 0.5 * (modal_force_ampl / modal_stiff) *
			(std::sin(modal_omega_n * (time_t - modal_force_starttime)) -
				(modal_omega_n * (time_t - modal_force_starttime)) * std::cos(modal_omega_n * (time_t - modal_force_starttime)));
	}
}

void pulse_analysis_solver::map_pulse_analysis_results(pulse_analysis_result_store& pulse_response_result,
	pulse_nodes_list_store& pulse_result_nodes,
	pulse_elementline_list_store& pulse_result_lineelements,
	const int& number_of_time_steps,
	const nodes_list_store& model_nodes,
	const elementline_list_store& model_lineelements,
	const std::unordered_map<int, pulse_node_result>& node_results)
{
	// Map the pulse analysis results
	// map the node results
	pulse_result_nodes.clear_data();

	for (auto& nd_m : model_nodes.nodeMap)
	{
		// Extract the model node
		node_store nd = nd_m.second;

		// Add to the pulse node results store
		pulse_result_nodes.add_result_node(nd.node_id, nd.node_pt, node_results.at(nd.node_id), number_of_time_steps);
	}

	// map the line results
	pulse_result_lineelements.clear_data();

	for (auto& ln_m : model_lineelements.elementlineMap)
	{
		// Extract the model lines
		elementline_store ln = ln_m.second;

		// Extract the pulse node store -> start node and end node
		pulse_node_store* startNode = &pulse_result_nodes.pulse_nodeMap[ln.startNode->node_id];
		pulse_node_store* endNode = &pulse_result_nodes.pulse_nodeMap[ln.endNode->node_id];

		// Add to the pulse element results store
		pulse_result_lineelements.add_pulse_elementline(ln.line_id, startNode, endNode);
	}

	//_________________________________________________________________________________________________________________
	double maximum_displacement = 0.0;

	for (auto& ln_m : pulse_result_lineelements.pulse_elementlineMap)
	{
		// get all the discretized line of every single line
		for (auto& h_ln : ln_m.second.discretized_bar_line_data)
		{
			//get all the two points
			// Point 1 displacement
			for (auto& pt1 : h_ln.pt1_modal_displ)
			{
				double displ1 = std::sqrt(std::pow(pt1.x, 2) + std::pow(pt1.y, 2));

				if (displ1 > maximum_displacement)
				{
					maximum_displacement = displ1;
				}
			}

			// Point 2 displacement
			for (auto& pt2 : h_ln.pt2_modal_displ)
			{
				double displ2 = std::sqrt(std::pow(pt2.x, 2) + std::pow(pt2.y, 2));

				if (displ2 > maximum_displacement)
				{
					maximum_displacement = displ2;
				}
			}
		}
	}

	// Set the maximim displacement
	pulse_result_nodes.max_node_displ = maximum_displacement;
	pulse_result_lineelements.max_line_displ = maximum_displacement;
}