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
	const std::unordered_map<int, material_data>& material_list,
	const bool& is_include_consistent_mass_matrix,
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
	// Number of loads (Exit if no load is present)
	if (model_loads.load_count == 0)
	{
		return;
	}

	// Create global stiffness, mass and load matrices
	//___________________________________________________________________________________
	// Create a node ID map (to create a nodes as ordered and numbered from 0,1,2...n)
	int i = 0;
	for (auto& nd : model_nodes.nodeMap)
	{
		nodeid_map[nd.first] = i;
		i++;
	}

	//___________________________________________________________________________________
	// Create a file to keep track of frequency response matrices
	std::ofstream output_file;
	output_file.open("pulse_analysis_results.txt");

	//____________________________________________________________________________________________________________________
	int numDOF = model_nodes.node_count * 3; // Number of degrees of freedom (3 DOFs per node (2 translation and 1 rotation)

	// Global Stiffness Matrix
	Eigen::MatrixXd globalStiffnessMatrix(numDOF, numDOF);
	globalStiffnessMatrix.setZero();

	get_global_stiffness_matrix(globalStiffnessMatrix,
		model_lineelements,
		material_list,
		model_constarints,
		output_file);


	//____________________________________________________________________________________________________________________
	// Global Point Mass Matrix
	Eigen::MatrixXd globalPointMassMatrix(numDOF, numDOF);
	globalPointMassMatrix.setZero();

	get_global_pointmass_matrix(globalPointMassMatrix,
		model_nodes,
		model_ptmass,
		output_file);

	//____________________________________________________________________________________________________________________
	// Global Consistent Mass Matrix
	Eigen::MatrixXd globalConsistentMassMatrix(numDOF, numDOF);
	globalConsistentMassMatrix.setZero();

	if (is_include_consistent_mass_matrix == true)
	{
		get_global_consistentmass_matrix(globalConsistentMassMatrix,
			model_lineelements,
			material_list,
			model_constarints,
			output_file);
	}

	//____________________________________________________________________________________________________________________
	// Global Consistent Mass Matrix
	Eigen::MatrixXd globalMassMatrix(numDOF, numDOF);
	globalMassMatrix.setZero();

	globalMassMatrix = globalPointMassMatrix + globalConsistentMassMatrix;

	//____________________________________________________________________________________________________________________
	// Global DOF Mass Matrix
	Eigen::MatrixXd globalDOFMatrix(numDOF, 1);
	globalDOFMatrix.setZero();

	// Determine the size of the reduced stiffness matrix based on the number of unconstrained degrees of freedom
	int reducedDOF = 0;

	get_global_dof_matrix(globalDOFMatrix,
		model_nodes,
		model_constarints,
		reducedDOF,
		output_file);

	//____________________________________________________________________________________________________________________
	// Create Reduced Global Mass and stiffness matrix
	Eigen::MatrixXd  reduced_globalStiffnessMatrix(reducedDOF, reducedDOF);
	reduced_globalStiffnessMatrix.setZero();

	// Reduced Global Mass matrix
	Eigen::MatrixXd reduced_globalMassMatrix(reducedDOF, reducedDOF);
	reduced_globalMassMatrix.setZero();

	get_reduced_global_matrices(reduced_globalStiffnessMatrix,
		reduced_globalMassMatrix,
		globalStiffnessMatrix,
		globalMassMatrix,
		globalDOFMatrix,
		numDOF,
		output_file);

	//____________________________________________________________________________________________________________________
	// Modal Decomposition

	// Reduced Eigen Vectors matrix
	Eigen::MatrixXd reduced_eigenVectorsMatrix(reducedDOF, reducedDOF);
	reduced_eigenVectorsMatrix.setZero();

	get_reduced_modal_vector_matrix(reduced_eigenVectorsMatrix, modal_results, reducedDOF, output_file);

	//____________________________________________________________________________________________________________________
	// Create modal matrices
	Eigen::VectorXd modalMass(reducedDOF);
	Eigen::VectorXd modalStiff(reducedDOF);

	get_modal_matrices(modalMass,
		modalStiff,
		reduced_eigenVectorsMatrix,
		reduced_globalMassMatrix,
		reduced_globalStiffnessMatrix,
		reducedDOF,
		output_file);

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
			model_lineelements,
			globalDOFMatrix,
			reduced_eigenVectorsMatrix,
			numDOF,
			reducedDOF);

		k++; // iterate load id
	}

	//____________________________________________________________________________________________________________________
	// Create the global support inclination matrix
	Eigen::MatrixXd globalSupportInclinationMatrix(numDOF, numDOF);
	globalSupportInclinationMatrix.setZero();


	get_globalSupportInclinationMatrix(globalSupportInclinationMatrix,
		model_nodes,
		model_constarints,
		numDOF,
		output_file);

	//____________________________________________________________________________________________________________________
	// Pulse Response
	std::unordered_map<int, pulse_node_result> node_results;
	int r_id = 0;

	for (double time_t = 0.0; time_t <= total_simulation_time; time_t = time_t + time_interval)
	{

		Eigen::MatrixXd displ_ampl_RespMatrix_reduced_b4eig_trans(reducedDOF, 1);
		displ_ampl_RespMatrix_reduced_b4eig_trans.setZero();

		for (int i = 0; i < reducedDOF; i++)
		{
			double total_displ_resp = 0.0;

			// get all the loads
			for (auto& pulse_load : pulse_loads)
			{
				// Go through all the force
				double at_force_displ_resp = 0.0;

				get_steady_state_pulse_soln(at_force_displ_resp,
					time_t,
					modalMass(i),
					modalStiff(i),
					pulse_load.modal_reducedLoadamplMatrix(i, 0),
					pulse_load.load_start_time,
					pulse_load.load_end_time);

				total_displ_resp = total_displ_resp + at_force_displ_resp;
			}

			// Add to the modal displ matrix
			displ_ampl_RespMatrix_reduced_b4eig_trans.coeffRef(i, 0) = total_displ_resp;
		}

		// Apply modal de-transformation
		Eigen::MatrixXd displ_ampl_RespMatrix_reduced(reducedDOF, 1);
		displ_ampl_RespMatrix_reduced.setZero();

		displ_ampl_RespMatrix_reduced = reduced_eigenVectorsMatrix * displ_ampl_RespMatrix_reduced_b4eig_trans;

		// Extend reduced modal displ matrix to global modal displ matrix
		Eigen::MatrixXd displ_ampl_RespMatrix_b4supp_trans(numDOF, 1);
		displ_ampl_RespMatrix_b4supp_trans.setZero();

		get_global_resp_matrix(displ_ampl_RespMatrix_b4supp_trans,
			displ_ampl_RespMatrix_reduced, globalDOFMatrix,
			numDOF,
			reducedDOF);

		// Apply support transformation
		Eigen::MatrixXd displ_ampl_RespMatrix(numDOF, 1);
		displ_ampl_RespMatrix.setZero();

		displ_ampl_RespMatrix = globalSupportInclinationMatrix * displ_ampl_RespMatrix_b4supp_trans;

		// Store the results to node results
		for (auto& nd_m : model_nodes.nodeMap)
		{
			// get the node id
			int nd_id = nd_m.second.node_id;
			int nd_index = nodeid_map[nd_id];

			// Node displacement response
			glm::vec3 node_displ = glm::vec3(displ_ampl_RespMatrix((nd_index * 3) + 0, 0),
				displ_ampl_RespMatrix((nd_index * 3) + 1, 0),
				displ_ampl_RespMatrix((nd_index * 3) + 2, 0));

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


void pulse_analysis_solver::get_global_stiffness_matrix(Eigen::MatrixXd& globalStiffnessMatrix,
	const elementline_list_store& model_lineelements,
	const std::unordered_map<int, material_data>& material_list,
	const nodeconstraint_list_store& model_constarints, std::ofstream& output_file)
{
	// Create global stiffness matrix
	for (auto& ln_m : model_lineelements.elementlineMap)
	{
		// Create the element stiffness matrix
		elementline_store ln = ln_m.second;
		material_data elementline_material = material_list.at(ln.material_id);

		// Create a matrix for element stiffness matrix
		Eigen::MatrixXd elementStiffnessMatrix(6, 6);
		elementStiffnessMatrix.setZero();

		get_element_stiffness_matrix(elementStiffnessMatrix, ln, elementline_material, model_constarints, output_file);

		// Get the Node ID
		int sn_id = nodeid_map[ln.startNode->node_id]; // get the ordered map of the start node ID
		int en_id = nodeid_map[ln.endNode->node_id]; // get the ordered map of the end node ID

		globalStiffnessMatrix.block<3, 3>(sn_id * 3, sn_id * 3) += elementStiffnessMatrix.block<3, 3>(0, 0);
		globalStiffnessMatrix.block<3, 3>(sn_id * 3, en_id * 3) += elementStiffnessMatrix.block<3, 3>(0, 3);
		globalStiffnessMatrix.block<3, 3>(en_id * 3, sn_id * 3) += elementStiffnessMatrix.block<3, 3>(3, 0);
		globalStiffnessMatrix.block<3, 3>(en_id * 3, en_id * 3) += elementStiffnessMatrix.block<3, 3>(3, 3);
	}
}

void pulse_analysis_solver::get_element_stiffness_matrix(Eigen::MatrixXd& elementStiffnessMatrix,
	const elementline_store& ln,
	const material_data& elementline_material,
	const nodeconstraint_list_store& model_constarints,
	std::ofstream& output_file)
{
	// Get element stiffness matrix
	// Compute the differences in x and y coordinates
	double dx = ln.endNode->node_pt.x - ln.startNode->node_pt.x;
	double dy = -1.0 * (ln.endNode->node_pt.y - ln.startNode->node_pt.y);

	// Compute the length of the truss element
	double eLength = std::sqrt((dx * dx) + (dy * dy));

	// Compute the direction cosines
	double Lcos = (dx / eLength);
	double Msin = (dy / eLength);

	//_________________________________________________________
	// Local -> Global transformation matrix
	Eigen::MatrixXd L_transformation_matrix(6, 6);
	L_transformation_matrix.setZero();

	L_transformation_matrix.row(0) = Eigen::RowVectorXd({ {Lcos, Msin, 0.0, 0.0, 0.0, 0.0 } });
	L_transformation_matrix.row(1) = Eigen::RowVectorXd({ { -Msin, Lcos, 0.0, 0.0, 0.0, 0.0} });
	L_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	L_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, Lcos, Msin, 0.0} });
	L_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, -Msin, Lcos, 0.0} });
	L_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });


	//_________________________________________________________
	// Local element stiffness matrix
	Eigen::MatrixXd local_element_stiffness_matrix(6, 6);
	local_element_stiffness_matrix.setZero();

	double k1 = (elementline_material.youngs_mod * elementline_material.cs_area) / eLength;
	double k2 = (1.0) / (eLength * eLength * eLength);
	double k3 = (1.0) / (eLength * eLength);
	double k4 = (1.0) / eLength;

	local_element_stiffness_matrix.row(0) = Eigen::RowVectorXd({ {k1, 0.0, 0.0, -1.0 * k1, 0.0, 0.0} });
	local_element_stiffness_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 12.0 * k2, 6.0 * k3, 0.0, -12.0 * k2, 6.0 * k3} });
	local_element_stiffness_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 6.0 * k3, 4.0 * k4, 0.0, -6.0 * k3, 2.0 * k4} });
	local_element_stiffness_matrix.row(3) = Eigen::RowVectorXd({ {-1.0 * k1, 0.0, 0.0, k1, 0.0, 0.0} });
	local_element_stiffness_matrix.row(4) = Eigen::RowVectorXd({ {0.0, -12.0 * k2, -6.0 * k3, 0.0, 12.0 * k2, -6.0 * k3} });
	local_element_stiffness_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 6.0 * k3, 2.0 * k4, 0.0, -6.0 * k3, 4.0 * k4} });

	//_________________________________________________________
	// Transformed element stiffness matrix
	Eigen::MatrixXd e_stiffness_matrix(6, 6);
	e_stiffness_matrix.setZero();

	e_stiffness_matrix = L_transformation_matrix.transpose() * local_element_stiffness_matrix * L_transformation_matrix;
	//_________________________________________________________

	// Transformation matrices to include support inclinatation
	Eigen::MatrixXd s_transformation_matrix(6, 6);
	s_transformation_matrix.setZero(); // support inclination transformation matrix

	int constraint_type;
	double constraint_angle_rad;
	double support_Lcos;
	double support_Msin;

	// Start node support inclination
	if (model_constarints.constraintMap.find(ln.startNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the start node
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {1.0, 0.0, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 1.0, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.startNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.startNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {support_Lcos, -support_Msin, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {support_Msin, support_Lcos, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	}

	// End node support inclination
	if (model_constarints.constraintMap.find(ln.endNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the end node
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 1.0, 0.0, 0.0} });
		s_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 1.0, 0.0} });
		s_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.endNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.endNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, support_Lcos, -support_Msin, 0.0} });
		s_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, support_Msin, support_Lcos, 0.0} });
		s_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });
	}

	// Calculate the element stiffness matrix
	elementStiffnessMatrix = s_transformation_matrix.transpose() * e_stiffness_matrix * s_transformation_matrix;

}

void pulse_analysis_solver::get_global_pointmass_matrix(Eigen::MatrixXd& globalPointMassMatrix,
	const nodes_list_store& model_nodes,
	const nodepointmass_list_store& model_ptmass,
	std::ofstream& output_file)
{
	// Create a global point mass matrix
	for (auto& nd_m : model_nodes.nodeMap)
	{
		// Get the node data
		node_store nd = nd_m.second;
		int nd_map = nodeid_map[nd.node_id]; // get the ordered map of the node ID

		if (model_ptmass.ptmassMap.find(nd.node_id) != model_ptmass.ptmassMap.end())
		{
			// Nodes have point mass
			nodepointmass_data ptm = model_ptmass.ptmassMap.at(nd.node_id);

			globalPointMassMatrix((nd_map * 3) + 0, (nd_map * 3) + 0) = ptm.ptmass_x;
			globalPointMassMatrix((nd_map * 3) + 1, (nd_map * 3) + 1) = ptm.ptmass_y;
			globalPointMassMatrix((nd_map * 3) + 2, (nd_map * 3) + 2) = 0.0;
		}
		else
		{
			// Nodes doesnt have point mass
			globalPointMassMatrix((nd_map * 3) + 0, (nd_map * 3) + 0) = 0.0;
			globalPointMassMatrix((nd_map * 3) + 1, (nd_map * 3) + 1) = 0.0;
			globalPointMassMatrix((nd_map * 3) + 2, (nd_map * 3) + 2) = 0.0;
		}
	}
}


void pulse_analysis_solver::get_global_consistentmass_matrix(Eigen::MatrixXd& globalConsistentMassMatrix,
	const elementline_list_store& model_lineelements,
	const std::unordered_map<int, material_data>& material_list,
	const nodeconstraint_list_store& model_constarints,
	std::ofstream& output_file)
{
	// Create global consistent mass matrix
	for (auto& ln_m : model_lineelements.elementlineMap)
	{
		// Create the element stiffness matrix
		elementline_store ln = ln_m.second;
		material_data elementline_material = material_list.at(ln.material_id);

		// Create a matrix for element stiffness matrix
		Eigen::MatrixXd elementConsistentMassMatrix(6, 6);
		elementConsistentMassMatrix.setZero();

		get_element_consistentmass_matrix(elementConsistentMassMatrix, ln, elementline_material, model_constarints, output_file);

		// Get the Node ID
		int sn_id = nodeid_map[ln.startNode->node_id]; // get the ordered map of the start node ID
		int en_id = nodeid_map[ln.endNode->node_id]; // get the ordered map of the end node ID

		globalConsistentMassMatrix.block<3, 3>(sn_id * 3, sn_id * 3) += elementConsistentMassMatrix.block<3, 3>(0, 0);
		globalConsistentMassMatrix.block<3, 3>(sn_id * 3, en_id * 3) += elementConsistentMassMatrix.block<3, 3>(0, 3);
		globalConsistentMassMatrix.block<3, 3>(en_id * 3, sn_id * 3) += elementConsistentMassMatrix.block<3, 3>(3, 0);
		globalConsistentMassMatrix.block<3, 3>(en_id * 3, en_id * 3) += elementConsistentMassMatrix.block<3, 3>(3, 3);
	}
}

void pulse_analysis_solver::get_element_consistentmass_matrix(Eigen::MatrixXd& elementConsistentMassMatrix,
	const elementline_store& ln,
	const material_data& elementline_material,
	const nodeconstraint_list_store& model_constarints,
	std::ofstream& output_file)
{
	// Create element consistent mass matrix
	// Compute the differences in x and y coordinates
	double dx = ln.endNode->node_pt.x - ln.startNode->node_pt.x;
	double dy = -1.0 * (ln.endNode->node_pt.y - ln.startNode->node_pt.y);

	// Compute the length of the frame element
	double eLength = std::sqrt((dx * dx) + (dy * dy));

	// Compute the direction cosines
	double Lcos = (dx / eLength);
	double Msin = (dy / eLength);

	//_________________________________________________________
	// Local -> Global transformation matrix
	Eigen::MatrixXd L_transformation_matrix(6, 6);
	L_transformation_matrix.setZero();

	L_transformation_matrix.row(0) = Eigen::RowVectorXd({ {Lcos, Msin, 0.0, 0.0, 0.0, 0.0} });
	L_transformation_matrix.row(1) = Eigen::RowVectorXd({ {-Msin, Lcos, 0.0, 0.0, 0.0, 0.0} });
	L_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	L_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, Lcos, Msin, 0.0} });
	L_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, -Msin, Lcos, 0.0} });
	L_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });

	//_________________________________________________________
	// Local element stiffness matrix
	Eigen::MatrixXd local_element_consistentmass_matrix(6, 6);
	local_element_consistentmass_matrix.setZero();

	double k1 = (elementline_material.mat_density * elementline_material.cs_area * eLength) / 420.0f;
	double k2 = k1 * eLength;
	double k3 = k1 * (eLength * eLength);

	local_element_consistentmass_matrix.row(0) = Eigen::RowVectorXd({ {140.0 * k1, 0.0, 0.0, 70.0 * k1, 0.0, 0.0} });
	local_element_consistentmass_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 156.0 * k1, 22.0 * k2, 0.0, 54.0 * k1, -13.0 * k2} });
	local_element_consistentmass_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 22.0 * k2, 4.0 * k3, 0.0, 13.0 * k2, -3.0 * k3} });
	local_element_consistentmass_matrix.row(3) = Eigen::RowVectorXd({ {70.0 * k1, 0.0, 0.0, 140.0 * k1, 0.0, 0.0} });
	local_element_consistentmass_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 54.0 * k1, 13.0 * k2, 0.0, 156.0 * k1, -22.0 * k2} });
	local_element_consistentmass_matrix.row(5) = Eigen::RowVectorXd({ {0.0, -13.0 * k2, -3.0 * k3, 0.0, -22.0 * k2, 4.0 * k3} });

	//_________________________________________________________
	// Transformed element stiffness matrix
	Eigen::MatrixXd e_consistentMass_matrix(6, 6);
	e_consistentMass_matrix.setZero();

	e_consistentMass_matrix = L_transformation_matrix.transpose() * local_element_consistentmass_matrix * L_transformation_matrix;
	//_________________________________________________________

	// Transformation matrices to include support inclinatation
	Eigen::MatrixXd s_transformation_matrix(6, 6);
	s_transformation_matrix.setZero(); // support inclination transformation matrix

	int constraint_type;
	double constraint_angle_rad;
	double support_Lcos;
	double support_Msin;

	// Start node support inclination
	if (model_constarints.constraintMap.find(ln.startNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the start node
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {1.0, 0.0, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 1.0, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.startNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.startNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {support_Lcos, -support_Msin, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {support_Msin, support_Lcos, 0.0, 0.0, 0.0, 0.0} });
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0} });
	}

	// End node support inclination
	if (model_constarints.constraintMap.find(ln.endNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the end node
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 1.0, 0.0, 0.0} });
		s_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 1.0, 0.0} });
		s_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.endNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.endNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, support_Lcos, -support_Msin, 0.0} });
		s_transformation_matrix.row(4) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, support_Msin, support_Lcos, 0.0 } });
		s_transformation_matrix.row(5) = Eigen::RowVectorXd({ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} });
	}

	// Calculate the element stiffness matrix
	elementConsistentMassMatrix = s_transformation_matrix.transpose() * e_consistentMass_matrix * s_transformation_matrix;
}

void pulse_analysis_solver::get_global_dof_matrix(Eigen::MatrixXd& globalDOFMatrix,
	const nodes_list_store& model_nodes,
	const nodeconstraint_list_store& model_constarints,
	int& reducedDOF,
	std::ofstream& output_file)
{
	// Create global DOF Matrix
	for (auto& nd_m : model_nodes.nodeMap)
	{
		// Get the node data
		node_store nd = nd_m.second;
		int nd_map = nodeid_map[nd.node_id]; // get the ordered map of the node ID

		if (model_constarints.constraintMap.find(nd.node_id) != model_constarints.constraintMap.end())
		{
			// Nodes have point mass
			constraint_data cd = model_constarints.constraintMap.at(nd.node_id);

			if (cd.constraint_type == 0)
			{
				// Fixed End

				globalDOFMatrix.coeffRef((nd_map * 3) + 0, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 1, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 2, 0) = 0.0;
			}
			else if (cd.constraint_type == 1)
			{
				// Fixed Roller end

				globalDOFMatrix.coeffRef((nd_map * 3) + 0, 0) = 1.0; // X is free to move
				globalDOFMatrix.coeffRef((nd_map * 3) + 1, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 2, 0) = 0.0;

				reducedDOF = reducedDOF + 1;
			}
			else if (cd.constraint_type == 2)
			{
				// Pin End

				globalDOFMatrix.coeffRef((nd_map * 3) + 0, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 1, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 2, 0) = 1.0; // XY Rotation is free to move

				reducedDOF = reducedDOF + 1;
			}
			else if (cd.constraint_type == 3)
			{
				// Pin Roller End

				globalDOFMatrix.coeffRef((nd_map * 3) + 0, 0) = 1.0; // X is free to move
				globalDOFMatrix.coeffRef((nd_map * 3) + 1, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 3) + 2, 0) = 1.0; // XY Rotation is free to move

				reducedDOF = reducedDOF + 2;
			}
		}
		else
		{
			// Nodes doesnt have Constraint
			globalDOFMatrix.coeffRef((nd_map * 3) + 0, 0) = 1.0;
			globalDOFMatrix.coeffRef((nd_map * 3) + 1, 0) = 1.0;
			globalDOFMatrix.coeffRef((nd_map * 3) + 2, 0) = 1.0;

			reducedDOF = reducedDOF + 3;
		}
	}
}

void pulse_analysis_solver::get_reduced_global_matrices(Eigen::MatrixXd& reduced_globalStiffnessMatrix,
	Eigen::MatrixXd& reduced_globalMassMatrix,
	const Eigen::MatrixXd& globalStiffnessMatrix,
	const Eigen::MatrixXd& globalMassMatrix,
	const Eigen::MatrixXd& globalDOFMatrix,
	const int& numDOF,
	std::ofstream& output_file)
{
	// Curtailment of Global stiffness and Global force matrix based on DOF
	// Get the reduced global stiffness matrix
	int r = 0;
	int s = 0;

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
			s = 0;
			for (int j = 0; j < numDOF; j++)
			{
				if (globalDOFMatrix(j, 0) == 0)
				{
					// constrained column index, so skip
					continue;
				}
				else
				{
					// Get the reduced matrices
					reduced_globalMassMatrix.coeffRef(r, s) = globalMassMatrix.coeffRef(i, j);
					reduced_globalStiffnessMatrix.coeffRef(r, s) = globalStiffnessMatrix.coeffRef(i, j);
					s++;
				}
			}

			// reduced_globalLoadMatrix.coeffRef(r, 0) = globalLoadMatrix(i, 0);
			r++;
		}
	}

	if (print_matrix == true)
	{
		// Print the Reduced Global Stiffness, Global Mass and Global Force matrix
		output_file << "Reduced Global Stiffness Matrix" << std::endl;
		output_file << reduced_globalStiffnessMatrix << std::endl;
		output_file << std::endl;

		output_file << "Reduced Global Mass Matrix" << std::endl;
		output_file << reduced_globalMassMatrix << std::endl;
		output_file << std::endl;
	}
}

void pulse_analysis_solver::get_modal_matrices(Eigen::VectorXd& modalMass,
	Eigen::VectorXd& modalStiff,
	const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
	const Eigen::MatrixXd& reduced_globalMassMatrix,
	const Eigen::MatrixXd& reduced_globalStiffnessMatrix,
	const int& reducedDOF,
	std::ofstream& output_file)
{
	// Get the modal matrices
	Eigen::MatrixXd modalMassMatrix(reducedDOF, reducedDOF);
	modalMassMatrix.setZero();

	modalMassMatrix = reduced_eigenVectorsMatrix.transpose() * reduced_globalMassMatrix * reduced_eigenVectorsMatrix;

	// Create modal stiffness matrix
	Eigen::MatrixXd modalStiffMatrix(reducedDOF, reducedDOF);
	modalStiffMatrix.setZero();

	modalStiffMatrix = reduced_eigenVectorsMatrix.transpose() * reduced_globalStiffnessMatrix * reduced_eigenVectorsMatrix;

	// Create the modal vectors
	modalMass = modalMassMatrix.diagonal();
	modalStiff = modalStiffMatrix.diagonal();

	if (print_matrix == true)
	{
		// Print the Modal Mass Matrix
		output_file << "Modal Mass Matrix" << std::endl;
		output_file << modalMassMatrix << std::endl;
		output_file << std::endl;

		// Print the Modal Stiffness matrix
		output_file << "Modal Stiffness Matrix" << std::endl;
		output_file << modalStiffMatrix << std::endl;
		output_file << std::endl;
	}
}

void pulse_analysis_solver::get_reduced_modal_vector_matrix(Eigen::MatrixXd& reduced_eigenVectorsMatrix,
	const modal_analysis_result_store& modal_results,
	int& reducedDOF,
	std::ofstream& output_file)
{

	// Get the global eigen vectors matrix
	for (int i = 0; i < reducedDOF; i++)
	{
		std::vector<double> eigen_vector_column = modal_results.eigen_vectors_reduced.at(i);
		for (int j = 0; j < reducedDOF; j++)
		{
			reduced_eigenVectorsMatrix.coeffRef(j, i) = eigen_vector_column[j];
		}
	}

	if (print_matrix == true)
	{
		// Print the Reduced Global Eigen Vector matrix
		output_file << "Reduced Global EigenVector Matrix" << std::endl;
		output_file << reduced_eigenVectorsMatrix << std::endl;
		output_file << std::endl;
	}
}


void pulse_analysis_solver::get_globalSupportInclinationMatrix(Eigen::MatrixXd& globalSupportInclinationMatrix,
	const nodes_list_store& model_nodes,
	const nodeconstraint_list_store& model_constarints,
	const int& numDOF,
	std::ofstream& output_file)
{
	// Create the global support inclination matrix
	int node_id = 0;

	// Transform the Nodal results with support inclination
	double constraint_angle_rad = 0.0;
	double support_Lcos = 0.0;
	double support_Msin = 0.0;

	for (auto& nd_m : model_nodes.nodeMap)
	{
		node_id = nd_m.first;
		int matrix_index = nodeid_map[node_id];


		if (model_constarints.constraintMap.find(node_id) == model_constarints.constraintMap.end())
		{
			// No constraint is in this node
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 0) = 1.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 1) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 2) = 0.0;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 0) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 1) = 1.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 2) = 0.0;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 0) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 1) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 2) = 1.0;
		}
		else
		{
			// Constraint present in this node
			constraint_angle_rad = (model_constarints.constraintMap.at(node_id).constraint_angle - 90.0f) * (m_pi / 180.0f); // Constrint angle in radians
			support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
			support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

			// Pin or Roller Support
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 0) = support_Lcos;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 1) = -1.0 * support_Msin;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 0, (matrix_index * 3) + 2) = 0.0;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 0) = support_Msin;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 1) = support_Lcos;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 1, (matrix_index * 3) + 2) = 0.0;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 0) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 1) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 3) + 2, (matrix_index * 3) + 2) = 1.0;
		}
	}

	if (print_matrix == true)
	{
		// Print the Global support inclination matrix
		output_file << "Global Support Inclination Matrix" << std::endl;
		output_file << globalSupportInclinationMatrix << std::endl;
		output_file << std::endl;
	}
}

void pulse_analysis_solver::create_pulse_load_matrices(pulse_load_data& pulse_loads,
	const load_data& ld,
	const elementline_list_store& model_lineelements,
	const Eigen::MatrixXd& globalDOFMatrix,
	const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
	const int& numDOF,
	const int& reducedDOF)
{
	// Create the global load amplitude matrix
	// Extract the line in which the load is applied
	elementline_store ln = model_lineelements.elementlineMap.at(ld.node_id);

	// Get the Matrix row ID
	int sn_id = nodeid_map[ln.startNode->node_id]; // get the ordered map of the start node ID
	int en_id = nodeid_map[ln.endNode->node_id]; // get the ordered map of the end node ID

	// Create a matrix for element load matrix
	Eigen::MatrixXd elementLoadMatrix(6, 1);
	elementLoadMatrix.setZero();

	get_element_load_matrix(elementLoadMatrix, ln, ld);

	// global load matrix
	Eigen::MatrixXd globalLoadamplMatrix(numDOF, 1);
	globalLoadamplMatrix.setZero();

	globalLoadamplMatrix.block<3, 1>(sn_id * 3, 0) += elementLoadMatrix.block<3, 1>(0, 0);
	globalLoadamplMatrix.block<3, 1>(en_id * 3, 0) += elementLoadMatrix.block<3, 1>(3, 0);

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

void pulse_analysis_solver::get_element_load_matrix(Eigen::MatrixXd& elementLoadMatrix,
	const elementline_store& ln,
	const load_data& ld)
{
	// Compute the differences in x and y coordinates
	double dx = ln.endNode->node_pt.x - ln.startNode->node_pt.x;
	double dy = -1.0 * (ln.endNode->node_pt.y - ln.startNode->node_pt.y);

	// Compute the length of the frame element
	double eLength = std::sqrt((dx * dx) + (dy * dy));

	// Compute the direction cosines
	double Lcos = (dx / eLength);
	double Msin = (dy / eLength);

	// Create element load mass matrix
		// Load's line id is equal to this line id
	double load_val = ld.load_value; // Load value
	double load_angle_rad = ld.load_angle * (m_pi / 180.0f);

	double f_x = load_val * std::cos(load_angle_rad);
	double f_y = load_val * std::sin(load_angle_rad);

	double loadHorizontal = (-f_x * Lcos) + (f_y * Msin);
	double loadVertical = (f_x * Msin) + (f_y * Lcos);

	double a = eLength *0.0;
	double b = eLength * 0.0;

	// Load parameter
	double phi_i = (loadVertical * b * (std::pow(eLength, 2) - std::pow(b, 2)) / (6 * eLength));
	double phi_j = (loadVertical * a * (std::pow(eLength, 2) - std::pow(a, 2)) / (6 * eLength));

	// Axial Load
	double Tai = loadHorizontal *0.0;
	double Taj = loadHorizontal * (1 - 0.0);

	// End Moments
	double Tmi = (((4 * phi_i) - (2 * phi_j)) / eLength);
	double Tmj = (((2 * phi_i) - (4 * phi_j)) / eLength);

	// Vertical load
	double Tfi = (((loadVertical * b) + Tmi + Tmj) / eLength);
	double Tfj = (((loadVertical * a) - Tmi - Tmj) / eLength);

	// Add to the element load matrix
	elementLoadMatrix.coeffRef(0, 0) += Tai;
	elementLoadMatrix.coeffRef(1, 0) += (-Tfi);
	elementLoadMatrix.coeffRef(2, 0) += (-Tmi);

	elementLoadMatrix.coeffRef(3, 0) += Taj;
	elementLoadMatrix.coeffRef(4, 0) += (-Tfj);
	elementLoadMatrix.coeffRef(5, 0) += (-Tmj);

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
		for (auto& h_ln : ln_m.second.hermite_line_data)
		{
			//get all the two points
			// Point 1 displacement
			for (auto& pt1_m : h_ln.pt1_modal_displ)
			{
				glm::vec2 pt1 = pt1_m.second; //get the end point 1 displacement
				double displ1 = std::sqrt(std::pow(pt1.x, 2) + std::pow(pt1.y, 2));

				if (displ1 > maximum_displacement)
				{
					maximum_displacement = displ1;
				}
			}

			// Point 2 displacement
			for (auto& pt2_m : h_ln.pt2_modal_displ)
			{
				glm::vec2 pt2 = pt2_m.second; //get the end point 2 displacement
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