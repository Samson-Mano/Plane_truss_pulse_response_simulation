#include "modal_analysis_solver.h"

modal_analysis_solver::modal_analysis_solver()
{
	// Empty constructor
}

modal_analysis_solver::~modal_analysis_solver()
{
	// Empty destructor
}

void modal_analysis_solver::modal_analysis_start(const nodes_list_store& model_nodes,
	const elementline_list_store& model_lineelements,
	const nodeconstraint_list_store& model_constarints,
	const nodepointmass_list_store& model_ptmass,
	const std::unordered_map<int, material_data>& material_list,
	const bool& is_include_consistent_mass_matrix,
	modal_analysis_result_store& modal_results,
	modal_nodes_list_store& modal_result_nodes,
	modal_elementline_list_store& modal_result_lineelements,
	bool& is_modal_analysis_complete)
{
	// Main solver call
	is_modal_analysis_complete = false;

	// Check the model
	// Number of nodes
	if (model_nodes.node_count == 0)
	{
		return;
	}

	// Number of elements
	if (model_lineelements.elementline_count == 0)
	{
		return;
	}

	// check whether masses are valid
	if (is_include_consistent_mass_matrix == false)
	{
		// No consistent mass
		// Check whether there are point masses
		if (model_ptmass.ptmass_count == 0)
		{
			// No Point mass, no consistent mass
			return;
		}
	}

	//____________________________________________

	// Create a node ID map (to create a nodes as ordered and numbered from 0,1,2...n)
	int i = 0;
	for (auto& nd : model_nodes.nodeMap)
	{
		nodeid_map[nd.first] = i;
		i++;
	}

	// Create a file to keep track of matrices
	std::ofstream output_file;
	output_file.open("modal_analysis_results.txt");

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
	// Solve generalized Eigen value matrix using Cholesky decomposition
	// Compute the Cholesky decomposition of Global Mass matrix
	// Generalized Symmetric Definite Eigenproblems 

	Eigen::LLT<Eigen::MatrixXd> llt;
	llt.compute(reduced_globalMassMatrix);


	if (llt.info() != Eigen::Success) {
		// Cholesky decomposition failed
		output_file << "Cholesky decomposition failed !!!!" << std::endl;
	}

	// Get the lower triangular matrix L
	Eigen::MatrixXd L_matrix = llt.matrixL();

	if (print_matrix == true)
	{
		// Print the Cholesky Decomposition L - Matrix
		output_file << "Cholesky Decomposition L - Matrix" << std::endl;
		output_file << L_matrix << std::endl;
		output_file << std::endl;
	}

	// Get the L^-1 inverse of L-matrix Lower triangular matrix
	Eigen::MatrixXd L_inv_matrix = L_matrix.inverse();

	if (print_matrix == true)
	{
		// Print the Inverse L - Matrix
		output_file << "L Inverse Matrix" << std::endl;
		output_file << L_inv_matrix << std::endl;
		output_file << std::endl;
	}

	//____________________________________________________________________________________________________________________
	//  Find the eigen value & eigen vector of eigen value problem Z_matrix
	//  Z_matrix = L_inv_matrix * stiff_matrix * L_inv_matrix^T

	Eigen::MatrixXd Z_matrix = L_inv_matrix * reduced_globalStiffnessMatrix * L_inv_matrix.transpose();

	if (print_matrix == true)
	{
		// Print the Inverse L - Matrix
		output_file << "Z Matrix" << std::endl;
		output_file << Z_matrix << std::endl;
		output_file << std::endl;
	}

	// Compute the eigenvalues and eigenvectors
	Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(Z_matrix);

	if (eigenSolver.info() != Eigen::Success) {
		// Eigenvalue problem failed to converge
		output_file << "Eigenvalue problem failed to converge !!!!! " << std::endl;
		return;
	}

	// Get the eigenvalues and eigenvectors
	Eigen::VectorXd eigenvalues = eigenSolver.eigenvalues().real(); // Real part of eigenvalues
	Eigen::MatrixXd eigenvectors_reduced = L_inv_matrix.transpose() * eigenSolver.eigenvectors().real(); // Real part of eigenvectors

	// sort the eigen value and eigen vector (ascending)
	sort_eigen_values_vectors(eigenvalues, eigenvectors_reduced, reducedDOF);

	// Normailize eigen vectors
	normalize_eigen_vectors(eigenvectors_reduced, reducedDOF);

	//____________________________________________________________________________________________________________________

	if (print_matrix == true)
	{
		// Eigenvalue problem failed to converge
		output_file << "Modal Analysis Success !!!!! " << std::endl;

		// Print the Eigen values
		output_file << "Eigen Values " << std::endl;
		output_file << eigenvalues << std::endl;
		output_file << std::endl;
	}

	// Convert the reduced eigenvectors to eigen vectors for the whole model (including the nodes with supports)
	Eigen::MatrixXd eigenvectors(numDOF, reducedDOF);
	eigenvectors.setZero();

	get_global_modal_vector_matrix(eigenvectors, eigenvectors_reduced, globalDOFMatrix, numDOF, reducedDOF, output_file);

	//____________________________________________________________________________________________________________________
	// Store the results

	is_modal_analysis_complete = true;

	// Clear the modal results
	modal_results.clear_data();

	// Add the eigen values and eigen vectors
	for (int i = 0; i < reducedDOF; i++)
	{
		double eigen_val = eigenvalues(i);
		std::vector<double> eigen_vec; // Eigen vectors of all nodes (including constrainded)
		std::vector<double> eigen_vec_reduced; // Eigen vectors of nodes (reduced)

		for (int j = 0; j < numDOF; j++)
		{
			eigen_vec.push_back(eigenvectors(j, i));
		}

		for (int j = 0; j < reducedDOF; j++)
		{
			eigen_vec_reduced.push_back(eigenvectors_reduced(j, i));
		}

		modal_results.add_eigen_data(i, eigen_val, eigen_vec, eigen_vec_reduced);
	}
	modal_results.add_node_map(nodeid_map);

	//____________________________________________________________________________________________
	// Convert the modal result to string  
	int num_of_rigidbody_modes = get_number_of_rigid_body_modes(model_nodes.node_count,globalDOFMatrix);
	std::vector<std::string> mode_result_str;
	int k = 0;

	for (int i = 0; i < modal_results.number_of_modes; i++)
	{
		// Convert the angular frequency^2 to Frequency Hz
		double freq = 0.0;

		if (modal_results.eigen_values[i] > 0.0001)
		{
			// Result is not zero
			freq = std::sqrt(modal_results.eigen_values[i]) / (2.0 * m_pi);
		}
	
		// Modal results
		std::stringstream ss;
		ss << std::fixed << std::setprecision(2) << freq;

		if (freq == 0)   //(i < num_of_rigidbody_modes)
		{
			mode_result_str.push_back("Rigid body mode " + std::to_string(i + 1) + " = " + ss.str() + " Hz");
		}
		else
		{
			mode_result_str.push_back("Mode " + std::to_string(k + 1) + " = " + ss.str() + " Hz");
			k++;
		}
	}

	modal_results.mode_result_str = mode_result_str;

	//_____________________________________________________________________________________________

	// Add the modal analysis results to node & element
	// Clear the modal node and modal element results
	modal_result_nodes.clear_data();
	modal_result_lineelements.clear_data();

	map_modal_analysis_results(model_nodes,
		model_lineelements,
		model_constarints,
		modal_results,
		modal_result_nodes,
		modal_result_lineelements,
		output_file);

	//____________________________________________________________________________________________________________________

	output_file.close();

}

void modal_analysis_solver::get_global_stiffness_matrix(Eigen::MatrixXd& globalStiffnessMatrix,
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

	if (print_matrix == true)
	{
		// Print the Global Stiffness matrix
		output_file << "Global Stiffness Matrix" << std::endl;
		output_file << globalStiffnessMatrix << std::endl;
		output_file << std::endl;
	}
}

void modal_analysis_solver::get_element_stiffness_matrix(Eigen::MatrixXd& elementStiffnessMatrix,
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

void modal_analysis_solver::get_global_pointmass_matrix(Eigen::MatrixXd& globalPointMassMatrix,
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
			globalPointMassMatrix((nd_map * 3) + 2, (nd_map * 3) + 2) = ptm.ptmass_xy;
		}
		else
		{
			// Nodes doesnt have point mass
			globalPointMassMatrix((nd_map * 3) + 0, (nd_map * 3) + 0) = 0.0;
			globalPointMassMatrix((nd_map * 3) + 1, (nd_map * 3) + 1) = 0.0;
			globalPointMassMatrix((nd_map * 3) + 2, (nd_map * 3) + 2) = 0.0;
		}
	}

	if (print_matrix == true)
	{
		// Print the Global Force matrix
		output_file << "Global Point Mass Matrix" << std::endl;
		output_file << globalPointMassMatrix << std::endl;
		output_file << std::endl;
	}
}




void modal_analysis_solver::get_global_consistentmass_matrix(Eigen::MatrixXd& globalConsistentMassMatrix,
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

	if (print_matrix == true)
	{
		// Print the Global Stiffness matrix
		output_file << "Global Consistent Mass Matrix" << std::endl;
		output_file << globalConsistentMassMatrix << std::endl;
		output_file << std::endl;
	}
}

void modal_analysis_solver::get_element_consistentmass_matrix(Eigen::MatrixXd& elementConsistentMassMatrix,
	const elementline_store& ln,
	const material_data& elementline_material,
	const nodeconstraint_list_store& model_constarints,
	std::ofstream& output_file)
{
	// Create element consistent mass matrix
	// Get element stiffness matrix
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


void modal_analysis_solver::get_global_dof_matrix(Eigen::MatrixXd& globalDOFMatrix,
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

	if (print_matrix == true)
	{
		// Print the Global Force matrix
		output_file << "Global DOF Matrix" << std::endl;
		output_file << globalDOFMatrix << std::endl;
		output_file << std::endl;
	}
}

void modal_analysis_solver::get_reduced_global_matrices(Eigen::MatrixXd& reduced_globalStiffnessMatrix,
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
			r++;
		}
	}

	if (print_matrix == true)
	{
		// Print the Reduced Global Stiffness and Reduced Global Force matrix
		output_file << "Reduced Global Stiffness Matrix" << std::endl;
		output_file << reduced_globalStiffnessMatrix << std::endl;
		output_file << std::endl;

		output_file << "Reduced Global Mass Matrix" << std::endl;
		output_file << reduced_globalMassMatrix << std::endl;
		output_file << std::endl;
	}
}

void modal_analysis_solver::sort_eigen_values_vectors(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors, const int& m_size)
{
	int p = 0;
	int q = 0;
	int i = 0;

	double swap_temp = 0.0;

	// sort the eigen value and eigen vector (ascending)
	for (p = 0; p < m_size; p++)
	{
		for (q = p + 1; q < m_size; q++)
		{
			if (eigenvalues(p) > eigenvalues(q))
			{
				swap_temp = eigenvalues(p);
				eigenvalues(p) = eigenvalues(q);
				eigenvalues(q) = swap_temp;
				for (i = 0; i < m_size; i++)
				{
					swap_temp = eigenvectors(i, p);
					eigenvectors(i, p) = eigenvectors(i, q);
					eigenvectors(i, q) = swap_temp;
				}
			}
		}
	}
}

void modal_analysis_solver::normalize_eigen_vectors(Eigen::MatrixXd& eigenvectors,
	const int& m_size)
{
	// Normalize eigen vectors
	int p = 0;
	int q = 0;

	// loop throught each column
	for (p = 0; p < m_size; p++)
	{
		double max_modal_vector = 0.0;

		// Loop through each row
		for (q = 0; q < m_size; q++)
		{
			if (std::abs(eigenvectors(q, p)) > max_modal_vector)
			{
				// Max modal vector in the column (for particular mode)
				max_modal_vector = std::abs(eigenvectors(q, p));
			}
		}

		// Normalize the column using maximum modal vector
		for (q = 0; q < m_size; q++)
		{
			eigenvectors(q, p) = eigenvectors(q, p) / max_modal_vector;

			// Round the eigen vectors to 6 digit precision after normalizing
			eigenvectors(q, p) = std::round(eigenvectors(q, p) * 1000000) / 1000000;
		}
	}
}

Eigen::MatrixXd modal_analysis_solver::convert_vector_to_1Dmatrix(const std::vector<double>& vec)
{
	// Convert Vector to 1D Column matrix
	int vec_size = static_cast<int>(vec.size());

	Eigen::MatrixXd mat(vec_size, 1);

	for (int i = 0; i < vec_size; ++i)
	{
		mat.coeffRef(i, 0) = vec.at(i);
	}
	return mat;
}

int modal_analysis_solver::get_number_of_rigid_body_modes(int num_of_nodes,
	const Eigen::MatrixXd& globalDOFMatrix)
{
	// Set all three DOF free 
	bool x_free = true;
	bool y_free = true;
	bool xy_free = true;

	// get the number of rigid body modes
	for (int i = 0; i < num_of_nodes; i++)
	{
		// x DOF
		if (globalDOFMatrix(((i * 3) + 0), 0) == 0)
		{
			x_free = false;
		}

		// y DOF
		if (globalDOFMatrix(((i * 3) + 1), 0) == 0)
		{
			y_free = false;
		}


		// z DOF
		if (globalDOFMatrix(((i * 3) + 2), 0) == 0)
		{
			xy_free = false;
		}
	}

	// Get the total free
	int rigid_body_modes = 0;

	if (x_free == true)
	{
		// No constraint at x direction
		rigid_body_modes++;
	}

	if (y_free == true)
	{
		// No constraint at y direction
		rigid_body_modes++;
	}

	if (xy_free == true)
	{
		// No constraint at xy rotation
		// rigid_body_modes++;
	}
	

	return rigid_body_modes;
}


void modal_analysis_solver::get_global_modal_vector_matrix(Eigen::MatrixXd& eigenvectors,
	const Eigen::MatrixXd& eigenvectors_reduced,
	const Eigen::MatrixXd& globalDOFMatrix,
	const int& numDOF,
	int& reducedDOF,
	std::ofstream& output_file)
{
	// Global eigen vector Matrix
	// Loop throug the Degree of freedom of indices

	// J loops through number of modes (along the column)
	for (int j = 0; j < reducedDOF; j++)
	{
		int s = 0;
		// i loops through the number of nodes (along the row)
		for (int i = 0; i < numDOF; i++)
		{
			if (globalDOFMatrix(i, 0) == 0)
			{
				// constrained row index, so Displacement is Zero
				eigenvectors.coeffRef(i, j) = 0;
			}
			else
			{
				// Un constrained row index, so Displacement is Zero
				eigenvectors.coeffRef(i, j) = eigenvectors_reduced.coeffRef(s, j);
				s++;
			}
		}
	}

	if (print_matrix == true)
	{
		// Print the Global Displacement matrix
		output_file << "Global Eigen Vector Matrix" << std::endl;
		output_file << eigenvectors << std::endl;
		output_file << std::endl;
	}
}

void modal_analysis_solver::map_modal_analysis_results(const nodes_list_store& model_nodes,
	const elementline_list_store& model_lineelements,
	const nodeconstraint_list_store& model_constarints,
	const modal_analysis_result_store& modal_results,
	modal_nodes_list_store& modal_result_nodes,
	modal_elementline_list_store& modal_result_lineelements,
	std::ofstream& output_file)
{
	// Map the modal analysis results to modal members
	// Create a global support inclination matrix
	int numDOF = model_nodes.node_count * 3; // Number of degrees of freedom (2 DOFs per node)
	int node_id = 0;

	Eigen::MatrixXd globalSupportInclinationMatrix(numDOF, numDOF);
	globalSupportInclinationMatrix.setZero();

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

	//___________________________________________________________________________________________________________
	// Transform the  Eigen Vector with support inclination
	Eigen::MatrixXd globalEigenVector_transformed(numDOF, modal_results.number_of_modes);

	for (int i = 0; i < modal_results.number_of_modes; i++)
	{
		// Get the modal vector matrix for this particular mode
		Eigen::MatrixXd globalEigenVector = convert_vector_to_1Dmatrix(modal_results.eigen_vectors.at(i));

		// Transform the global Modal Vector w.r.t support inclination
		globalEigenVector_transformed.col(i) = globalSupportInclinationMatrix * globalEigenVector;
	}

	//___________________________________________________________________________________________________________
	// Add to the result nodes
	for (auto& nd_m : model_nodes.nodeMap)
	{
		node_id = nd_m.first;
		int matrix_index = nodeid_map[node_id];

		// Modal analysis results
		std::unordered_map<int, glm::vec3> node_modal_displ;

		for (int i = 0; i < modal_results.number_of_modes; i++)
		{
			// get the appropriate modal displacement of this particular node
			glm::vec3 modal_displ = glm::vec3(globalEigenVector_transformed((matrix_index * 3) + 0, i),
				globalEigenVector_transformed((matrix_index * 3) + 1, i),
				globalEigenVector_transformed((matrix_index * 3) + 2, i));

			// add to modal result of this node
			node_modal_displ.insert({ i,modal_displ });
		}

		// Create the modal analysis result node
		glm::vec2 node_pt = model_nodes.nodeMap.at(node_id).node_pt;
		modal_result_nodes.add_result_node(node_id, node_pt, node_modal_displ);
	}

	// Add the modal line element result
	for (auto& ln_m : model_lineelements.elementlineMap)
	{
		elementline_store ln = ln_m.second;

		modal_result_lineelements.add_modal_elementline(ln.line_id,
			&modal_result_nodes.modal_nodeMap[ln.startNode->node_id],
			&modal_result_nodes.modal_nodeMap[ln.endNode->node_id]);
	}

	// Find the maximum displacement for individual modes
	std::unordered_map<int, double> max_node_displ;
	std::unordered_map<int, double> min_node_displ;

	for (int i = 0; i < modal_results.number_of_modes; i++)
	{
		// Go through all the line points at this particular mode
		double max_displ = 0.0;
		double min_displ = INT32_MAX;

		for (auto& ln_m : modal_result_lineelements.modal_elementlineMap)
		{
			modal_elementline_store ln = ln_m.second;

			// Get all the modal line data
			for (auto& mln_hld : ln.hermite_line_data)
			{
				// loop through all the mode line end points
				glm::vec2 mln_pt1 = mln_hld.pt1_modal_displ[i];
				glm::vec2 mln_pt2 = mln_hld.pt2_modal_displ[i];

				// Check all the point
				// Maximum
				if (max_displ < std::abs(std::sqrt((mln_pt1.x * mln_pt1.x) + (mln_pt1.y * mln_pt1.y))))
				{
					// pt1
					max_displ = std::abs(std::sqrt((mln_pt1.x * mln_pt1.x) + (mln_pt1.y * mln_pt1.y)));
				}

				if (max_displ < std::abs(std::sqrt((mln_pt2.x * mln_pt2.x) + (mln_pt2.y * mln_pt2.y))))
				{
					// pt2
					max_displ = std::abs(std::sqrt((mln_pt2.x * mln_pt2.x) + (mln_pt2.y * mln_pt2.y)));
				}
				//____________________________________________________________________________________________
				// Minimum
				if (min_displ > std::abs(std::sqrt((mln_pt1.x * mln_pt1.x) + (mln_pt1.y * mln_pt1.y))))
				{
					// pt1
					min_displ = std::abs(std::sqrt((mln_pt1.x * mln_pt1.x) + (mln_pt1.y * mln_pt1.y)));
				}

				if (min_displ > std::abs(std::sqrt((mln_pt2.x * mln_pt2.x) + (mln_pt2.y * mln_pt2.y))))
				{
					// pt2
					min_displ = std::abs(std::sqrt((mln_pt2.x * mln_pt2.x) + (mln_pt2.y * mln_pt2.y)));
				}
			}

		}

		// Add to the maximum displacement of this mode
		max_node_displ.insert({ i,max_displ });
		min_node_displ.insert({ i,min_displ });
	}

	// Set the maximum modal displacement
	modal_result_nodes.max_node_displ.clear();
	modal_result_nodes.max_node_displ = max_node_displ;
	modal_result_nodes.min_node_displ = min_node_displ;

	modal_result_lineelements.max_node_displ.clear();
	modal_result_lineelements.max_node_displ = max_node_displ;
	modal_result_lineelements.min_node_displ = min_node_displ;
}

