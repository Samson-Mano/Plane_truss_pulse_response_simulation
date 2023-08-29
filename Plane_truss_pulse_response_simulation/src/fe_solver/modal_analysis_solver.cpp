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
	Eigen::initParallel();  // Initialize Eigen's thread pool

	stopwatch.start();
	std::stringstream stopwatch_elapsed_str;
	stopwatch_elapsed_str << std::fixed << std::setprecision(6);

	std::cout << "Modal analysis started" << std::endl;

	// Create a node ID map (to create a nodes as ordered and numbered from 0,1,2...n)
	int i = 0;
	for (auto& nd : model_nodes.nodeMap)
	{
		nodeid_map[nd.first] = i;
		i++;
	}

	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Node maping completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Create a file to keep track of matrices
	std::ofstream output_file;
	output_file.open("modal_analysis_results.txt");

	//____________________________________________________________________________________________________________________
	numDOF = model_nodes.node_count * 2; // Number of degrees of freedom (2 DOFs per node (2 translation)

	// Global Stiffness Matrix
	globalStiffnessMatrix.resize(numDOF, numDOF);
	globalStiffnessMatrix.setZero();

	get_global_stiffness_matrix(globalStiffnessMatrix,
		model_lineelements,
		material_list,
		model_constarints,
		output_file);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Global stiffness matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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
	// Global Mass Matrix
	globalMassMatrix.resize(numDOF, numDOF);
	globalMassMatrix.setZero();

	globalMassMatrix = globalPointMassMatrix + globalConsistentMassMatrix;

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Global mass matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	//____________________________________________________________________________________________________________________
	// Global DOF Mass Matrix
	globalDOFMatrix.resize(numDOF, 1);
	globalDOFMatrix.setZero();

	// Determine the size of the reduced stiffness matrix based on the number of unconstrained degrees of freedom
	reducedDOF = 0;

	get_global_dof_matrix(globalDOFMatrix,
		model_nodes,
		model_constarints,
		reducedDOF,
		output_file);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Global DOF matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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


	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Global stiffness/mass matrces are reduced at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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


	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Cholesky decomposition L-matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Get the L^-1 inverse of L-matrix Lower triangular matrix
	Eigen::MatrixXd L_inv_matrix = L_matrix.inverse();

	if (print_matrix == true)
	{
		// Print the Inverse L - Matrix
		output_file << "L Inverse Matrix" << std::endl;
		output_file << L_inv_matrix << std::endl;
		output_file << std::endl;
	}

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Inverse of lower triangle matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Generalized Eigen value problem Z-matrix created at " << stopwatch_elapsed_str.str() << " secs" << std::endl;
	std::cout << "Size of the Z-matrix is " << reducedDOF << " x " << reducedDOF << std::endl;


	// Compute the eigenvalues and eigenvectors
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(Z_matrix);

	if (eigenSolver.info() != Eigen::Success) {
		// Eigenvalue problem failed to converge
		output_file << "Eigenvalue problem failed to converge !!!!! " << std::endl;
		return;
	}

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Eigen value problem solved at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Get the eigenvalues and eigenvectors
	Eigen::VectorXd eigenvalues = eigenSolver.eigenvalues(); // Eigenvalues
	Eigen::MatrixXd eigenvectors_reduced = L_inv_matrix.transpose() * eigenSolver.eigenvectors(); // Eigenvectors

	// sort the eigen value and eigen vector (ascending)
	sort_eigen_values_vectors(eigenvalues, eigenvectors_reduced, reducedDOF);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Eigen values and Eigen vectors are sorted at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Normailize eigen vectors
	normalize_eigen_vectors(eigenvectors_reduced, reducedDOF);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Eigen vectors are normalized at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Eigen vectors globalized at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

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
	int num_of_rigidbody_modes = get_number_of_rigid_body_modes(model_nodes.node_count, globalDOFMatrix);
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

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Modal results are storage completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	//_____________________________________________________________________________________________
	// Create global support inclination matrix
	globalSupportInclinationMatrix.resize(numDOF, numDOF);
	globalSupportInclinationMatrix.setZero();

	get_globalSupportInclinationMatrix(globalSupportInclinationMatrix,
		model_nodes,
		model_constarints,
		numDOF,
		output_file);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Global support inclination matrix completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	//_____________________________________________________________________________________________

	// Add the modal analysis results to node & element
	// Clear the modal node and modal element results
	modal_result_nodes.clear_data();
	modal_result_lineelements.clear_data();

	map_modal_analysis_results(model_nodes,
		model_lineelements,
		model_constarints,
		modal_results,
		globalSupportInclinationMatrix,
		modal_result_nodes,
		modal_result_lineelements,
		output_file);


	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Modal analysis results maped to nodes and elements at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	//____________________________________________________________________________________________________________________
	// Modal Decomposition

	// Reduced Eigen Vectors matrix
	reduced_eigenVectorsMatrix.resize(reducedDOF, reducedDOF);
	reduced_eigenVectorsMatrix.setZero();

	get_reduced_modal_vector_matrix(reduced_eigenVectorsMatrix, modal_results, reducedDOF, output_file);

	//____________________________________________________________________________________________________________________
	// Create modal matrices
	modalMass.resize(reducedDOF);
	modalStiff.resize(reducedDOF);

	get_modal_matrices(modalMass,
		modalStiff,
		reduced_eigenVectorsMatrix,
		reduced_globalMassMatrix,
		reduced_globalStiffnessMatrix,
		reducedDOF,
		output_file);


	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Modal mass and stiffness storage completed at " << stopwatch_elapsed_str.str() << " secs" << std::endl;
	std::cout << "Modal analysis complete " << std::endl;

	//____________________________________________________________________________________________________________________
	stopwatch.stop();

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
		Eigen::MatrixXd elementStiffnessMatrix(4, 4);
		elementStiffnessMatrix.setZero();

		get_element_stiffness_matrix(elementStiffnessMatrix, ln, elementline_material, model_constarints, output_file);

		// Get the Node ID
		int sn_id = nodeid_map[ln.startNode->node_id]; // get the ordered map of the start node ID
		int en_id = nodeid_map[ln.endNode->node_id]; // get the ordered map of the end node ID

		globalStiffnessMatrix.block<2, 2>(sn_id * 2, sn_id * 2) += elementStiffnessMatrix.block<2, 2>(0, 0);
		globalStiffnessMatrix.block<2, 2>(sn_id * 2, en_id * 2) += elementStiffnessMatrix.block<2, 2>(0, 2);
		globalStiffnessMatrix.block<2, 2>(en_id * 2, sn_id * 2) += elementStiffnessMatrix.block<2, 2>(2, 0);
		globalStiffnessMatrix.block<2, 2>(en_id * 2, en_id * 2) += elementStiffnessMatrix.block<2, 2>(2, 2);
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
	Eigen::MatrixXd L_transformation_matrix(2, 4);
	L_transformation_matrix.setZero();

	L_transformation_matrix.row(0) = Eigen::RowVectorXd({ {Lcos, Msin, 0.0, 0.0 } });
	L_transformation_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 0.0,  Lcos, Msin} });

	//_________________________________________________________
	// Local element stiffness matrix
	Eigen::MatrixXd local_element_stiffness_matrix(2, 2);
	local_element_stiffness_matrix.setZero();

	double k1 = (elementline_material.youngs_mod * elementline_material.cs_area) / eLength;

	local_element_stiffness_matrix.row(0) = Eigen::RowVectorXd({ {		 k1, -1.0 * k1 } });
	local_element_stiffness_matrix.row(1) = Eigen::RowVectorXd({ {-1.0 * k1,		k1 } });

	//_________________________________________________________
	// Transformed element stiffness matrix
	Eigen::MatrixXd e_stiffness_matrix(4, 4);
	e_stiffness_matrix.setZero();

	e_stiffness_matrix = L_transformation_matrix.transpose() * local_element_stiffness_matrix * L_transformation_matrix;
	//_________________________________________________________

	// Transformation matrices to include support inclinatation
	Eigen::MatrixXd s_transformation_matrix(4, 4);
	s_transformation_matrix.setZero(); // support inclination transformation matrix

	int constraint_type;
	double constraint_angle_rad;
	double support_Lcos;
	double support_Msin;

	// Start node support inclination
	if (model_constarints.constraintMap.find(ln.startNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the start node
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {1.0, 0.0, 0.0, 0.0 } });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 1.0, 0.0, 0.0 } });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.startNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.startNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {support_Lcos, -support_Msin, 0.0, 0.0} });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {support_Msin, support_Lcos, 0.0, 0.0 } });
	}

	// End node support inclination
	if (model_constarints.constraintMap.find(ln.endNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the end node
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ { 0.0, 0.0, 1.0, 0.0 } });
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ { 0.0, 0.0, 0.0, 1.0 } });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.endNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.endNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ {0.0, 0.0, support_Lcos, -support_Msin } });
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ {0.0, 0.0, support_Msin, support_Lcos } });
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

			globalPointMassMatrix((nd_map * 2) + 0, (nd_map * 2) + 0) = ptm.ptmass_x;
			globalPointMassMatrix((nd_map * 2) + 1, (nd_map * 2) + 1) = ptm.ptmass_y;
		}
		else
		{
			// Nodes doesnt have point mass
			globalPointMassMatrix((nd_map * 2) + 0, (nd_map * 2) + 0) = 0.0;
			globalPointMassMatrix((nd_map * 2) + 1, (nd_map * 2) + 1) = 0.0;
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
		Eigen::MatrixXd elementConsistentMassMatrix(4, 4);
		elementConsistentMassMatrix.setZero();

		get_element_consistentmass_matrix(elementConsistentMassMatrix, ln, elementline_material, model_constarints, output_file);

		// Get the Node ID
		int sn_id = nodeid_map[ln.startNode->node_id]; // get the ordered map of the start node ID
		int en_id = nodeid_map[ln.endNode->node_id]; // get the ordered map of the end node ID

		globalConsistentMassMatrix.block<2, 2>(sn_id * 2, sn_id * 2) += elementConsistentMassMatrix.block<2, 2>(0, 0);
		globalConsistentMassMatrix.block<2, 2>(sn_id * 2, en_id * 2) += elementConsistentMassMatrix.block<2, 2>(0, 2);
		globalConsistentMassMatrix.block<2, 2>(en_id * 2, sn_id * 2) += elementConsistentMassMatrix.block<2, 2>(2, 0);
		globalConsistentMassMatrix.block<2, 2>(en_id * 2, en_id * 2) += elementConsistentMassMatrix.block<2, 2>(2, 2);
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
	Eigen::MatrixXd L_transformation_matrix(2, 4);
	L_transformation_matrix.setZero();

	L_transformation_matrix.row(0) = Eigen::RowVectorXd({ { Lcos, Msin, 0.0, 0.0 } });
	L_transformation_matrix.row(1) = Eigen::RowVectorXd({ { 0.0, 0.0, Lcos, Msin } });

	//_________________________________________________________
	// Local element stiffness matrix
	Eigen::MatrixXd local_element_consistentmass_matrix(2, 2);
	local_element_consistentmass_matrix.setZero();

	double k1 = (elementline_material.mat_density * elementline_material.cs_area * eLength) / 2.0f;

	local_element_consistentmass_matrix.row(0) = Eigen::RowVectorXd({ { k1, 0.0 } });
	local_element_consistentmass_matrix.row(1) = Eigen::RowVectorXd({ { 0.0, k1 } });

	//_________________________________________________________
	// Transformed element stiffness matrix
	Eigen::MatrixXd e_consistentMass_matrix(4, 4);
	e_consistentMass_matrix.setZero();

	e_consistentMass_matrix = L_transformation_matrix.transpose() * local_element_consistentmass_matrix * L_transformation_matrix;
	//_________________________________________________________

	// Transformation matrices to include support inclinatation
	Eigen::MatrixXd s_transformation_matrix(4, 4);
	s_transformation_matrix.setZero(); // support inclination transformation matrix

	int constraint_type;
	double constraint_angle_rad;
	double support_Lcos;
	double support_Msin;

	// Start node support inclination
	if (model_constarints.constraintMap.find(ln.startNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the start node
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {1.0, 0.0, 0.0, 0.0 } });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {0.0, 1.0, 0.0, 0.0 } });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.startNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.startNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(0) = Eigen::RowVectorXd({ {support_Lcos, -support_Msin, 0.0, 0.0 } });
		s_transformation_matrix.row(1) = Eigen::RowVectorXd({ {support_Msin, support_Lcos, 0.0, 0.0  } });
	}

	// End node support inclination
	if (model_constarints.constraintMap.find(ln.endNode->node_id) == model_constarints.constraintMap.end())
	{
		// No constraint at the end node
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ { 0.0, 0.0, 1.0, 0.0 } });
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ { 0.0, 0.0, 0.0, 1.0 } });
	}
	else
	{
		constraint_type = model_constarints.constraintMap.at(ln.endNode->node_id).constraint_type; // Constrint type (0 - pin support, 1 - roller support)
		constraint_angle_rad = (model_constarints.constraintMap.at(ln.endNode->node_id).constraint_angle - 90.0) * (m_pi / 180.0f); // Constrint angle in radians
		support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
		support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

		// Pin or Roller Support
		s_transformation_matrix.row(2) = Eigen::RowVectorXd({ { 0.0, 0.0, support_Lcos, -support_Msin } });
		s_transformation_matrix.row(3) = Eigen::RowVectorXd({ { 0.0, 0.0, support_Msin, support_Lcos  } });
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
				// Pin End

				globalDOFMatrix.coeffRef((nd_map * 2) + 0, 0) = 0.0;
				globalDOFMatrix.coeffRef((nd_map * 2) + 1, 0) = 0.0;
			}
			else if (cd.constraint_type == 1)
			{
				// Pin Roller end

				globalDOFMatrix.coeffRef((nd_map * 2) + 0, 0) = 1.0; // X is free to move
				globalDOFMatrix.coeffRef((nd_map * 2) + 1, 0) = 0.0;

				reducedDOF = reducedDOF + 1;
			}
		}
		else
		{
			// Nodes doesnt have Constraint
			globalDOFMatrix.coeffRef((nd_map * 2) + 0, 0) = 1.0;
			globalDOFMatrix.coeffRef((nd_map * 2) + 1, 0) = 1.0;

			reducedDOF = reducedDOF + 2;
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
		if (globalDOFMatrix(((i * 2) + 0), 0) == 0)
		{
			x_free = false;
		}

		// y DOF
		if (globalDOFMatrix(((i * 2) + 1), 0) == 0)
		{
			y_free = false;
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

void modal_analysis_solver::get_reduced_modal_vector_matrix(Eigen::MatrixXd& reduced_eigenVectorsMatrix,
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

void modal_analysis_solver::get_modal_matrices(Eigen::VectorXd& modalMass,
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


void modal_analysis_solver::get_globalSupportInclinationMatrix(Eigen::MatrixXd& globalSupportInclinationMatrix,
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
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 0, (matrix_index * 2) + 0) = 1.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 0, (matrix_index * 2) + 1) = 0.0;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 1, (matrix_index * 2) + 0) = 0.0;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 1, (matrix_index * 2) + 1) = 1.0;
		}
		else
		{
			// Constraint present in this node
			constraint_angle_rad = (model_constarints.constraintMap.at(node_id).constraint_angle - 90.0f) * (m_pi / 180.0f); // Constrint angle in radians
			support_Lcos = std::cos(constraint_angle_rad); // cosine of support inclination
			support_Msin = std::sin(constraint_angle_rad); // sine of support inclination

			// Pin or Roller Support
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 0, (matrix_index * 2) + 0) = support_Lcos;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 0, (matrix_index * 2) + 1) = -1.0 * support_Msin;

			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 1, (matrix_index * 2) + 0) = support_Msin;
			globalSupportInclinationMatrix.coeffRef((matrix_index * 2) + 1, (matrix_index * 2) + 1) = support_Lcos;
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

void modal_analysis_solver::map_modal_analysis_results(const nodes_list_store& model_nodes,
	const elementline_list_store& model_lineelements,
	const nodeconstraint_list_store& model_constarints,
	const modal_analysis_result_store& modal_results,
	const Eigen::MatrixXd& globalSupportInclinationMatrix,
	modal_nodes_list_store& modal_result_nodes,
	modal_elementline_list_store& modal_result_lineelements,
	std::ofstream& output_file)
{
	// Map the modal analysis results to modal members
	// Create a global support inclination matrix
	int numDOF = model_nodes.node_count * 2; // Number of degrees of freedom (2 DOFs per node)
	int node_id = 0;

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

	std::stringstream stopwatch_elapsed_str;
	stopwatch_elapsed_str << std::fixed << std::setprecision(6);

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Transformed eigen vector with support inclination at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	//___________________________________________________________________________________________________________
	// Add to the result nodes
	// Create a matrix to hold the vector lengths of modal displacements
	Eigen::MatrixXd modal_displ_matrix(modal_results.number_of_modes, model_nodes.nodeMap.size());

	for (auto& nd_m : model_nodes.nodeMap)
	{
		node_id = nd_m.first;
		int matrix_index = nodeid_map[node_id];

		// Modal analysis results
		std::unordered_map<int, glm::vec2> node_modal_displ;

		for (int i = 0; i < modal_results.number_of_modes; i++)
		{
			// get the appropriate modal displacement of this particular node
			glm::vec2 modal_displ = glm::vec2(globalEigenVector_transformed((matrix_index * 2) + 0, i),
				globalEigenVector_transformed((matrix_index * 2) + 1, i));

			// Calculate the vector length of the modal displacement
			double vector_length_squared = (modal_displ.x * modal_displ.x) + (modal_displ.y * modal_displ.y);

			// Populate the matrix entry
			modal_displ_matrix(i, matrix_index) = vector_length_squared;

			// add to modal result of this node
			node_modal_displ.insert({ i,modal_displ });
		}

		// Create the modal analysis result node
		glm::vec2 node_pt = model_nodes.nodeMap.at(node_id).node_pt;
		modal_result_nodes.add_result_node(node_id, node_pt, node_modal_displ);
	}

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Results mapped to model nodes at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Add the modal line element result
	for (auto& ln_m : model_lineelements.elementlineMap)
	{
		elementline_store ln = ln_m.second;

		modal_result_lineelements.add_modal_elementline(ln.line_id,
			&modal_result_nodes.modal_nodeMap[ln.startNode->node_id],
			&modal_result_nodes.modal_nodeMap[ln.endNode->node_id]);
	}

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Results mapped to model elements at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Find the maximum displacement for individual modes
	std::unordered_map<int, double> max_node_displ;
	std::unordered_map<int, double> min_node_displ;

	// Find the maximum and minimum of the modal displacement
	Eigen::VectorXd max_values_per_row = modal_displ_matrix.rowwise().maxCoeff();
	Eigen::VectorXd min_values_per_row = modal_displ_matrix.rowwise().minCoeff();

	for (int i = 0; i < modal_results.number_of_modes; i++)
	{
		// Add to the maximum displacement of this mode
		max_node_displ[i] = std::sqrt(max_values_per_row(i));
		min_node_displ[i] = std::sqrt(min_values_per_row(i));
	}

	stopwatch_elapsed_str.str("");
	stopwatch_elapsed_str.clear();
	stopwatch_elapsed_str << stopwatch.elapsed();
	std::cout << "Maximum and minimum modal displacement found at " << stopwatch_elapsed_str.str() << " secs" << std::endl;

	// Set the maximum modal displacement
	modal_result_nodes.max_node_displ.clear();
	modal_result_nodes.max_node_displ = max_node_displ;
	modal_result_nodes.min_node_displ = min_node_displ;

	modal_result_lineelements.max_node_displ.clear();
	modal_result_lineelements.max_node_displ = max_node_displ;
	modal_result_lineelements.min_node_displ = min_node_displ;
}

