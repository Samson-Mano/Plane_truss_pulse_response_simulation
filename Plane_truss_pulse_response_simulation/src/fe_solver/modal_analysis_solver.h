#pragma once
#include <iostream>
#include <fstream>

// FE Objects
#include "../geometry_store/fe_objects/nodes_list_store.h"
#include "../geometry_store/fe_objects/elementline_list_store.h"
#include "../geometry_store/fe_objects/nodeconstraint_list_store.h"
#include "../geometry_store/fe_objects/nodepointmass_list_store.h"

// FE Results Modal Analysis
#include "../geometry_store/modal_result_objects/modal_analysis_result_store.h"
#include "../geometry_store/modal_result_objects/modal_nodes_list_store.h"
#include "../geometry_store/modal_result_objects/modal_elementline_store.h"

// Stop watch
#include "../events_handler/Stopwatch_events.h"

#pragma warning(push)
#pragma warning (disable : 26451)
#pragma warning (disable : 26495)
#pragma warning (disable : 6255)
#pragma warning (disable : 6294)
#pragma warning (disable : 26813)
#pragma warning (disable : 26454)

#include <Eigen/Dense>
#include <Eigen/Sparse>
// Define the sparse matrix type for the reduced global stiffness matrix
typedef Eigen::SparseMatrix<double> SparseMatrix;
#pragma warning(pop)

class modal_analysis_solver
{
public:
	const double m_pi = 3.14159265358979323846;
	bool print_matrix = false;
	Stopwatch_events stopwatch;

	int numDOF = 0;
	int reducedDOF = 0;
	std::unordered_map<int, int> nodeid_map;
	Eigen::MatrixXd globalStiffnessMatrix; // global stiffness matrix
	Eigen::MatrixXd globalMassMatrix; // global mass matrix
	Eigen::MatrixXd globalDOFMatrix; // global DOF matrix
	Eigen::MatrixXd reduced_eigenVectorsMatrix; // Reduced Eigen Vectors matrix
	Eigen::MatrixXd globalSupportInclinationMatrix; // Global support inclination matrix
	Eigen::VectorXd modalMass; // Modal mass matrix
	Eigen::VectorXd modalStiff; // Modal stiffness matrix


	modal_analysis_solver();
	~modal_analysis_solver();
	void modal_analysis_start(const nodes_list_store& model_nodes,
		const elementline_list_store& model_lineelements,
		const nodeconstraint_list_store& model_constarints,
		const nodepointmass_list_store& model_ptmass,
		const std::unordered_map<int, material_data>& material_list,
		const bool& is_include_consistent_mass_matrix,
		modal_analysis_result_store& modal_results,
		modal_nodes_list_store& modal_result_nodes, 
		modal_elementline_list_store& modal_result_lineelements,
		bool& is_modal_analysis_complete);
private:

	void get_global_stiffness_matrix(Eigen::MatrixXd& globalStiffnessMatrix,
		const elementline_list_store& model_lineelements,
		const std::unordered_map<int, material_data>& material_list, 
		const nodeconstraint_list_store& model_constarints,
		std::ofstream& output_file);

	
	void get_element_stiffness_matrix(Eigen::MatrixXd& elementStiffnessMatrix,
		const elementline_store& ln,
		const material_data& elementline_material,
		const nodeconstraint_list_store& model_constarints,
		std::ofstream& output_file);

	void get_global_pointmass_matrix(Eigen::MatrixXd& globalPointMassMatrix,
		const nodes_list_store& model_nodes,
		const nodepointmass_list_store& model_ptmass,
		std::ofstream& output_file);

	void get_global_consistentmass_matrix(Eigen::MatrixXd& globalConsistentMassMatrix,
		const elementline_list_store& model_lineelements,
		const std::unordered_map<int, material_data>& material_list,
		const nodeconstraint_list_store& model_constarints,
		std::ofstream& output_file);

	void get_element_consistentmass_matrix(Eigen::MatrixXd& elementConsistentMassMatrix,
		const elementline_store& ln,
		const material_data& elementline_material,
		const nodeconstraint_list_store& model_constarints,
		std::ofstream& output_file);

	void get_global_dof_matrix(Eigen::MatrixXd& globalDOFMatrix,
		const nodes_list_store& model_nodes,
		const nodeconstraint_list_store& model_constarints,
		int& reducedDOF,
		std::ofstream& output_file);

	void get_reduced_global_matrices(Eigen::MatrixXd& reduced_globalStiffnessMatrix,
		Eigen::MatrixXd& reduced_globalMassMatrix, 
		const Eigen::MatrixXd& globalStiffnessMatrix,
		const Eigen::MatrixXd& globalMassMatrix,
		const Eigen::MatrixXd& globalDOFMatrix,
		const int& numDOF,
		std::ofstream& output_file);

	void get_global_modal_vector_matrix(Eigen::MatrixXd& eigenvectors,
		const Eigen::MatrixXd& eigenvectors_reduced,
		const Eigen::MatrixXd& globalDOFMatrix,
		const int& numDOF,
		int& reducedDOF,
		std::ofstream& output_file);

	void get_reduced_modal_vector_matrix(Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const modal_analysis_result_store& modal_results,
		int& reducedDOF,
		std::ofstream& output_file);

	void get_modal_matrices(Eigen::VectorXd& modalMass,
		Eigen::VectorXd& modalStiff,
		const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const Eigen::MatrixXd& reduced_globalMassMatrix,
		const Eigen::MatrixXd& reduced_globalStiffnessMatrix,
		const int& reducedDOF,
		std::ofstream& output_file);

	void sort_eigen_values_vectors(Eigen::VectorXd& eigenvalues,
		Eigen::MatrixXd& eigenvectors,
		const int& m_size);

	void normalize_eigen_vectors(Eigen::MatrixXd& eigenvectors,
		const int& m_size);

	Eigen::MatrixXd convert_vector_to_1Dmatrix(const std::vector<double>& vec);

	int get_number_of_rigid_body_modes(int num_of_nodes,
		const Eigen::MatrixXd& globalDOFMatrix);

	void get_globalSupportInclinationMatrix(Eigen::MatrixXd& globalSupportInclinationMatrix,
		const nodes_list_store& model_nodes,
		const nodeconstraint_list_store& model_constarints,
		const int& numDOF,
		std::ofstream& output_file);

	void map_modal_analysis_results(const nodes_list_store& model_nodes,
		const elementline_list_store& model_lineelements,
		const nodeconstraint_list_store& model_constarints,
		const modal_analysis_result_store& modal_results,
		const Eigen::MatrixXd& globalSupportInclinationMatrix,
		modal_nodes_list_store& modal_result_nodes,
		modal_elementline_list_store& modal_result_lineelements,
		std::ofstream& output_file);

};


