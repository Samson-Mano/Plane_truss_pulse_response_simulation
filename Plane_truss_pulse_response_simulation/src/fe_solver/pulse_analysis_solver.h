#pragma once
#include <iostream>
#include <fstream>

// FE Objects
#include "../geometry_store/fe_objects/nodes_list_store.h"
#include "../geometry_store/fe_objects/elementline_list_store.h"
#include "../geometry_store/fe_objects/nodeconstraint_list_store.h"
#include "../geometry_store/fe_objects/nodeload_list_store.h"
#include "../geometry_store/fe_objects/nodepointmass_list_store.h"

// FE Results Modal Analysis
#include "../geometry_store/modal_result_objects/modal_analysis_result_store.h"

// FE Results Freq Analysis
#include "../geometry_store/pulse_result_objects/pulse_analysis_result_store.h"
#include "../geometry_store/pulse_result_objects/pulse_nodes_list_store.h"
#include "../geometry_store/pulse_result_objects/pulse_elementline_store.h"

#pragma warning(push)
#pragma warning (disable : 26451)
#pragma warning (disable : 26495)
#pragma warning (disable : 6255)
#pragma warning (disable : 26813)
#pragma warning (disable : 26454)

#include <Eigen/Dense>
#include <Eigen/Sparse>
// Define the sparse matrix type for the reduced global stiffness matrix
typedef Eigen::SparseMatrix<double> SparseMatrix;
#pragma warning(pop)

struct pulse_load_data
{
	int load_id = 0;
	double load_start_time = 0.0;
	double load_end_time = 0.0;
	Eigen::MatrixXd modal_reducedLoadamplMatrix;
	Eigen::MatrixXd reducedLoadamplMatrix;
	Eigen::MatrixXd globalLoadamplMatrix;
};

class pulse_analysis_solver
{
public:	
	const double m_pi = 3.14159265358979323846;
	const double epsilon = 0.000001;
	bool print_matrix = true;

	pulse_analysis_solver();
	~pulse_analysis_solver();
	void pulse_analysis_start(const nodes_list_store& model_nodes,
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
		bool& is_pulse_analysis_complete);
private:
	std::unordered_map<int, int> nodeid_map;

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

	void get_modal_matrices(Eigen::VectorXd& modalMass,
		Eigen::VectorXd& modalStiff,
		const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const Eigen::MatrixXd& reduced_globalMassMatrix,
		const Eigen::MatrixXd& reduced_globalStiffnessMatrix,
		const int& reducedDOF,
		std::ofstream& output_file);

	void get_reduced_modal_vector_matrix(Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const modal_analysis_result_store& modal_results,
		int& reducedDOF,
		std::ofstream& output_file);

	void get_globalSupportInclinationMatrix(Eigen::MatrixXd& globalSupportInclinationMatrix,
		const nodes_list_store& model_nodes,
		const nodeconstraint_list_store& model_constarints,
		const int& numDOF,
		std::ofstream& output_file);

	void create_pulse_load_matrices(pulse_load_data& pulse_loads,
		const load_data& ld,
		const elementline_list_store& model_lineelements,
		const Eigen::MatrixXd& globalDOFMatrix,
		const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const int& numDOF, 
		const int& reducedDOF);

	void get_element_load_matrix(Eigen::MatrixXd& elementLoadMatrix,
		const elementline_store& ln,
		const load_data& ld);

	void get_global_resp_matrix(Eigen::MatrixXd& displ_ampl_RespMatrix_b4supp_trans,
		const Eigen::MatrixXd& displ_ampl_RespMatrix_reduced,
		const Eigen::MatrixXd& globalDOFMatrix,
		const int& numDOF,
		const int& reducedDOF);

	void get_steady_state_pulse_soln(double& steady_state_displ_resp,
		const double& time_t,
		const double& modal_mass,
		const double& modal_stiff,
		const double& modal_force_ampl,
		const double& modal_force_starttime,
		const double& modal_force_endtime);

	void map_pulse_analysis_results(pulse_analysis_result_store& pulse_response_result,
		pulse_nodes_list_store& pulse_result_nodes,
		pulse_elementline_list_store& pulse_result_lineelements,
		const int& number_of_time_steps,
		const nodes_list_store& model_nodes,
		const elementline_list_store& model_lineelements,
		const std::unordered_map<int, pulse_node_result>& node_results);
};
