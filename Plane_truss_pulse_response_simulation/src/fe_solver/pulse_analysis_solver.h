#pragma once
#include <iostream>
#include <fstream>

// Modal analysis solver
#include "../fe_solver/modal_analysis_solver.h"

// FE Objects
#include "../geometry_store/fe_objects/nodeload_list_store.h"
#include "../geometry_store/fe_objects/nodeinlcond_list_store.h"

// FE Results Freq Analysis
#include "../geometry_store/pulse_result_objects/pulse_analysis_result_store.h"
#include "../geometry_store/pulse_result_objects/pulse_nodes_list_store.h"
#include "../geometry_store/pulse_result_objects/pulse_elementline_store.h"

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
	bool print_matrix = false;

	pulse_analysis_solver();
	~pulse_analysis_solver();
	void pulse_analysis_start(const nodes_list_store& model_nodes,
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
		bool& is_pulse_analysis_complete);
private:
	std::unordered_map<int, int> nodeid_map;

	void create_initial_condition_matrices(Eigen::MatrixXd& modal_reducedInitialDisplacementMatrix,
		Eigen::MatrixXd& modal_reducedInitialVelocityMatrix,
		const nodeinlcond_list_store& model_inlcond,
		const nodes_list_store& model_nodes,
		const Eigen::MatrixXd& globalDOFMatrix,
		const Eigen::MatrixXd& reduced_eigenVectorsMatrix,
		const int& numDOF,
		const int& reducedDOF);

	void create_pulse_load_matrices(pulse_load_data& pulse_loads,
		const load_data& ld,
		const nodes_list_store& model_nodes,
		const Eigen::MatrixXd& globalDOFMatrix,
		const Eigen::MatrixXd& reduced_eigenVectorsMatrix_transpose,
		const int& numDOF, 
		const int& reducedDOF);

	void get_reduced_global_matrix(Eigen::MatrixXd& reducedglobalMatrix,
		const Eigen::MatrixXd& globalMatrix,
		const Eigen::MatrixXd& globalDOFMatrix,
		const int& numDOF,
		const int& reducedDOF);

	void get_global_resp_matrix(Eigen::MatrixXd& displ_ampl_RespMatrix_b4supp_trans,
		const Eigen::MatrixXd& displ_ampl_RespMatrix_reduced,
		const Eigen::MatrixXd& globalDOFMatrix,
		const int& numDOF,
		const int& reducedDOF);

	void get_steady_state_initial_condition_soln(double& steady_state_displ_resp,
		const double& time_t,
		const double& modal_mass,
		const double& modal_stiff,
		const double& modal_initial_displacement,
		const double& modal_initial_velocity);

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
