#pragma once
#include "geom_parameters.h"
// File system
#include <fstream>
#include <sstream>
#include <iomanip>

// Window includes
#include "../tool_window/constraint_window.h"
#include "../tool_window/load_window.h"
#include "../tool_window/material_window.h"
#include "../tool_window/options_window.h"
#include "../tool_window/pointmass_window.h"
#include "../tool_window/solver_window.h"
#include "../tool_window/modal_analysis_window.h"
#include "../tool_window/pulse_response_window.h"

// Solver
#include "../fe_solver/modal_analysis_solver.h"
#include "../fe_solver/pulse_analysis_solver.h"

// FE Objects
#include "fe_objects/nodes_list_store.h"
#include "fe_objects/elementline_list_store.h"
#include "fe_objects/nodeconstraint_list_store.h"
#include "fe_objects/nodeload_list_store.h"
#include "fe_objects/nodepointmass_list_store.h"

// FE Results Modal Analysis
#include "modal_result_objects/modal_analysis_result_store.h"
#include "modal_result_objects/modal_nodes_list_store.h"
#include "modal_result_objects/modal_elementline_store.h"

// FE Results Result Objects Pulse response
#include "pulse_result_objects/pulse_analysis_result_store.h"
#include "pulse_result_objects/pulse_elementline_store.h"
#include "pulse_result_objects/pulse_nodes_list_store.h"

class geom_store
{
public:
	const double m_pi = 3.14159265358979323846;
	bool is_geometry_set = false;

	// Main Variable to strore the geometry parameters
	geom_parameters geom_param;

	geom_store();
	~geom_store();
	void init(options_window* op_window, material_window* mat_window, 
		modal_analysis_window* sol_modal_window, pulse_response_window* sol_pulse_window);
	void fini();

	// Reading and writing the geometry file
	void read_varai2d(std::ifstream& input_file);
	void read_dxfdata(std::ostringstream& input_data);
	void read_rawdata(std::ifstream& input_file);
	void write_rawdata(std::ofstream& output_file);

	// Functions to control the origin
	void update_WindowDimension(const int& window_width, const int& window_height);
	void update_model_matrix();
	void update_model_zoomfit();
	void update_model_pan(glm::vec2& transl);
	void update_model_zoom(double& z_scale);
	void update_model_transperency(bool is_transparent);

	// Function to add/ remove loads, constraints, lumped mass and material properties
	void set_nodal_constraint(glm::vec2 mouse_click_loc, int& constraint_type, double& constraint_angle, bool is_add);
	void set_member_load(glm::vec2 mouse_click_loc, double& load_start_time, double& load_end_time,
		double& load_value, double& load_angle, bool is_add);
	void set_elementline_material(glm::vec2 mouse_click_loc);
	void set_nodal_pointmass(glm::vec2 mouse_click_loc, double& pt_mass_x, double& pt_mass_y, 
		double& pt_mass_xy, bool is_add);

	// Functions to paint the geometry and results
	void paint_geometry();
private:
	// main variables to store the geometry (Nodes, lines, loads, supports)
	// Geometry objects
	nodes_list_store model_nodes;
	elementline_list_store model_lineelements;
	nodeconstraint_list_store model_constarints;
	nodeload_list_store model_loads;
	nodepointmass_list_store model_ptmass;

	// Modal Analysis results
	modal_analysis_result_store modal_results;
	modal_nodes_list_store modal_result_nodes;
	modal_elementline_list_store modal_result_lineelements;

	// Pulse analysis results
	pulse_analysis_result_store pulse_response_result;
	pulse_nodes_list_store pulse_result_nodes;
	pulse_elementline_list_store pulse_result_lineelements;

	// View options ptr and Material window ptr
	options_window* op_window = nullptr;
	material_window* mat_window = nullptr;
	modal_analysis_window* sol_modal_window = nullptr;
	pulse_response_window* sol_pulse_window = nullptr;

	// Analysis 
	bool is_modal_analysis_complete = false;
	bool is_freq_analysis_complete = false;
	bool is_pulse_analysis_complete = false;

	// Create geometry private functions
	void create_geometry(nodes_list_store& model_nodes,
		elementline_list_store& model_lineelements,
		nodeconstraint_list_store& model_constarints,
		nodeload_list_store& model_loads, 
		nodepointmass_list_store& model_ptmass);

	std::pair<glm::vec2, glm::vec2> findMinMaxXY(const std::unordered_map<int, node_store>& model_nodes);
	glm::vec2 findGeometricCenter(const std::unordered_map<int, node_store>& model_nodes);
	void update_delete_material(int& del_material_id);

	void paint_model(); // Paint the model
	void paint_modal_analysis(); // Paint the modal analysis results
	void paint_pulse_analysis(); // Paint the pulse response analysis results
};
