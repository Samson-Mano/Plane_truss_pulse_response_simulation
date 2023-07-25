#pragma once
#include "pulse_nodes_list_store.h"
#include "../geometry_objects/dynamic_line_list_store.h"

struct pulse_line_points
{
	int split_line_id = 0; // line id of the individual hermitian interpolation line 
	// Point coordinate
	glm::vec2 pt1 = glm::vec2(0);
	glm::vec2 pt2 = glm::vec2(0);

	// Point displacements
	std::vector<glm::vec2> pt1_modal_displ;
	std::vector<glm::vec2> pt2_modal_displ;
};

struct pulse_elementline_store
{
	int line_id = 0; // ID of the line
	pulse_node_store* startNode = nullptr; // start node
	pulse_node_store* endNode = nullptr; // end node

	// Line pulse displacement data
	std::vector<pulse_line_points> discretized_bar_line_data;
};

class pulse_elementline_list_store
{
public:
	const int interpolation_count = 3;
	unsigned int pulse_elementline_count = 0;
	std::unordered_map<int, pulse_elementline_store> pulse_elementlineMap; // Create an unordered_map to store lines with ID as key
	double max_line_displ = 0.0; // Maximum line displacement

	pulse_elementline_list_store();
	~pulse_elementline_list_store();
	void init(geom_parameters* geom_param_ptr);
	void clear_data();
	void add_pulse_elementline(int& line_id, pulse_node_store* startNode, pulse_node_store* endNode);
	std::vector<pulse_line_points> set_line_bar_interpolation(const int& interpolation_count, pulse_node_store* startNode, pulse_node_store* endNode);
	double linear_bar_element_interpolation(double q1, double q2, double s);
	double hermite_beam_element_interpolation(double v1, double theta1, double v2, double theta2, double s);

	void set_buffer();
	void paint_pulse_elementlines(const int& dyn_index);
	void update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale);
private:
	geom_parameters* geom_param_ptr = nullptr;
	dynamic_line_list_store pulse_element_lines;

	glm::vec3 getContourColor(float value);
};
