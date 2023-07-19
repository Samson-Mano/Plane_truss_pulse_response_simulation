#pragma once
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <unordered_map>
#include "../geom_parameters.h"
#include "../geometry_buffers/gBuffers.h"
#include "../geometry_objects/label_list_store.h"
#include "../geometry_objects/point_list_store.h"

struct modal_node_store
{
	int node_id = 0;
	glm::vec2 node_pt = glm::vec2(0);

	// Modal results (mode number, (x, y, rot_xy))
	std::unordered_map<int, glm::vec3> node_modal_displ;
};

class modal_nodes_list_store
{
public:
	unsigned int node_count = 0;
	std::unordered_map<int, modal_node_store> modal_nodeMap; // Create an unordered_map to store nodes with ID as key
	std::unordered_map<int, double> max_node_displ; // Stores the maximum nodal displacement for the whole model
	std::unordered_map<int, double> min_node_displ; // Stores the minimum nodal displacement for the whole model

	modal_nodes_list_store();
	~modal_nodes_list_store();
	void init(geom_parameters* geom_param_ptr);
	void clear_data();
	void add_result_node(int& node_id, glm::vec2& node_pt, std::unordered_map<int, glm::vec3> node_modal_displ);
	void set_buffer(int selected_mode);
	void paint_modal_nodes();
	void paint_label_mode_vectors();
	void update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale);

private:
	geom_parameters* geom_param_ptr = nullptr;
	label_list_store modal_node_vector_labels;
	point_list_store modal_node_points;

	glm::vec3 getContourColor(float value);
};
