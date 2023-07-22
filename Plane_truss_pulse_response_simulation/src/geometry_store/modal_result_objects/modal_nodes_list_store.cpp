#include "modal_nodes_list_store.h"

modal_nodes_list_store::modal_nodes_list_store()
{
	// Empty constructor
}

modal_nodes_list_store::~modal_nodes_list_store()
{
	// Empty destructor
}

void modal_nodes_list_store::init(geom_parameters* geom_param_ptr)
{
	// Set the geometry parameters
	this->geom_param_ptr = geom_param_ptr;

	// Set the geometry parameters for the labels (and clear the labels)
	modal_node_points.init(geom_param_ptr);
	modal_node_vector_labels.init(geom_param_ptr);

	// Clear the results
	node_count = 0;
	modal_nodeMap.clear();
	max_node_displ.clear();
}

void modal_nodes_list_store::clear_data()
{
	// Clear the results
	node_count = 0;
	modal_nodeMap.clear();
	max_node_displ.clear();
}

void modal_nodes_list_store::add_result_node(int& node_id, glm::vec2& node_pt, std::unordered_map<int, glm::vec2> node_modal_displ)
{
	// Add result nodes
	modal_node_store temp_node;
	temp_node.node_id = node_id;
	temp_node.node_pt = node_pt;
	temp_node.node_modal_displ = node_modal_displ;

	// Check whether the node_id is already there
	if (modal_nodeMap.find(node_id) != modal_nodeMap.end())
	{
		// Node ID already exist (do not add)
		return;
	}

	// Insert to the nodes
	modal_nodeMap.insert({ node_id, temp_node });
	node_count++;
}

void modal_nodes_list_store::set_buffer(int selected_mode)
{
	// Clear existing modal line
	modal_node_points.clear_points();
	modal_node_vector_labels.clear_labels();

	double mode_max_displ = max_node_displ[selected_mode];
	double mode_min_displ = min_node_displ[selected_mode];

	// Add the lines
	// Loop throug every line element
	int i = 0;
	for (auto& nd_m : modal_nodeMap)
	{
		modal_node_store nd = nd_m.second;

		glm::vec2 pt_displ = glm::vec2(nd.node_modal_displ[selected_mode].x / mode_max_displ, nd.node_modal_displ[selected_mode].y / mode_max_displ);

		// Find the displacment value
		double pt_displ_value = std::sqrt(std::pow(nd.node_modal_displ[selected_mode].x, 2) +
			std::pow(nd.node_modal_displ[selected_mode].y, 2));

		// scale the value with mode max displacement
		double pt_displ_scale = (pt_displ_value - mode_min_displ) / (mode_max_displ - mode_min_displ);

		glm::vec3 pt_contour_color = getContourColor(static_cast<float>(1.0 - pt_displ_scale));

		// Add all the points
		modal_node_points.add_point(nd.node_id, nd.node_pt, pt_displ, pt_contour_color, true);


		//__________________________ Add the node labels
		// Mode x
		std::stringstream mode_x;
		mode_x << std::fixed << std::setprecision(3) << nd.node_modal_displ[selected_mode].x;

		// Mode y
		std::stringstream mode_y;
		mode_y << std::fixed << std::setprecision(3) << nd.node_modal_displ[selected_mode].y;


		std::string temp_str = "(" + mode_x.str() + ", " + mode_y.str() + " )";
		
		modal_node_vector_labels.add_text(temp_str, nd.node_pt, pt_displ, pt_contour_color, 0.0f, true, true);
	}

	// Set the buffer
	modal_node_points.set_buffer();
	modal_node_vector_labels.set_buffer();
}

void modal_nodes_list_store::paint_modal_nodes()
{
	modal_node_points.paint_points();
}

void modal_nodes_list_store::paint_label_mode_vectors()
{
	modal_node_vector_labels.paint_text();
}

void modal_nodes_list_store::update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale)
{
	modal_node_points.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
	modal_node_vector_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
}


glm::vec3 modal_nodes_list_store::getContourColor(float value)
{
	// return the contour color based on the value (0 to 1)
	glm::vec3 color;
	float r, g, b;

	// Rainbow color map
	float hue = value * 5.0f; // Scale the value to the range of 0 to 5
	float c = 1.0f;
	float x = c * (1.0f - glm::abs(glm::mod(hue / 2.0f, 1.0f) - 1.0f));
	float m = 0.0f;

	if (hue >= 0 && hue < 1) {
		r = c;
		g = x;
		b = m;
	}
	else if (hue >= 1 && hue < 2) {
		r = x;
		g = c;
		b = m;
	}
	else if (hue >= 2 && hue < 3) {
		r = m;
		g = c;
		b = x;
	}
	else if (hue >= 3 && hue < 4) {
		r = m;
		g = x;
		b = c;
	}
	else {
		r = x;
		g = m;
		b = c;
	}

	color = glm::vec3(r, g, b);
	return color;


}