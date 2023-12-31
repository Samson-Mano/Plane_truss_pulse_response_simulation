#include "pulse_nodes_list_store.h"

pulse_nodes_list_store::pulse_nodes_list_store()
{
	// Empty constructor
}

pulse_nodes_list_store::~pulse_nodes_list_store()
{
	// Empty destructor
}

void pulse_nodes_list_store::init(geom_parameters* geom_param_ptr)
{
	// Initialize
	// Set the geometry parameters
	this->geom_param_ptr = geom_param_ptr;

	// Set the geometry parameters for the points
	pulse_node_points.init(geom_param_ptr);

	// Clear the results
	node_count = 0;
	pulse_nodeMap.clear();
	pulse_node_points.clear_points();
	max_node_displ = 0.0; // Maximum nodal displacement
}

void pulse_nodes_list_store::clear_data()
{
	// Clear the results
	node_count = 0;
	pulse_nodeMap.clear();
	pulse_node_points.clear_points();
	max_node_displ = 0.0; // Maximum nodal displacement
}

void pulse_nodes_list_store::add_result_node(int& node_id, glm::vec2& node_pt, pulse_node_result node_pulse_result, const int& number_of_time_steps)
{
	// Add the result nodes
	pulse_node_store temp_pulse_result;

	temp_pulse_result.node_id = node_id; // Add the node ID
	temp_pulse_result.node_pt = node_pt; // Add the node point
	temp_pulse_result.node_pulse_result = node_pulse_result; // Node point results
	temp_pulse_result.number_of_timesteps = number_of_time_steps; // Number of time steps

	// Check whether the node_id is already there
	if (pulse_nodeMap.find(node_id) != pulse_nodeMap.end())
	{
		// Node ID already exist (do not add)
		return;
	}

	// Add to the pulse nodeMap
	pulse_nodeMap.insert({ node_id,temp_pulse_result });
	node_count++;
}

void pulse_nodes_list_store::set_buffer()
{
	// Clear the points
	pulse_node_points.clear_points();

	//__________________________ Add the Dynamic points
	for (auto& nd_m : pulse_nodeMap)
	{
		pulse_node_store nd = nd_m.second; // get the node data

		std::vector<glm::vec2> point_offset; // point offset
		std::vector<glm::vec3> point_color; // point color

		for(int i = 0; i <static_cast<int>(nd.node_pulse_result.node_pulse_displ.size());i++)
		{
			// Create the node offset
			// Point displacement
			double pt_displ = std::sqrt(std::pow(nd.node_pulse_result.node_pulse_displ[i].x, 2) +
				std::pow(nd.node_pulse_result.node_pulse_displ[i].y, 2));

			// Distance ratio
			double dist_ratio = pt_displ / max_node_displ;

			// Create the point offset
			glm::vec2 temp_point_offset = glm::vec2(nd.node_pulse_result.node_pulse_displ[i].x / max_node_displ,
				nd.node_pulse_result.node_pulse_displ[i].y / max_node_displ);

			point_offset.push_back(temp_point_offset);

			// Create the node color vectors
			glm::vec3 temp_point_color = getContourColor(static_cast<float>(1.0 - dist_ratio));

			point_color.push_back(temp_point_color);
		}

		// Add to the pulse points
		pulse_node_points.add_point(nd.node_id, nd.node_pt, point_offset, point_color);
	}

	// Set buffer
	pulse_node_points.set_buffer();
}

void pulse_nodes_list_store::paint_pulse_nodes(const int& dyn_index)
{
	// Paint the points
	pulse_node_points.paint_points(dyn_index);
}

void pulse_nodes_list_store::update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale)
{
	// Pulse node points update geometry 
	pulse_node_points.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
}

glm::vec3 pulse_nodes_list_store::getContourColor(float value)
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