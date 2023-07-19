#include "elementline_list_store.h"

elementline_list_store::elementline_list_store()
{
	// Empty constructor
}

elementline_list_store::~elementline_list_store()
{
	// Empty destructor
}

void elementline_list_store::init(geom_parameters* geom_param_ptr)
{
	// Set the geometry parameters
	this->geom_param_ptr = geom_param_ptr;

	// Set the geometry parameters for the labels (and clear the labels)
	element_lines.init(geom_param_ptr);
	line_id_labels.init(geom_param_ptr);
	line_length_labels.init(geom_param_ptr);
	line_material_id_labels.init(geom_param_ptr);

	// Clear the lines
	elementline_count = 0;
	elementlineMap.clear();
}

void elementline_list_store::add_elementline(int& line_id, node_store* startNode, node_store* endNode,int& material_id)
{
	// Add the line to the list
	elementline_store temp_line;
	temp_line.line_id = line_id;
	temp_line.startNode = startNode;
	temp_line.endNode = endNode;
	temp_line.material_id = material_id;

	// Check whether the node_id is already there
	if (elementlineMap.find(line_id) != elementlineMap.end())
	{
		// Element ID already exist (do not add)
		return;
	}

	// Insert to the lines
	elementlineMap.insert({ line_id, temp_line });
	elementline_count++;

	//__________________________ Add the node points
	glm::vec3 temp_color = geom_param_ptr->geom_colors.line_color;
	glm::vec2 start_node_pt = (*startNode).node_pt;
	glm::vec2 end_node_pt = (*endNode).node_pt;

	//__________________________ Add the lines
	element_lines.add_line(line_id, start_node_pt, end_node_pt,
		glm::vec2(0), glm::vec2(0), temp_color, temp_color, false);

	//__________________________ Add the node labels
	std::string temp_str = "[" + std::to_string(line_id) + "]";
	glm::vec2 line_mid_pt = (start_node_pt + end_node_pt)/2.0f;

	// Calculate the angle between the line segment and the x-axis
	float line_angle = atan2(end_node_pt.y - start_node_pt.y, end_node_pt.x - start_node_pt.x);

	line_id_labels.add_text(temp_str, line_mid_pt, glm::vec2(0), temp_color, line_angle, true, false);

	// Add the node coordinate label
	float line_length = glm::distance(start_node_pt, end_node_pt);

	std::stringstream ss;
	ss << std::fixed << std::setprecision(geom_param_ptr->length_precision) << line_length;

	temp_str = ss.str();

	line_length_labels.add_text(temp_str, line_mid_pt, glm::vec2(0), temp_color, line_angle, false, false);
}

void elementline_list_store::set_buffer()
{
	// Set the buffers for the Model
	element_lines.set_buffer();
	line_id_labels.set_buffer();
	line_length_labels.set_buffer();
	update_material_id_labels();
}

void elementline_list_store::paint_elementlines()
{
	// Paint the model lines
	element_lines.paint_lines();
}

void elementline_list_store::paint_label_line_ids()
{
	// Paint the line id labels
	line_id_labels.paint_text();
}

void elementline_list_store::paint_label_line_lengths()
{
	// Paint the line length labels
	line_length_labels.paint_text();
}

void elementline_list_store::paint_lines_material_id()
{
	// Paint line material ID
	line_material_id_labels.paint_text();
}

int elementline_list_store::is_line_hit(glm::vec2& loc)
{
	// Return the line id of line which is clicked
	// Covert mouse location to screen location
	int max_dim = geom_param_ptr->window_width > geom_param_ptr->window_height ? geom_param_ptr->window_width : geom_param_ptr->window_height;

	// Transform the mouse location to openGL screen coordinates
	glm::vec2 screenPt = glm::vec2(2.0f * ((loc.x - (geom_param_ptr->window_width * 0.5f)) / max_dim),
		2.0f * (((geom_param_ptr->window_height * 0.5f) - loc.y) / max_dim));


	// Nodal location
	glm::mat4 scaling_matrix = glm::mat4(1.0) * static_cast<float>(geom_param_ptr->zoom_scale);
	scaling_matrix[3][3] = 1.0f;

	glm::mat4 scaledModelMatrix = scaling_matrix * geom_param_ptr->modelMatrix;

	// Loop through all nodes in map and update min and max values
	for (auto it = elementlineMap.begin(); it != elementlineMap.end(); ++it)
	{
		elementline_store elementline = it->second;

		glm::vec2 s_node = elementline.startNode->node_pt;
	 	glm::vec2 e_node = elementline.endNode->node_pt;

		glm::vec4 s_node_finalPosition = scaledModelMatrix * glm::vec4(s_node.x, s_node.y, 0, 1.0f) * geom_param_ptr->panTranslation;
		glm::vec4 e_node_finalPosition = scaledModelMatrix * glm::vec4(e_node.x, e_node.y, 0, 1.0f) * geom_param_ptr->panTranslation;

		// S & E Point 
		glm::vec2 spt = glm::vec2(s_node_finalPosition.x, s_node_finalPosition.y);
		glm::vec2 ept = glm::vec2(e_node_finalPosition.x, e_node_finalPosition.y);

		float threshold = 8 * geom_param_ptr->node_circle_radii;

		if (isClickPointOnLine(screenPt, spt, ept, threshold) == true)
		{
			// Return the Id of the line if hit == true
			return it->first;
		}
	}

	return -1;
}

bool elementline_list_store::isClickPointOnLine(const glm::vec2& clickPoint, const glm::vec2& lineStart, const glm::vec2& lineEnd, float threshold)
{
	glm::vec2 lineDirection = lineEnd - lineStart;
	float lineLengthSq = glm::dot(lineDirection, lineDirection);

	glm::vec2 clickToLineStart = clickPoint - lineStart;
	float dotProduct = glm::dot(clickToLineStart, lineDirection);

	// Calculate the normalized projection of clickToLineStart onto the line
	glm::vec2 projection = (dotProduct / lineLengthSq) * lineDirection;

	// Calculate the squared normal distance between the click point and the line
	float normalDistanceSq = glm::dot(clickToLineStart - projection, clickToLineStart - projection);

	// Check if the click point is within the line segment's bounding box
	if (dotProduct >= 0.0f && dotProduct <= lineLengthSq)
	{
		// Check if the normal distance is less than or equal to the threshold
		if (normalDistanceSq <= threshold * threshold)
		{
			return true; // Click point is on the line segment
		}
	}

	return false; // Click point is not on the line segment
}

void elementline_list_store::update_material_id_labels()
{
	line_material_id_labels.clear_labels();

	// Update the material ids
	glm::vec3 temp_color;
	std::string temp_str;

	for (auto it = elementlineMap.begin(); it != elementlineMap.end(); ++it)
	{
		elementline_store elementline = it->second;

		// Add the line id
		glm::vec2 start_pt = elementline.startNode->node_pt;
		glm::vec2 end_pt = elementline.endNode->node_pt;

		// Calculate the midpoint of the line segment
		glm::vec2 line_mid_pt = glm::vec2((start_pt.x + end_pt.x) * 0.5f, (start_pt.y + end_pt.y) * 0.5f);

		float line_angle = atan2(end_pt.y - start_pt.y, end_pt.x - start_pt.x);

		// Add the material ID
		temp_color = material_window::get_standard_color(elementline.material_id);
		temp_str = "                 M = " + std::to_string(elementline.material_id);
		line_material_id_labels.add_text(temp_str, line_mid_pt, glm::vec2(0), temp_color, line_angle, true, false);
	}

	line_material_id_labels.set_buffer();
}

void elementline_list_store::update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale)
{
	// Update model openGL uniforms
	element_lines.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
	line_id_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
	line_length_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
	line_material_id_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
}
