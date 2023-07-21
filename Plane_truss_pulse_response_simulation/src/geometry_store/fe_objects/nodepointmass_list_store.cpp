#include "nodepointmass_list_store.h"

nodepointmass_list_store::nodepointmass_list_store()
{
	// Empty constructor
}

nodepointmass_list_store::~nodepointmass_list_store()
{
	// Empty destructor
}

void nodepointmass_list_store::init(geom_parameters* geom_param_ptr)
{
	// Set the geometry parameters
	this->geom_param_ptr = geom_param_ptr;
	ptmass_value_labels.init(geom_param_ptr);

	// Create the shader and Texture for the drawing the constraints
	std::filesystem::path shadersPath = geom_param_ptr->resourcePath;

	ptmass_shader.create_shader((shadersPath.string() + "/resources/shaders/ptmass_vert_shader.vert").c_str(),
		(shadersPath.string() + "/resources/shaders/ptmass_frag_shader.frag").c_str());

	// Set texture uniform
	ptmass_shader.setUniform("u_Texture", 0);

	// Load the texture
	ptmass_texture.LoadTexture((shadersPath.string() + "/resources/images/pic_3D_circle.png").c_str());

	// Clear the Point mass
	ptmass_count = 0;
	max_ptmass_value = 0.0;
	ptmassMap.clear();
}

void nodepointmass_list_store::add_pointmass(int& node_id, glm::vec2& ptmass_loc, glm::vec2 ptmass_defl, 
	double& ptmass_x, double& ptmass_y, bool is_offset)
{
	// Add the point mass
	nodepointmass_data temp_ptmass_data;
	temp_ptmass_data.node_id = node_id;
	temp_ptmass_data.ptmass_loc = ptmass_loc;
	temp_ptmass_data.ptmass_defl = ptmass_defl;
	temp_ptmass_data.ptmass_x = ptmass_x;
	temp_ptmass_data.ptmass_y = ptmass_y;
	temp_ptmass_data.is_offset = is_offset;

	// Insert the constarint data to unordered map
	// Searching for node_id
	if (ptmassMap.find(node_id) != ptmassMap.end())
	{
		// Node is already have constraint
		// so remove the constraint
		ptmassMap[node_id] = temp_ptmass_data;

		return;
	}

	// Insert the constraint to nodes
	ptmassMap.insert({ node_id, temp_ptmass_data });
	ptmass_count++;
}

void nodepointmass_list_store::delete_pointmass(int& node_id)
{
	// Delete the point mass
	if (ptmass_count != 0)
	{
		// Remove the point mass data to unordered map
		// Searching for node_id
		// Check there is already a point mass in the found node
		if (ptmassMap.find(node_id) != ptmassMap.end())
		{
			// Node is already have point mass
			// so remove the point mass
			ptmassMap.erase(node_id);

			// Update the buffer
			set_buffer();

			// adjust the point mass count
			ptmass_count--;
		}
	}
}

void nodepointmass_list_store::set_buffer()
{
	// Set the buffer for point mass
	if (ptmass_count == 0)
	{
		// No point mass to paint
		return;
	}

	// Find the maximum mass values
	max_ptmass_value = 0.0;
	// Set the load lables
	ptmass_value_labels.clear_labels();

	// Find the load maximum
	for (auto& ptmx : ptmassMap)
	{
		nodepointmass_data ptm = ptmx.second;

		if (max_ptmass_value < ptm.ptmass_x)
		{
			// mass x
			max_ptmass_value = ptm.ptmass_x;
		}

		if (max_ptmass_value < ptm.ptmass_y)
		{
			// mass y
			max_ptmass_value = ptm.ptmass_y;
		}

		//__________________________________________________________________________

		glm::vec3 temp_color = geom_param_ptr->geom_colors.ptmass_color;
		std::string	temp_str = "";

		if (ptm.ptmass_x != 0)
		{
			std::stringstream ss;
			ss << std::fixed << std::setprecision(geom_param_ptr->load_precision) << ptm.ptmass_x;

			temp_str = temp_str + "m_x = " + ss.str();
		}
		//__________________________________________________________________________
		if (ptm.ptmass_y != 0)
		{
			if (temp_str != "")
			{
				// Add a comma if value is already there
				temp_str = temp_str + ", ";
			}

			std::stringstream ss;
			ss << std::fixed << std::setprecision(geom_param_ptr->load_precision) << ptm.ptmass_y;

			temp_str = temp_str + "m_y = " + ss.str();
		}

		ptmass_value_labels.add_text(temp_str, ptm.ptmass_loc, glm::vec2(0), temp_color, 0.0, true, false);
	}

	ptmass_value_labels.set_buffer();

	//__________________________________________________________________________

	unsigned int ptmass_vertex_count = 4 * 12 * ptmass_count;
	float* ptmass_vertices = new float[ptmass_vertex_count];

	unsigned int ptmass_indices_count = 6 * ptmass_count;
	unsigned int* ptmass_indices = new unsigned int[ptmass_indices_count];

	unsigned int ptmass_v_index = 0;
	unsigned int ptmass_i_index = 0;

	for (auto& ptmx : ptmassMap)
	{
		nodepointmass_data ptm = ptmx.second;

		// Add the texture buffer
		get_constraint_buffer(ptm, ptmass_vertices, ptmass_v_index, ptmass_indices, ptmass_i_index);
	}

	VertexBufferLayout ptmass_layout;
	ptmass_layout.AddFloat(2);  // Position
	ptmass_layout.AddFloat(2);  // Center
	ptmass_layout.AddFloat(2);  // Defl
	ptmass_layout.AddFloat(3);  // Color
	ptmass_layout.AddFloat(2);  // Texture co-ordinate
	ptmass_layout.AddFloat(1);  // Is Deflection

	unsigned int ptmass_vertex_size = ptmass_vertex_count * sizeof(float);

	// Create the Constraint buffers
	ptmass_buffer.CreateBuffers(ptmass_vertices, ptmass_vertex_size,
		ptmass_indices, ptmass_indices_count, ptmass_layout);

	// Delete the Dynamic arrays
	delete[] ptmass_vertices;
	delete[] ptmass_indices;
}

void nodepointmass_list_store::paint_pointmass()
{
	// Paint the point mass
	ptmass_texture.Bind();
	ptmass_shader.Bind();
	ptmass_buffer.Bind();
	glDrawElements(GL_TRIANGLES, 6 * ptmass_count, GL_UNSIGNED_INT, 0);
	ptmass_buffer.UnBind();
	ptmass_shader.UnBind();
	ptmass_texture.UnBind();

}

void nodepointmass_list_store::paint_pointmass_label()
{
	// Paint the label of point mass
	ptmass_value_labels.paint_text();
}

void nodepointmass_list_store::update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale)
{
	// Update the geometry matrices of point mass labels
	ptmass_value_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);

	if (set_modelmatrix == true)
	{
		// set the model matrix
		ptmass_shader.setUniform("geom_scale", static_cast<float>(geom_param_ptr->geom_scale));
		ptmass_shader.setUniform("transparency", 1.0f);

		ptmass_shader.setUniform("modelMatrix", geom_param_ptr->modelMatrix, false);
	}

	if (set_pantranslation == true)
	{
		// set the pan translation
		ptmass_shader.setUniform("panTranslation", geom_param_ptr->panTranslation, false);
	}

	if (set_zoomtranslation == true)
	{
		// set the zoom translation
		ptmass_shader.setUniform("zoomscale", static_cast<float>(geom_param_ptr->zoom_scale));
	}

	if (set_transparency == true)
	{
		// set the alpha transparency
		ptmass_shader.setUniform("transparency", static_cast<float>(geom_param_ptr->geom_transparency));
	}

	if (set_deflscale == true)
	{
		// set the deflection scale
		// ptmass_shader.setUniform("deflscale", static_cast<float>(geom_param_ptr->defl_scale));
	}
}

void nodepointmass_list_store::get_constraint_buffer(nodepointmass_data& ptm, float* ptmass_vertices, unsigned int& ptmass_v_index,
	unsigned int* ptmass_indices, unsigned int& ptmass_i_index)
{
	// Constraint color
	glm::vec3 ptmass_color = geom_param_ptr->geom_colors.ptmass_color;
	double mass_val = std::max(ptm.ptmass_x, ptm.ptmass_y);
	float corner_size = static_cast<float>(-3.4 * (mass_val / max_ptmass_value) * (geom_param_ptr->node_circle_radii / geom_param_ptr->geom_scale));

	// Set the Point mass vertices Corner 1 Top Left
	ptmass_vertices[ptmass_v_index + 0] = ptm.ptmass_loc.x - corner_size;
	ptmass_vertices[ptmass_v_index + 1] = ptm.ptmass_loc.y + corner_size;

	// Set the node center
	ptmass_vertices[ptmass_v_index + 2] = ptm.ptmass_loc.x;
	ptmass_vertices[ptmass_v_index + 3] = ptm.ptmass_loc.y;

	// Defl value
	ptmass_vertices[ptmass_v_index + 4] = ptm.ptmass_defl.x;
	ptmass_vertices[ptmass_v_index + 5] = ptm.ptmass_defl.y;

	// Set the node color
	ptmass_vertices[ptmass_v_index + 6] = ptmass_color.x;
	ptmass_vertices[ptmass_v_index + 7] = ptmass_color.y;
	ptmass_vertices[ptmass_v_index + 8] = ptmass_color.z;

	// Set the Texture co-ordinates
	ptmass_vertices[ptmass_v_index + 9] = 0.0;
	ptmass_vertices[ptmass_v_index + 10] = 0.0;

	// Set the defl_value
	ptmass_vertices[ptmass_v_index + 11] = static_cast<float>(ptm.is_offset);

	// Increment
	ptmass_v_index = ptmass_v_index + 12;


	// Set the Point mass vertices Corner 2 Top Right
	ptmass_vertices[ptmass_v_index + 0] = ptm.ptmass_loc.x + corner_size;
	ptmass_vertices[ptmass_v_index + 1] = ptm.ptmass_loc.y + corner_size;

	// Set the node center
	ptmass_vertices[ptmass_v_index + 2] = ptm.ptmass_loc.x;
	ptmass_vertices[ptmass_v_index + 3] = ptm.ptmass_loc.y;

	// Defl value
	ptmass_vertices[ptmass_v_index + 4] = ptm.ptmass_defl.x;
	ptmass_vertices[ptmass_v_index + 5] = ptm.ptmass_defl.y;

	// Set the node color
	ptmass_vertices[ptmass_v_index + 6] = ptmass_color.x;
	ptmass_vertices[ptmass_v_index + 7] = ptmass_color.y;
	ptmass_vertices[ptmass_v_index + 8] = ptmass_color.z;

	// Set the Texture co-ordinates
	ptmass_vertices[ptmass_v_index + 9] = 1.0;
	ptmass_vertices[ptmass_v_index + 10] = 0.0;

	// Set the defl_value
	ptmass_vertices[ptmass_v_index + 11] = static_cast<float>(ptm.is_offset);

	// Increment
	ptmass_v_index = ptmass_v_index + 12;


	// Set the Point Mass vertices Corner 3 Bot Right
	ptmass_vertices[ptmass_v_index + 0] = ptm.ptmass_loc.x + corner_size;
	ptmass_vertices[ptmass_v_index + 1] = ptm.ptmass_loc.y - corner_size;

	// Set the node center
	ptmass_vertices[ptmass_v_index + 2] = ptm.ptmass_loc.x;
	ptmass_vertices[ptmass_v_index + 3] = ptm.ptmass_loc.y;

	// Defl value
	ptmass_vertices[ptmass_v_index + 4] = ptm.ptmass_defl.x;
	ptmass_vertices[ptmass_v_index + 5] = ptm.ptmass_defl.y;

	// Set the node color
	ptmass_vertices[ptmass_v_index + 6] = ptmass_color.x;
	ptmass_vertices[ptmass_v_index + 7] = ptmass_color.y;
	ptmass_vertices[ptmass_v_index + 8] = ptmass_color.z;

	// Set the Texture co-ordinates
	ptmass_vertices[ptmass_v_index + 9] = 1.0;
	ptmass_vertices[ptmass_v_index + 10] = 1.0;

	// Set the defl_value
	ptmass_vertices[ptmass_v_index + 11] = static_cast<float>(ptm.is_offset);

	// Increment
	ptmass_v_index = ptmass_v_index + 12;


	// Set the Constraint vertices Corner 4 Bot Left
	ptmass_vertices[ptmass_v_index + 0] = ptm.ptmass_loc.x - corner_size;
	ptmass_vertices[ptmass_v_index + 1] = ptm.ptmass_loc.y - corner_size;

	// Set the node center
	ptmass_vertices[ptmass_v_index + 2] = ptm.ptmass_loc.x;
	ptmass_vertices[ptmass_v_index + 3] = ptm.ptmass_loc.y;

	// Defl value
	ptmass_vertices[ptmass_v_index + 4] = ptm.ptmass_defl.x;
	ptmass_vertices[ptmass_v_index + 5] = ptm.ptmass_defl.y;

	// Set the node color
	ptmass_vertices[ptmass_v_index + 6] = ptmass_color.x;
	ptmass_vertices[ptmass_v_index + 7] = ptmass_color.y;
	ptmass_vertices[ptmass_v_index + 8] = ptmass_color.z;

	// Set the Texture co-ordinates
	ptmass_vertices[ptmass_v_index + 9] = 0.0;
	ptmass_vertices[ptmass_v_index + 10] = 1.0;

	// Set the defl_value
	ptmass_vertices[ptmass_v_index + 11] = static_cast<float>(ptm.is_offset);

	// Increment
	ptmass_v_index = ptmass_v_index + 12;


	//______________________________________________________________________

	// Set the Quad indices
	unsigned int t_id = ((ptmass_i_index / 6) * 4);

	// Triangle 0,1,2
	ptmass_indices[ptmass_i_index + 0] = t_id + 0;
	ptmass_indices[ptmass_i_index + 1] = t_id + 1;
	ptmass_indices[ptmass_i_index + 2] = t_id + 2;

	// Triangle 2,3,0
	ptmass_indices[ptmass_i_index + 3] = t_id + 2;
	ptmass_indices[ptmass_i_index + 4] = t_id + 3;
	ptmass_indices[ptmass_i_index + 5] = t_id + 0;

	// Increment
	ptmass_i_index = ptmass_i_index + 6;


}
