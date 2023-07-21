#include "nodeinlcond_list_store.h"

nodeinlcond_list_store::nodeinlcond_list_store()
{
	// Empty constructor
}

nodeinlcond_list_store::~nodeinlcond_list_store()
{
	// Empty destructor
}

void nodeinlcond_list_store::init(geom_parameters* geom_param_ptr)
{
	// Set the geometry parameters
	this->geom_param_ptr = geom_param_ptr;
	inl_condition_displ_labels.init(geom_param_ptr);
	inl_condition_velo_labels.init(geom_param_ptr);

	// Clear the intial condition
	inlcond_count = 0;
	inlcondMap.clear();
}

void nodeinlcond_list_store::add_inlcondition(int& node_id, glm::vec2& inlcond_loc,
	double& inl_disp_x, double& inl_disp_y, double& inl_velo_x, double& inl_velo_y)
{
	// Add the initial condition to the particular node
	nodeinl_condition_data temp_inl_condition_data;
	temp_inl_condition_data.node_id = node_id;
	temp_inl_condition_data.inlcond_loc = inlcond_loc;
	temp_inl_condition_data.inl_displacement_x = inl_disp_x;
	temp_inl_condition_data.inl_displacement_y = inl_disp_y;
	temp_inl_condition_data.inl_velocity_x = inl_velo_x;
	temp_inl_condition_data.inl_velocity_y = inl_velo_y;

	// Insert the inital condition data to unordered map
	// Searching for node_id
	if (inlcondMap.find(node_id) != inlcondMap.end())
	{
		// Node is already have constraint
		// so remove the constraint
		inlcondMap[node_id] = temp_inl_condition_data;

		return;
	}

	// Insert the constraint to nodes
	inlcondMap.insert({ node_id, temp_inl_condition_data });
	inlcond_count++;
}

void nodeinlcond_list_store::delete_inlcondition(int& node_id)
{
	// Delete the initial condition in this node
	if (inlcond_count != 0)
	{
		// Remove the intial condition data to unordered map
		// Searching for node_id
		// Check there is already a initial conditon in the found node
		if (inlcondMap.find(node_id) != inlcondMap.end())
		{
			// Node is already have initial condition
			// so remove the intial condition
			inlcondMap.erase(node_id);

			// Update the buffer
			set_buffer();

			// adjust the initial condition count
			inlcond_count--;
		}
	}
}

void nodeinlcond_list_store::set_buffer()
{
	// Set the buffer for initial condition label
	if (inlcond_count == 0)
	{
		// No initial condition label to paint
		return;
	}

	// Clear the labels
	inl_condition_displ_labels.clear_labels();
	inl_condition_velo_labels.clear_labels();

	// Fill the label data
	for (auto& inl_data_m : inlcondMap)
	{
		nodeinl_condition_data inl_data = inl_data_m.second;

		glm::vec3 temp_color = geom_param_ptr->geom_colors.inlcond_displ_color;
		std::string	temp_str = "";
		std::stringstream ss_displ_x;
		std::stringstream ss_displ_y;
		
		// Displacement initial condition
		if (inl_data.inl_displacement_x != 0 || inl_data.inl_displacement_y != 0)
		{
			// initial displacement x
			ss_displ_x << std::fixed << std::setprecision(geom_param_ptr->inlcond_precision) << inl_data.inl_displacement_x;
			// initial displacement y
			ss_displ_y << std::fixed << std::setprecision(geom_param_ptr->inlcond_precision) << inl_data.inl_displacement_y;

			temp_str = "(" + ss_displ_x.str() + ", " + ss_displ_y.str() + ")";

			// Add to the diplacement initial condition label
			inl_condition_displ_labels.add_text(temp_str, inl_data.inlcond_loc, glm::vec2(0), temp_color, 0.0, true, false);
		}


		temp_color = geom_param_ptr->geom_colors.inlcond_velo_color;
		temp_str = "";
		std::stringstream ss_velo_x;
		std::stringstream ss_velo_y;

		// Velocity initial condition
		if (inl_data.inl_velocity_x != 0 || inl_data.inl_velocity_y != 0)
		{
			// initial velocity x
			ss_velo_x << std::fixed << std::setprecision(geom_param_ptr->inlcond_precision) << inl_data.inl_velocity_x;
			// initial velocity y
			ss_velo_y << std::fixed << std::setprecision(geom_param_ptr->inlcond_precision) << inl_data.inl_velocity_y;

			temp_str = "(" + ss_velo_x.str() + ", " + ss_velo_y.str() + ")";

			// Add to the velocity initial condition label (Paint below)
			inl_condition_velo_labels.add_text(temp_str, inl_data.inlcond_loc, glm::vec2(0), temp_color, 0.0, false, false);
		}
	}

	// Set the buffer
	inl_condition_displ_labels.set_buffer();
	inl_condition_velo_labels.set_buffer();
}

void nodeinlcond_list_store::paint_inlcondition_label()
{
	// Paint the labels of intial condition
	inl_condition_displ_labels.paint_text();
	inl_condition_velo_labels.paint_text();
}

void nodeinlcond_list_store::update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale)
{
	// Update the geometry matrices of initial condition label
	inl_condition_displ_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
	inl_condition_velo_labels.update_opengl_uniforms(set_modelmatrix, set_pantranslation, set_zoomtranslation, set_transparency, set_deflscale);
}
