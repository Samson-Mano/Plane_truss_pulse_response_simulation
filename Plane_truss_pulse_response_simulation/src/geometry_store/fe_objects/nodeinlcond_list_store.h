#pragma once
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <unordered_map>
#include "../geom_parameters.h"
#include "../geometry_buffers/gBuffers.h"
#include "../geometry_objects/label_list_store.h"

struct nodeinl_condition_data
{
	int node_id = 0;
	glm::vec2 inlcond_loc = glm::vec2(0);
	double inl_displacement_x = 0.0; // initial displacement x
	double inl_displacement_y = 0.0; // initial displacement y
	double inl_velocity_x = 0.0; // initial velocity x
	double inl_velocity_y = 0.0; // initial velocity y
};

class nodeinlcond_list_store
{
public:
	unsigned int inlcond_count = 0;
	 std::unordered_map<int, nodeinl_condition_data> inlcondMap;

	nodeinlcond_list_store();
	~nodeinlcond_list_store();
	void init(geom_parameters* geom_param_ptr);
	void add_inlcondition(int& node_id, glm::vec2& inlcond_loc, double& inl_disp_x, double& inl_disp_y, double& inl_velo_x, double& inl_velo_y);
	void delete_inlcondition(int& node_id);
	void set_buffer();
	void paint_inlcondition_label();
	void update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale);

private:
	geom_parameters* geom_param_ptr = nullptr;
	label_list_store inl_condition_displ_labels;
	label_list_store inl_condition_velo_labels;
};
