#pragma once
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <unordered_map>
#include "../geom_parameters.h"
#include "../geometry_buffers/gBuffers.h"
#include "../geometry_objects/label_list_store.h"

struct nodepointmass_data
{
	int node_id = 0;
	glm::vec2 ptmass_loc = glm::vec2(0);
	glm::vec2 ptmass_defl = glm::vec2(0);
	double ptmass_x = 0.0;
	double ptmass_y = 0.0;
	double ptmass_xy = 0.0;
	bool is_offset = 0.0;
};


class nodepointmass_list_store
{
public:
	const double ptmass_quad_size = 0.04;
	unsigned int ptmass_count = 0;
	std::unordered_map<int, nodepointmass_data> ptmassMap;

	nodepointmass_list_store();
	~nodepointmass_list_store();
	void init(geom_parameters* geom_param_ptr);
	void add_pointmass(int& node_id, glm::vec2& ptmass_loc, glm::vec2 ptmass_defl, double& ptmass_x, double& ptmass_y, double& ptmass_xy, bool is_offset);
	void delete_pointmass(int& node_id);
	void set_buffer();
	void paint_pointmass();
	void paint_pointmass_label();
	void update_geometry_matrices(bool set_modelmatrix, bool set_pantranslation, bool set_zoomtranslation, bool set_transparency, bool set_deflscale);

private:
	geom_parameters* geom_param_ptr = nullptr;
	gBuffers ptmass_buffer;
	Shader ptmass_shader;
	Texture ptmass_texture;
	label_list_store ptmass_value_labels;
	double max_ptmass_value = 0.0;

	void get_constraint_buffer(nodepointmass_data& ptm, float* ptmass_vertices, unsigned int& ptmass_v_index,
		unsigned int* ptmass_indices, unsigned int& ptmass_i_index);
};
