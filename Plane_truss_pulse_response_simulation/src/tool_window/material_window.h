#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include "../ImGui/imgui.h"
#include "../ImGui/imgui_impl_glfw.h"
#include "../ImGui/imgui_impl_opengl3.h"

struct material_data
{
	unsigned int material_id = 0;
	std::string material_name = "";
	double youngs_mod = 0.0;
	double mat_density = 0.0;
	double cs_area = 0.0;
};

class material_window
{
public:
	bool is_show_window = false;
	bool is_assign_material = false;
	int execute_delete_materialid = -1;

	int selected_material_option = 0;
	std::unordered_map<int, material_data> material_list;

	material_window();
	~material_window();
	void init();
	void render_window();
	static glm::vec3 get_standard_color(int color_index);
private:
	int selected_list_option = 0;
	int get_unique_material_id();
};
