#pragma once
#include <iostream>
#include "../ImGui/stb_image.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/vec2.hpp>
#include "../ImGui/imgui.h"
#include "../ImGui/imgui_impl_glfw.h"
#include "../ImGui/imgui_impl_opengl3.h"
#include "../ImGui/stb_implement.h"

struct cnst_image_data
{
	int image_width = 0;
	int image_height = 0;
	unsigned int image_texture_ID = 0;
	bool is_loaded = false;
};

class constraint_window
{
public:
	bool is_show_window = false;
	bool is_add_constraint = false;

	// Constraint type 0 - fixed support, 1 - fixed roller, 2 - pin support, 3 - pin roller
	int constraint_type = 0; // Constraint type
	double constraint_angle = 90.0; // Constraint angle

	constraint_window();
	~constraint_window();
	void init(); // initialize bind images
	void render_window();
private:
	cnst_image_data cnst_image;

	void draw_support();
	void LoadTextureFromFile(const char* filename, cnst_image_data& cnst_image);
	bool get_image_min_max_coord(ImVec2& window_pos, ImVec2& window_size, ImVec2& img_pos_top_left,
		ImVec2& img_pos_top_right, ImVec2& img_pos_bot_right, ImVec2& img_pos_bot_left, double& orientation);
};
