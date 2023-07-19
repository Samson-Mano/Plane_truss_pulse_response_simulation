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

struct load_image_data
{
	int image_width = 0;
	int image_height = 0;
	unsigned int image_texture_ID = 0;
	bool is_loaded = false;
};


class load_window
{
public:
	bool is_show_window = false;
	bool is_add_load = false;

	double load_amplitude = 100.0; // load value
	double load_param = 0.5; // Load parameter 
	double load_start_time = 0.2;
	double load_end_time = 0.6;
	double load_angle = 90.0; // load angle


	load_window();
	~load_window();
	void init(); // initialize bind images
	void render_window();
private:
	load_image_data load_image;

	void draw_load();
	void LoadTextureFromFile(const char* filename, load_image_data& load_image);
	bool get_image_min_max_coord(ImVec2& window_pos, ImVec2& window_size, ImVec2& img_pos_top_left,
		ImVec2& img_pos_top_right, ImVec2& img_pos_bot_right, ImVec2& img_pos_bot_left, double& orientation);
};
