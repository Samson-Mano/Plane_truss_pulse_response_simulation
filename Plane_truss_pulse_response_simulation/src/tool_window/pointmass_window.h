#pragma once
#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/vec2.hpp>
#include "../ImGui/imgui.h"
#include "../ImGui/imgui_impl_glfw.h"
#include "../ImGui/imgui_impl_opengl3.h"
#include "../ImGui/stb_implement.h"

class pointmass_window
{
public:
	bool is_show_window = false;
	bool is_add_pointmass = false;
	double mass_x = 100;
	double mass_y = 100;
	double mass_xy = 100;


	pointmass_window();
	~pointmass_window();
	void init();
	void render_window();
private:

};
