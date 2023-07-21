#pragma once
#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/vec2.hpp>
#include "../ImGui/imgui.h"
#include "../ImGui/imgui_impl_glfw.h"
#include "../ImGui/imgui_impl_opengl3.h"
#include "../ImGui/stb_implement.h"

class inlcondition_window
{
public:
	bool is_show_window = false;
	bool is_add_initial_condition = false;
	// Initial displacement
	double initial_displacement_x = 0.0;
	double initial_displacement_y = 0.0;
	// Initial velocity
	double initial_velocity_x = 0.0;
	double initial_velocity_y = 0.0;
	
	inlcondition_window();
	~inlcondition_window();
	void init();
	void render_window();
private:

};
