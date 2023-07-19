#pragma once
#include <iostream>
#include <vector>
#include "../ImGui/imgui.h"
#include "../ImGui/imgui_impl_glfw.h"
#include "../ImGui/imgui_impl_opengl3.h"
#include "../ImGui/stb_implement.h"
#include <chrono>


class Stopwatch
{
public:
	void reset_time();
	double current_elapsed() const;

private:
	std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime = std::chrono::high_resolution_clock::time_point();
	// std::chrono::time_point<std::chrono::high_resolution_clock> m_endTime;
};

class pulse_response_window
{
public:
	bool is_show_window = false;
	bool execute_pulse_analysis = false; // Main solver run event flag
	bool execute_open = false; // Solver window open event flag
	bool execute_close = false; // Closing of solution window event flag

	// Inputs for response calculation
	double total_simulation_time = 10.0;
	double time_interval = 0.1;
	double damping_ratio = 0.01;
	
	// Modal analysis Results
	double modal_first_frequency = 0.0;
	double modal_end_frequency = 0.0;
	int number_of_modes = 0;

	// Pulse response analysis results
	bool pulse_response_analysis_complete = false;
	
	bool show_undeformed_model = true; // show undeformed model 

	// Animation control
	bool animate_play = true;
	bool animate_pause = false;
	double deformation_scale_max = 10.0;
	double animation_speed = 1.0;
	// double normailzed_defomation_scale = 0.0;
	// double deformation_scale = 0.0;

	// Time step control
	double time_interval_atrun = 0.0; // Value of time interval used in the pulse response 
	int time_step_count = 0;
	int time_step = 0;

	pulse_response_window();
	~pulse_response_window();
	void init();
	void render_window();
private:
	Stopwatch stopwatch;
};

