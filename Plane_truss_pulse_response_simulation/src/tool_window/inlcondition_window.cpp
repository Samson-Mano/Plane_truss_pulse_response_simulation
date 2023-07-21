#include "inlcondition_window.h"

inlcondition_window::inlcondition_window()
{
	// Empty constructor
}

inlcondition_window::~inlcondition_window()
{
	// Empty destructor
}

void inlcondition_window::init()
{
}

void inlcondition_window::render_window()
{
	if (is_show_window == false)
		return;

	ImGui::Begin("Initial Condition");

	ImGui::Text("Initial Displacement");
	// ________________________________________________________________________________________________________________________________
	// Input box to give input via text
	static bool displx_input_mode = false;
	static char displx_str[16] = ""; // buffer to store input Displacement X string
	static float displx_input = 0; // buffer to store input Displacement X value

	// Button to switch to input mode
	if (!displx_input_mode)
	{
		if (ImGui::Button("X Displacmenet"))
		{
			displx_input_mode = true;
			snprintf(displx_str, 16, "%.1f", initial_displacement_x); // set the buffer to current Displacement X
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Displacement X", displx_str, IM_ARRAYSIZE(displx_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			initial_displacement_x = atof(displx_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			displx_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", initial_displacement_x);

	// ________________________________________________________________________________________________________________________________
	// Input box to give input via text
	static bool disply_input_mode = false;
	static char disply_str[16] = ""; // buffer to store input Displacement Y string
	static float disply_input = 0; // buffer to store input Displacement Y value

	// Button to switch to input mode
	if (!disply_input_mode)
	{
		if (ImGui::Button("Y Displacement"))
		{
			disply_input_mode = true;
			snprintf(disply_str, 16, "%.1f", initial_displacement_y); // set the buffer to current Displacement Y
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Displacement Y", disply_str, IM_ARRAYSIZE(disply_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			initial_displacement_y = atof(disply_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			disply_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", initial_displacement_y);

	ImGui::Spacing();
	ImGui::Spacing();

	ImGui::Text("Initial Velocity");

	// ________________________________________________________________________________________________________________________________
	// Input box to give input via text
	static bool velox_input_mode = false;
	static char velox_str[16] = ""; // buffer to store input Velocity X string
	static float velox_input = 0; // buffer to store input Velocity X value

	// Button to switch to input mode
	if (!velox_input_mode)
	{
		if (ImGui::Button("X Velocity"))
		{
			velox_input_mode = true;
			snprintf(velox_str, 16, "%.1f", initial_velocity_x); // set the buffer to current Velocity X
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Velocity X", velox_str, IM_ARRAYSIZE(velox_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			initial_velocity_x = atof(velox_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			velox_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", initial_velocity_x);

	// ________________________________________________________________________________________________________________________________
	// Input box to give input via text
	static bool veloy_input_mode = false;
	static char veloy_str[16] = ""; // buffer to store input Velocity Y string
	static float veloy_input = 0; // buffer to store input Velocity Y value

	// Button to switch to input mode
	if (!veloy_input_mode)
	{
		if (ImGui::Button("Y Velocity"))
		{
			veloy_input_mode = true;
			snprintf(veloy_str, 16, "%.1f", initial_velocity_y); // set the buffer to current Velocity Y
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Velocity Y", veloy_str, IM_ARRAYSIZE(veloy_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			initial_velocity_y = atof(veloy_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			veloy_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", initial_velocity_y);

	ImGui::Spacing();
	ImGui::Spacing();

	// ________________________________________________________________________________________________________________________________
	// Add Initial Condition
	if (is_add_initial_condition == true)
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
	}
	else
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
	}

	if (ImGui::Button("Add Intial Condition"))
	{
		is_add_initial_condition = !is_add_initial_condition;
	}

	// Text to inform user whether Add initial condition is added or not
	ImGui::SameLine(); // move cursor to the same line
	if (is_add_initial_condition == true)
	{
		ImGui::Text("Click on the node to add the initial condition");
	}
	else
	{
		ImGui::Text("<= Click here to start adding point mass");
	}

	ImGui::PopStyleColor(1); // Pop the style color after using it

	ImGui::Spacing();
	ImGui::Spacing();
	//_________________________________________________________________________________________

	// Close button
	if (ImGui::Button("Close"))
	{
		is_show_window = false; // set the flag to close the window
		is_add_initial_condition = false; // clear the add initial condition flag
	}

	ImGui::End();
}
