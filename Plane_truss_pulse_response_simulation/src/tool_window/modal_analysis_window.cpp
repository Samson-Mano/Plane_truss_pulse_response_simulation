#include "modal_analysis_window.h"

modal_analysis_window::modal_analysis_window()
{
	// Empty Constructor
}

modal_analysis_window::~modal_analysis_window()
{
	// Empty Destructor
}

void modal_analysis_window::init()
{
	is_show_window = false;
	execute_modal_analysis = false;
	execute_open = false;
	execute_close = false;
	selected_modal_option = 0;
	mode_result_str.clear();
}

void modal_analysis_window::render_window()
{
	if (is_show_window == false)
		return;

	ImGui::Begin("Modal Analysis Solver");

	// Add check boxes to Include consistent mass matrix
	ImGui::Checkbox("Include Consistent Mass", &is_include_consistent_mass_matrix);
	ImGui::Spacing();

	// Add a Modal Analysis button
	if (ImGui::Button("Modal Analysis"))
	{
		execute_modal_analysis = true;
	}

	// Add the modal analysis result list box
	std::vector<const char*> items_cstr;
	for (const auto& item:mode_result_str) 
	{
		items_cstr.push_back(item.c_str());
	}

	ImGui::ListBox("Natural Frequency", &selected_modal_option, items_cstr.data(), static_cast<int>( items_cstr.size()));

	// Capture if selection changed or not
	if (selected_modal_option != selection_change_flag)
	{
		selection_change_flag = selected_modal_option;
		is_mode_selection_changed = true;
	}

	ImGui::Spacing();
	//_________________________________________________________________________________________


	// Add check boxes to show the Deformed model
	ImGui::Checkbox("Show Model", &show_undeformed_model);
	ImGui::Checkbox("Show Result values", &show_result_text_values);


	ImGui::Spacing();
	//_________________________________________________________________________________________

	// Input box to give input via text
	static bool defscale_input_mode = false;
	static char defscale_str[16] = ""; // buffer to store input deformation scale string
	static double defscale_input = 0; // buffer to store input deformation scale value

	// Button to switch to input mode
	if (!defscale_input_mode)
	{
		if (ImGui::Button("Deformation Scale"))
		{
			defscale_input_mode = true;
			snprintf(defscale_str, 16, "%.1f", deformation_scale_max); // set the buffer to current deformation scale value
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Deformation Scale", defscale_str, IM_ARRAYSIZE(defscale_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			defscale_input = atoi(defscale_str);
			// set the load value to input value
			deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			defscale_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", deformation_scale_max);

	// Slider for Deflection scale
	float deformation_scale_flt = static_cast<float>(deformation_scale_max);

	ImGui::Text("Deformation Scale");
	ImGui::SameLine();
	ImGui::SliderFloat(".", &deformation_scale_flt, 0.0f, 100.0f, "%.1f");
	deformation_scale_max = deformation_scale_flt;

	//Set the deformation scale
	normailzed_defomation_scale = 1.0f;
	deformation_scale = deformation_scale_max;

	ImGui::Spacing();
	//_________________________________________________________________________________________

	if (ImGui::CollapsingHeader("Animate"))
	{
		// Animate the solution
		// Start a horizontal layout
		ImGui::BeginGroup();

		// Play button active
		if (animate_play == true)
		{
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
		}
		else
		{
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
		}

		// Add the Play button
		if (ImGui::Button("Play"))
		{
			// Handle Play button click
			animate_play = !animate_play;
			animate_pause = false;
		}
		ImGui::PopStyleColor(1);  // Revert back to default style

		// Add some spacing between buttons
		ImGui::SameLine();


		// Pause button active
		if (animate_pause == true)
		{
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
		}
		else
		{
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
		}

		// Add the Pause button
		if (ImGui::Button("Pause"))
		{
			// Handle Pause button click
			animate_pause = !animate_pause;
			animate_play = false;
		}
		ImGui::PopStyleColor(1);  // Revert back to default style

		// Add some spacing between buttons
		ImGui::SameLine();

		// Add the Stop button
		if (ImGui::Button("Stop"))
		{
			// Handle Stop button click
			animate_play = false;
			animate_pause = false;
			time_val = 0.0f;
		}

		// Animation speed control
		// Input box to give input via text
		static bool animation_speed_input_mode = false;
		static char animation_speed_str[16] = ""; // buffer to store input deformation scale string
		static float animation_speed_input = 0; // buffer to store input deformation scale value

		// Button to switch to input mode
		if (!animation_speed_input_mode)
		{
			if (ImGui::Button("Animation Speed"))
			{
				animation_speed_input_mode = true;
				snprintf(animation_speed_str, 16, "%.1f", animation_speed); // set the buffer to current deformation scale value
			}
		}
		else // input mode
		{
			// Text box to input value
			ImGui::SetNextItemWidth(60.0f);
			if (ImGui::InputText("##Animation Speed", animation_speed_str, IM_ARRAYSIZE(animation_speed_str), ImGuiInputTextFlags_CharsDecimal))
			{
				// convert the input string to int
				animation_speed_input = static_cast<float>(atof(animation_speed_str));
				// set the load value to input value
				animation_speed = animation_speed_input;
			}

			// Button to switch back to slider mode
			ImGui::SameLine();
			if (ImGui::Button("OK"))
			{
				animation_speed_input_mode = false;
			}
		}

		// Text for Animation speed value
		ImGui::SameLine();
		ImGui::Text(" %.1f", animation_speed);

		// Display the frame rate
		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
			1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

		// End the horizontal layout
		ImGui::EndGroup();
	}

	ImGui::Spacing();
	//_________________________________________________________________________________________

		// Close button
	if (ImGui::Button("Close"))
	{
		execute_close = true;
		execute_open = false;
		is_show_window = false; // set the flag to close the window
	}

	ImGui::End();


	// Set the animation data
	if (animate_play == true)
	{
		if (time_val > 6.283185307f)
		{
			time_val = 0.0f;
		}

		// Animation is playing 
		normailzed_defomation_scale = std::sin(time_val * animation_speed); // Varies between 0 and 1
		deformation_scale = normailzed_defomation_scale * deformation_scale_max;
		time_val = time_val + 0.0002f;
	}
	else if (animate_pause == true)
	{
		normailzed_defomation_scale = std::sin(time_val * animation_speed) ; // Varies between 0 and 1
		deformation_scale = normailzed_defomation_scale * deformation_scale_max;
	}

}
