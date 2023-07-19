#include "pointmass_window.h"

pointmass_window::pointmass_window()
{
	// Empty constructor
}

pointmass_window::~pointmass_window()
{
	// Empty destructor
}

void pointmass_window::init()
{

}

void pointmass_window::render_window()
{
	if (is_show_window == false)
		return;

	ImGui::Begin("Point Mass");

	// ________________________________________________________________________________________________________________________________
	// Input box to give input via text
	static bool massx_input_mode = false;
	static char massx_str[16] = ""; // buffer to store input Mass X string
	static float massx_input = 0; // buffer to store input Mass X value

	// Button to switch to input mode
	if (!massx_input_mode)
	{
		if (ImGui::Button("Mass X"))
		{
			massx_input_mode = true;
			snprintf(massx_str, 16, "%.1f", mass_x); // set the buffer to current Mass X
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Mass X", massx_str, IM_ARRAYSIZE(massx_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			mass_x = atof(massx_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			massx_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", mass_x);


	// ________________________________________________________________________________________________________________________________
		// Input box to give input via text
	static bool massy_input_mode = false;
	static char massy_str[16] = ""; // buffer to store input Mass Y string
	static float massy_input = 0; // buffer to store input Mass Y value

	// Button to switch to input mode
	if (!massy_input_mode)
	{
		if (ImGui::Button("Mass Y"))
		{
			massy_input_mode = true;
			snprintf(massy_str, 16, "%.1f", mass_y); // set the buffer to current Mass Y
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Mass Y", massy_str, IM_ARRAYSIZE(massy_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			mass_y = atof(massy_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			massy_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", mass_y);

	// ________________________________________________________________________________________________________________________________
		// Input box to give input via text
	static bool massxy_input_mode = false;
	static char massxy_str[16] = ""; // buffer to store input mass xy string
	static float massxy_input = 0; // buffer to store input mass xy value

	// Button to switch to input mode
	if (!massxy_input_mode)
	{
		if (ImGui::Button("Mass XY"))
		{
			massxy_input_mode = true;
			snprintf(massxy_str, 16, "%.1f", mass_xy); // set the buffer to current Mass XY
		}
	}
	else // input mode
	{
		// Text box to input value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##Mass XY", massxy_str, IM_ARRAYSIZE(massxy_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			mass_xy = atof(massxy_str);
			// set the load value to input value
			// deformation_scale_max = defscale_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			massxy_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text(" %.1f", mass_xy);

	// ________________________________________________________________________________________________________________________________



	ImGui::Spacing();

	// Add Point Mass
	if (is_add_pointmass == true)
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
	}
	else
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
	}

	if (ImGui::Button("Add Point Mass"))
	{
		is_add_pointmass = !is_add_pointmass;
	}

	// Text to inform user whether Add constrain is added or not
	ImGui::SameLine(); // move cursor to the same line
	if (is_add_pointmass == true)
	{
		ImGui::Text("Click on the node to add the point mass");
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
		is_add_pointmass = false; // clear the add point mass flag
	}

	ImGui::End();
}
