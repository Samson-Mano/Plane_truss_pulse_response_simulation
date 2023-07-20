#include "load_window.h"


load_window::load_window()
{
	// Empty constructor
}

load_window::~load_window()
{
	// Empty destructor
}

void load_window::init()
{
	// Load the texture from a file
	std::string img_path = "./resources/images/load_img.png";

	// Bind the constraint image
	LoadTextureFromFile(img_path.c_str(), load_image);
}

void load_window::render_window()
{
	if (is_show_window == false)
		return;

	ImGui::Begin("Loads");

	// Add Load input controls
	// Input box to give input via text
	static bool loadval_input_mode = false;
	static char load_str[16] = ""; // buffer to store input load string
	static float load_input = static_cast<float>(load_amplitude); // buffer to store input load value

	// Button to switch to input mode
	if (!loadval_input_mode)
	{
		if (ImGui::Button("Input Load"))
		{
			loadval_input_mode = true;
			snprintf(load_str, 16, "%.2f", load_input); // set the buffer to current load value
		}
	}
	else // input mode
	{
		// Text box to input load value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##InputLoad", load_str, IM_ARRAYSIZE(load_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			load_input = static_cast<float>(atof(load_str));
			// set the load value to input value
			load_amplitude = load_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			loadval_input_mode = false;
		}
	}

	// Text for load value
	ImGui::SameLine();
	ImGui::Text("Load = %.2f", load_input);

	//_________________________________________________________________________________________
	// Input box to give input via text
	static bool loadstarttime_input_mode = false;
	static char loadstarttime_str[16] = ""; // buffer to store input load string
	static float loadstarttime_input = static_cast<float>(load_start_time); // buffer to store input load start time

	// Button to switch to input mode
	if (!loadstarttime_input_mode)
	{
		if (ImGui::Button("Input Start Time"))
		{
			loadstarttime_input_mode = true;
			snprintf(loadstarttime_str, 16, "%.3f", loadstarttime_input); // set the buffer to current load start time
		}
	}
	else // input mode
	{
		// Text box to input load start time
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##InputStartTime", loadstarttime_str, IM_ARRAYSIZE(loadstarttime_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			loadstarttime_input = static_cast<float>(atof(loadstarttime_str));
			// set the load start time to input value
			load_start_time = loadstarttime_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			loadstarttime_input_mode = false;
		}
	}

	// Text for load start time
	ImGui::SameLine();
	ImGui::Text("Start Time = %.3f", loadstarttime_input);

	//_________________________________________________________________________________________
	// Input box to give input via text
	static bool loadendtime_input_mode = false;
	static char loadendtime_str[16] = ""; // buffer to store input load string
	static float loadendtime_input = static_cast<float>(load_end_time); // buffer to store input load End Time

	// Button to switch to input mode
	if (!loadendtime_input_mode)
	{
		if (ImGui::Button("Input End Time"))
		{
			loadendtime_input_mode = true;
			snprintf(loadendtime_str, 16, "%.3f", loadendtime_input); // set the buffer to current load End Time
		}
	}
	else // input mode
	{
		// Text box to input load End Time
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##InputEndTime", loadendtime_str, IM_ARRAYSIZE(loadendtime_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to int
			loadendtime_input = static_cast<float>(atof(loadendtime_str));
			// set the load End Time to input value
			if (loadendtime_input > load_start_time)
			{
				load_end_time = loadendtime_input;
			}
			else
			{
				loadendtime_input = static_cast<float>(load_end_time);
			}
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			loadendtime_input_mode = false;
		}
	}

	// Text for load End Time
	ImGui::SameLine();
	ImGui::Text("End Time = %.3f", loadendtime_input);

	//_________________________________________________________________________________________
	// Display the Load frequency
	
	ImGui::Text("Load frequency = %.3f Hz",
		1.0f / (loadendtime_input - loadstarttime_input));


	//_________________________________________________________________________________________
	// Input box to give input via text
	static bool input_mode = false;
	static char angle_str[8] = ""; // buffer to store input angle string
	static float angle_input = static_cast<float>(load_angle); // buffer to store input angle value

	// Button to switch to input mode
	if (!input_mode)
	{
		if (ImGui::Button("Input Angle"))
		{
			input_mode = true;
			snprintf(angle_str, 8, "%.2f", angle_input); // set the buffer to current angle
		}
	}
	else // input mode
	{
		// Text box to input angle value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##InputAngle", angle_str, IM_ARRAYSIZE(angle_str), ImGuiInputTextFlags_CharsDecimal)) // ImGuiInputTextFlags_CharsDecimal
		{
			// convert the input string to float
			angle_input = static_cast<float>(atof(angle_str));
			// limit the value to 0 - 360 range
			angle_input = fmaxf(0.0f, fminf(angle_input, 360.0f));
			// set the angle to input value
			load_angle = angle_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			input_mode = false;
		}
	}

	// Slider for angle
	angle_input = static_cast<float>(load_angle);

	ImGui::Text("Angle");
	ImGui::SameLine();
	ImGui::SliderFloat("Degrees", &angle_input, 0.0f, 360.0f, "%.2f");

	load_angle = angle_input;
	//_________________________________________________________________________________________

	// add some vertical spacing
	ImGui::Spacing();

	// Add Constraint
	if (is_add_load == true)
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
	}
	else
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
	}

	if (ImGui::Button("Add Load"))
	{
		is_add_load = !is_add_load;
	}

	// Text to inform user whether Add constrain is added or not
	ImGui::SameLine(); // move cursor to the same line
	if (is_add_load == true)
	{
		ImGui::Text("Click on the elements to add load");
	}
	else
	{
		ImGui::Text("<= Click here to start adding load");
	}

	ImGui::PopStyleColor(1);  // Revert back to default style

	//_________________________________________________________________________________________


	// Close button
	if (ImGui::Button("Close"))
	{
		is_add_load = false;
		is_show_window = false; // set the flag to close the window
	}
	//_________________________________________________________________________________________


	// Render reference image
	draw_load();

	ImGui::End();
}

void load_window::draw_load()
{
	// Draw the pin support image
	if (load_image.is_loaded == true)
	{
		ImVec2 window_pos = ImGui::GetCursorScreenPos();
		ImVec2 window_size = ImGui::GetContentRegionAvail();
		ImDrawList* draw_list = ImGui::GetWindowDrawList();


		draw_list->AddRect(window_pos, ImVec2(window_pos.x + window_size.x, window_pos.y + window_size.y), ImColor(255, 255, 255));

		ImVec2 img_pos_top_left(0, 0);
		ImVec2 img_pos_top_right(0, 0);
		ImVec2 img_pos_bot_right(0, 0);
		ImVec2 img_pos_bot_left(0, 0);

		double orientation = load_angle;

		if (load_amplitude < 0)
		{
			orientation = orientation + 180.0;
		}

		bool draw_img = get_image_min_max_coord(window_pos, window_size, img_pos_top_left, img_pos_top_right, img_pos_bot_right, img_pos_bot_left, orientation);

		if (draw_img == true)
		{
			draw_list->AddImageQuad((void*)(intptr_t)load_image.image_texture_ID, img_pos_top_left, img_pos_top_right, img_pos_bot_right, img_pos_bot_left);
		}
	}
}

void load_window::LoadTextureFromFile(const char* filename, load_image_data& load_image)
{
	// Private function to Load Texture From File
	// Simple helper function to load an image into a OpenGL texture with common settings

	// Load from file only once
	if (!load_image.is_loaded)
	{
		// Load image data from file
		unsigned char* image_data = stbi_load(filename, &load_image.image_width, &load_image.image_height, NULL, 4);
		if (image_data == NULL)
		{
			load_image.is_loaded = false;
			return;
		}

		// Create OpenGL texture
		glGenTextures(1, &load_image.image_texture_ID);
		glBindTexture(GL_TEXTURE_2D, load_image.image_texture_ID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, load_image.image_width, load_image.image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
		stbi_image_free(image_data);

		load_image.is_loaded = true;
	}
}

bool load_window::get_image_min_max_coord(ImVec2& window_pos, ImVec2& window_size, ImVec2& img_pos_top_left,
	ImVec2& img_pos_top_right, ImVec2& img_pos_bot_right, ImVec2& img_pos_bot_left, double& orientation)
{
	// Get the image position based on rotation angle
	if (window_size.x < 50 || window_size.y < 50)
		return false;

	// Rectangle origin
	glm::vec2 rect_window_origin = glm::vec2((window_pos.x + window_pos.x + window_size.x) * 0.5f,
		(window_pos.y + window_pos.y + window_size.y) * 0.5f);

	float rect_min_size = std::min(window_size.x, window_size.y);

	const float size_factor = 0.415f;

	// Corners of the image with origin (0,0)
	glm::vec2 top_left = glm::vec2(-(rect_min_size * size_factor), (rect_min_size * size_factor));
	glm::vec2 top_right = glm::vec2((rect_min_size * size_factor), (rect_min_size * size_factor));
	glm::vec2 bot_right = glm::vec2((rect_min_size * size_factor), -(rect_min_size * size_factor));
	glm::vec2 bot_left = glm::vec2(-(rect_min_size * size_factor), -(rect_min_size * size_factor));

	float radians = ((static_cast<float>(orientation) + 90.0f) * 3.14159365f) / 180.0f; // convert degrees to radians
	float cos_theta = cos(radians);
	float sin_theta = sin(radians);

	// Rotated point of the corners
	ImVec2 rotated_pt_top_left = ImVec2((top_left.x * cos_theta) - (top_left.y * sin_theta),
		(top_left.x * sin_theta) + (top_left.y * cos_theta));

	ImVec2 rotated_pt_top_right = ImVec2((top_right.x * cos_theta) - (top_right.y * sin_theta),
		(top_right.x * sin_theta) + (top_right.y * cos_theta));

	ImVec2 rotated_pt_bot_right = ImVec2((bot_right.x * cos_theta) - (bot_right.y * sin_theta),
		(bot_right.x * sin_theta) + (bot_right.y * cos_theta));

	ImVec2 rotated_pt_bot_left = ImVec2((bot_left.x * cos_theta) - (bot_left.y * sin_theta),
		(bot_left.x * sin_theta) + (bot_left.y * cos_theta));

	// Set the corner points after rotating and translating with the origin
	img_pos_top_left = ImVec2(rect_window_origin.x + rotated_pt_top_left.x, rect_window_origin.y + rotated_pt_top_left.y);
	img_pos_top_right = ImVec2(rect_window_origin.x + rotated_pt_top_right.x, rect_window_origin.y + rotated_pt_top_right.y);
	img_pos_bot_right = ImVec2(rect_window_origin.x + rotated_pt_bot_right.x, rect_window_origin.y + rotated_pt_bot_right.y);
	img_pos_bot_left = ImVec2(rect_window_origin.x + rotated_pt_bot_left.x, rect_window_origin.y + rotated_pt_bot_left.y);

	return true;
}
