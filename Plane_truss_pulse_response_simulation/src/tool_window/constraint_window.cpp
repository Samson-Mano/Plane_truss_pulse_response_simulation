#include "constraint_window.h"


constraint_window::constraint_window()
{
	// Empty constructor
}

constraint_window::~constraint_window()
{
	// Empty destructor
}

void constraint_window::init()
{
	// Load the texture from a file
	std::string img_path = "./resources/images/frame_supports.png";

	// Bind the constraint image
	LoadTextureFromFile(img_path.c_str(), cnst_image);
}

void constraint_window::render_window()
{
	if (is_show_window == false)
		return;

	ImGui::Begin("Constraints");

	// Add constraint input controls
// Option to select the types of support
// Define an array of options
	const int options_count = 4;
	const char* options[] = { "Fixed end support", "Fixed roller support",
							  "Pin end support", "Pin roller support" };

	// Define a string to hold the label for the popup select button
	std::string popupLabel = "Support: ";

	if (ImGui::Button((popupLabel + options[constraint_type]).c_str()))
	{
		ImGui::OpenPopup("Select an option");
	}

	if (ImGui::BeginPopup("Select an option"))
	{
		ImGui::Text("- Constraint Type -");
		ImGui::Separator();

		for (int i = 0; i < options_count; i++)
		{
			if (ImGui::Selectable(options[i], constraint_type == i))
			{
				constraint_type = i;
			}
		}

		ImGui::EndPopup();
	}
	//_________________________________________________________________________________________

	// Input box to give input via text
	static bool input_mode = false;
	static char angle_str[8] = ""; // buffer to store input angle string
	static float angle_input = static_cast<float>(constraint_angle); // buffer to store input angle value

	// Button to switch to input mode
	if (!input_mode)
	{
		if (ImGui::Button("Input Angle"))
		{
			input_mode = true;
			snprintf(angle_str, 8, "%.1f", angle_input); // set the buffer to current angle
		}
	}
	else // input mode
	{
		// Text box to input angle value
		ImGui::SetNextItemWidth(60.0f);
		if (ImGui::InputText("##InputAngle", angle_str, IM_ARRAYSIZE(angle_str), ImGuiInputTextFlags_CharsDecimal))
		{
			// convert the input string to float
			angle_input = static_cast<float>(atof(angle_str));
			// limit the value to 0 - 360 range
			angle_input = fmaxf(0.0f, fminf(angle_input, 360.0f));
			// set the angle to input value
			constraint_angle = angle_input;
		}

		// Button to switch back to slider mode
		ImGui::SameLine();
		if (ImGui::Button("OK"))
		{
			input_mode = false;
		}
	}

	angle_input = static_cast<float>(constraint_angle);
	// Slider for angle
	ImGui::Text("Angle");
	ImGui::SameLine();
	ImGui::SliderFloat("Degrees", &angle_input, 0.0f, 360.0f, "%.1f");

	constraint_angle = angle_input;
	//_________________________________________________________________________________________

	// add some vertical spacing
	ImGui::Spacing();

	// Add Constraint
	if (is_add_constraint == true)
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.8f, 0.4f, 1.0f)); // brighter green color
	}
	else
	{
		ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5f, 0.5f, 0.5f, 1.0f)); // default color
	}

	if (ImGui::Button("Add Constraint"))
	{
		is_add_constraint = !is_add_constraint;
	}

	// Text to inform user whether Add constrain is added or not
	ImGui::SameLine(); // move cursor to the same line
	if (is_add_constraint == true)
	{
		ImGui::Text("Click on the nodes to add constraint");
	}
	else
	{
		ImGui::Text("<= Click here to start adding constraint");
	}

	ImGui::PopStyleColor(1);  // Revert back to default style

	//_________________________________________________________________________________________


	// Close button
	if (ImGui::Button("Close"))
	{
		is_add_constraint = false;
		is_show_window = false; // set the flag to close the window
	}
	//_________________________________________________________________________________________


	// Render reference image
	draw_support();

	ImGui::End();
}

void constraint_window::draw_support()
{
	// Draw the fixed end support image
	if (cnst_image.is_loaded == true)
	{
		ImVec2 window_pos = ImGui::GetCursorScreenPos();
		ImVec2 window_size = ImGui::GetContentRegionAvail();
		ImDrawList* draw_list = ImGui::GetWindowDrawList();


		draw_list->AddRect(window_pos, ImVec2(window_pos.x + window_size.x, window_pos.y + window_size.y), ImColor(255, 255, 255));

		// Initialize zero
		ImVec2 img_pos_top_left(0, 0);
		ImVec2 img_pos_top_right(0, 0);
		ImVec2 img_pos_bot_right(0, 0);
		ImVec2 img_pos_bot_left(0, 0);

		bool draw_img = get_image_min_max_coord(window_pos, window_size, img_pos_top_left, img_pos_top_right,
			img_pos_bot_right, img_pos_bot_left, constraint_angle);

		if (draw_img == true)
		{
			ImVec2 uv_top_left;
			ImVec2 uv_top_right;
			ImVec2 uv_bot_right;
			ImVec2 uv_bot_left;

			if (constraint_type == 0)
			{
				// Draw fixed support
				uv_top_left = ImVec2(0.0f, 0.5f);
				uv_top_right = ImVec2(0.5f, 0.5f);
				uv_bot_right = ImVec2(0.5f, 1.0f);
				uv_bot_left = ImVec2(0.0f, 1.0f);
			}
			else if (constraint_type == 1)
			{
				// Draw fixed roller support
				uv_top_left = ImVec2(0.5f, 0.5f);
				uv_top_right = ImVec2(1.0f, 0.5f);
				uv_bot_right = ImVec2(1.0f, 1.0f);
				uv_bot_left = ImVec2(0.5f, 1.0f);
			}
			else if (constraint_type == 2)
			{
				// Draw pin support
				uv_top_left = ImVec2(0.0f, 0.0f);
				uv_top_right = ImVec2(0.5f, 0.0f);
				uv_bot_right = ImVec2(0.5f, 0.5f);
				uv_bot_left = ImVec2(0.0f, 0.5f);
			}
			else if (constraint_type == 3)
			{
				// Draw pin roller support
				uv_top_left = ImVec2(0.5f, 0.0f);
				uv_top_right = ImVec2(1.0f, 0.0f);
				uv_bot_right = ImVec2(1.0f, 0.5f);
				uv_bot_left = ImVec2(0.5f, 0.5f);
			}


			draw_list->AddImageQuad((void*)(intptr_t)cnst_image.image_texture_ID, img_pos_top_left, img_pos_top_right,
				img_pos_bot_right, img_pos_bot_left,
				uv_top_left,uv_top_right,uv_bot_right,uv_bot_left);
		}
	}
}

void constraint_window::LoadTextureFromFile(const char* filename, cnst_image_data& cnst_image)
{
	// Private function to Load Texture From File
	// Simple helper function to load an image into a OpenGL texture with common settings

	// Load from file only once
	if (!cnst_image.is_loaded)
	{
		// Load image data from file
		unsigned char* image_data = stbi_load(filename, &cnst_image.image_width, &cnst_image.image_height, NULL, 4);
		if (image_data == NULL)
		{
			cnst_image.is_loaded = false;
			return;
		}

		// Create OpenGL texture
		glGenTextures(1, &cnst_image.image_texture_ID);
		glBindTexture(GL_TEXTURE_2D, cnst_image.image_texture_ID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, cnst_image.image_width, cnst_image.image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
		stbi_image_free(image_data);

		cnst_image.is_loaded = true;
	}
}

bool constraint_window::get_image_min_max_coord(ImVec2& window_pos, ImVec2& window_size, ImVec2& img_pos_top_left,
	ImVec2& img_pos_top_right, ImVec2& img_pos_bot_right, ImVec2& img_pos_bot_left, double& orientation)
{
	// Get the image position based on rotation angle
	if (window_size.x < 50 || window_size.y < 50)
		return false;

	double rect_min_size = std::min(window_size.x, window_size.y);
	const double size_factor = 0.415f;


	// Rectangle origin
	glm::vec2 rect_window_origin= glm::vec2((window_pos.x + window_pos.x + window_size.x) * 0.5f,
		(window_pos.y + window_pos.y + window_size.y) * 0.5f);

	// Corners of the image with origin (0,0)
	glm::vec2 top_left = glm::vec2(-(rect_min_size * size_factor), (rect_min_size * size_factor));
	glm::vec2 top_right = glm::vec2((rect_min_size * size_factor), (rect_min_size * size_factor));
	glm::vec2 bot_right = glm::vec2((rect_min_size * size_factor), -(rect_min_size * size_factor));
	glm::vec2 bot_left = glm::vec2(-(rect_min_size * size_factor), -(rect_min_size * size_factor));

	//___________________________________________________________________________________________________________
	float radians = ((static_cast<float>(orientation)+90.0f) * 3.14159365f) / 180.0f; // convert degrees to radians
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

