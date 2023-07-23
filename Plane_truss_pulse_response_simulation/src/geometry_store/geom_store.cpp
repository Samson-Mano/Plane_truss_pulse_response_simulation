#include "geom_store.h"

geom_store::geom_store()
{
	// Empty constructor
}

geom_store::~geom_store()
{
	// Empty Destructor
}

void geom_store::init(options_window* op_window, material_window* mat_window,
	modal_analysis_window* sol_modal_window, pulse_response_window* sol_pulse_window)
{
	// Initialize
	// Initialize the geometry parameters
	geom_param.init();

	is_geometry_set = false;
	is_modal_analysis_complete = false;

	// Initialize the model nodes and lines
	model_nodes.init(&geom_param);
	model_lineelements.init(&geom_param);
	model_constarints.init(&geom_param);
	model_loads.init(&geom_param);
	model_ptmass.init(&geom_param);
	model_inlcond.init(&geom_param);

	// Initialize the modal analysis result nodes and lines
	modal_results.clear_data();
	modal_result_nodes.init(&geom_param);
	modal_result_lineelements.init(&geom_param);

	// Initialize the pulse analysis result nodes and lines
	pulse_response_result.clear_results();
	pulse_result_nodes.init(&geom_param);
	pulse_result_lineelements.init(&geom_param);

	// Add the window pointers
	this->op_window = op_window;
	this->mat_window = mat_window;
	this->sol_modal_window = sol_modal_window;
	this->sol_pulse_window = sol_pulse_window;
}

void geom_store::fini()
{
	// Deinitialize
	is_geometry_set = false;

}

void geom_store::read_varai2d(std::ifstream& input_file)
{
	// Read the varai2D
	// Read the entire file into a string
	std::string file_contents((std::istreambuf_iterator<char>(input_file)),
		std::istreambuf_iterator<char>());

	// Split the string into lines
	std::istringstream iss(file_contents);
	std::string line;
	std::vector<std::string> lines;
	while (std::getline(iss, line))
	{
		lines.push_back(line);
	}

	int j = 0, i = 0;

	// Create a temporary variable to store the nodes
	nodes_list_store model_nodes;
	model_nodes.init(&geom_param);

	// Create a temporary variable to store the lines
	elementline_list_store model_lineelements;
	model_lineelements.init(&geom_param);

	// Process the lines
	while (j < lines.size())
	{
		std::cout << "Line: " << lines[j] << std::endl;
		// Check for the start of nodes input
		if (lines[j].find("[+] End Points") != std::string::npos)
		{
			int num_nodes;
			// Read the number of nodes
			std::stringstream ss(lines[j]);
			std::string token;
			std::getline(ss, token, ','); // Get the string "[+] End Points"
			std::getline(ss, token, ','); // Get the number of nodes as a string
			num_nodes = std::stoi(token) + j; // Convert the string to an integer

			// Read and store the nodes
			for (i = j; i < num_nodes; i++)
			{
				int node_id;
				double x, y;

				std::stringstream ss(lines[i + 1]);
				std::string token;

				std::getline(ss, token, ','); // read the node ID
				node_id = std::stoi(token);

				std::getline(ss, token, ','); // read the x-coordinate
				x = std::stod(token);

				std::getline(ss, token, ','); // read the y-coordinate
				y = std::stod(token);

				// Add to node store list
				glm::vec2 node_pt = glm::vec2(x, y);
				model_nodes.add_node(node_id, node_pt);

				j++;
			}
		}
		// Check for the start of lines input
		else if (lines[j].find("[+] Lines") != std::string::npos) {
			int num_lines;
			// Read the number of nodes
			std::stringstream ss(lines[j]);
			std::string token;
			std::getline(ss, token, ','); // Get the string "[+] Lines"
			std::getline(ss, token, ','); // Get the number of nodes as a string
			num_lines = std::stoi(token) + j; // Convert the string to an integer

			// Read and store the lines
			for (i = j; i < num_lines; i++)
			{
				int line_id, start_node_id, end_node_id;
				std::stringstream ss(lines[i + 1]);
				std::string token;

				std::getline(ss, token, ','); // read the line ID
				line_id = std::stoi(token);

				std::getline(ss, token, ','); // read the start node ID
				start_node_id = std::stoi(token);

				std::getline(ss, token, ','); // read the end node ID
				end_node_id = std::stoi(token);

				// Create lines_store object using references to startNode and endNode
				int mat_id = 0;
				model_lineelements.add_elementline(line_id, &model_nodes.nodeMap[start_node_id], &model_nodes.nodeMap[end_node_id], mat_id);

				j++;
			}
		}

		// iterate line
		j++;
	}

	if (model_nodes.node_count < 1 || model_lineelements.elementline_count < 1)
	{
		// No elements added
		return;
	}


	// add a default material to the material list
	double dia_m = 87.404; // mm
	const double m_pi = 3.14159265358979323846;
	double cs_area = (m_pi * dia_m * dia_m) / 4.0f;
	double second_moment_of_area = (m_pi * dia_m * dia_m * dia_m * dia_m) / 64.0f;

	material_data inpt_material;
	inpt_material.material_id = 0; // Get the material id
	inpt_material.material_name = "Default material"; //Default material name
	inpt_material.mat_density = 7.83 * std::pow(10, -9); // tons/mm3
	inpt_material.youngs_mod = 2.07 * std::pow(10, 5); //  MPa
	inpt_material.cs_area = cs_area; // mm2  (pi * D^2)/4

	// Add to materail list
	mat_window->material_list.clear();
	mat_window->material_list[inpt_material.material_id] = inpt_material;

	nodeconstraint_list_store model_constarints;
	model_constarints.init(&geom_param);

	nodeload_list_store model_loads;
	model_loads.init(&geom_param);

	nodepointmass_list_store model_ptmass;
	model_ptmass.init(&geom_param);

	nodeinlcond_list_store model_inlcond;
	model_inlcond.init(&geom_param);

	// Re-instantitize geom_store object using the nodeMap and lineMap
	create_geometry(model_nodes, model_lineelements, model_constarints, model_loads, model_ptmass, model_inlcond);
}

void geom_store::read_dxfdata(std::ostringstream& input_data)
{
	// Read the data from string
	std::string inputStr = input_data.str();
	std::stringstream ss(inputStr);

	std::string temp;
	std::vector<std::string> lines;
	while (std::getline(ss, temp))
	{
		lines.push_back(temp);
	}

	int j = 0, i = 0;

	// Create a temporary variable to store the nodes
	nodes_list_store model_nodes;
	model_nodes.init(&geom_param);

	// Create a temporary variable to store the lines
	elementline_list_store model_lineelements;
	model_lineelements.init(&geom_param);

	// Process the lines
	while (j < lines.size())
	{
		std::string line = lines[j];
		std::string type = line.substr(0, 4);  // Extract the first 4 characters of the line

		// Split the line into comma-separated fields
		std::istringstream iss(line);
		std::string field;
		std::vector<std::string> fields;
		while (std::getline(iss, field, ','))
		{
			fields.push_back(field);
		}

		if (type == "node")
		{
			// Read the nodes
			int node_id = std::stoi(fields[1]); // node ID
			double x = std::stod(fields[2]); // Node coordinate x
			double y = std::stod(fields[3]); // Node coordinate y

			// Add to node Map
			glm::vec2 node_pt = glm::vec2(x, y);
			model_nodes.add_node(node_id, node_pt);
		}
		else if (type == "line")
		{
			int line_id = std::stoi(fields[1]); // line ID
			int start_node_id = std::stoi(fields[2]); // line id start node
			int end_node_id = std::stoi(fields[3]); // line id end node
			int material_id = std::stoi(fields[4]); // materail ID of the line

			// Add to line Map (Note that Nodes needed to be added before the start of line addition !!!!)
			int mat_id = 0;
			model_lineelements.add_elementline(line_id, &model_nodes.nodeMap[start_node_id], &model_nodes.nodeMap[end_node_id], mat_id);
		}

		// Iterate line
		j++;
	}

	if (model_nodes.node_count < 1 || model_lineelements.elementline_count < 1)
	{
		// No elements added
		return;
	}


	// add a default material to the material list
	double dia_m = 87.404; // mm
	const double m_pi = 3.14159265358979323846;
	double cs_area = (m_pi * dia_m * dia_m) / 4.0f;

	material_data inpt_material;
	inpt_material.material_id = 0; // Get the material id
	inpt_material.material_name = "Default material"; //Default material name
	inpt_material.mat_density = 7.83 * std::pow(10, -9); // tons/mm3
	inpt_material.youngs_mod = 2.07 * std::pow(10, 5); //  MPa
	inpt_material.cs_area = cs_area; // mm2

	// Add to materail list
	mat_window->material_list.clear();
	mat_window->material_list[inpt_material.material_id] = inpt_material;

	nodeconstraint_list_store model_constarints;
	model_constarints.init(&geom_param);

	nodeload_list_store model_loads;
	model_loads.init(&geom_param);

	nodepointmass_list_store model_ptmass;
	model_ptmass.init(&geom_param);

	nodeinlcond_list_store model_inlcond;
	model_inlcond.init(&geom_param);

	// Re-instantitize geom_store object using the nodeMap and lineMap
	create_geometry(model_nodes, model_lineelements, model_constarints, model_loads, model_ptmass, model_inlcond);
}


void geom_store::read_rawdata(std::ifstream& input_file)
{
	// Read the Raw Data
	// Read the entire file into a string
	std::string file_contents((std::istreambuf_iterator<char>(input_file)),
		std::istreambuf_iterator<char>());

	// Split the string into lines
	std::istringstream iss(file_contents);
	std::string line;
	std::vector<std::string> lines;
	while (std::getline(iss, line))
	{
		lines.push_back(line);
	}

	int j = 0, i = 0;

	// Create a temporary variable to store the nodes
	nodes_list_store model_nodes;
	model_nodes.init(&geom_param);

	// Create a temporary variable to store the lines
	elementline_list_store model_lineelements;
	model_lineelements.init(&geom_param);

	// Constraint data store
	nodeconstraint_list_store model_constarints;
	model_constarints.init(&geom_param);

	// Load data store
	nodeload_list_store model_loads;
	model_loads.init(&geom_param);

	// Point mass data store
	nodepointmass_list_store model_ptmass;
	model_ptmass.init(&geom_param);

	// Initial condition data store
	nodeinlcond_list_store model_inlcond;
	model_inlcond.init(&geom_param);

	// Material data list
	std::unordered_map<int, material_data> mat_data;

	// Process the lines
	while (j < lines.size())
	{
		std::string line = lines[j];
		std::string type = line.substr(0, 4);  // Extract the first 4 characters of the line

		// Split the line into comma-separated fields
		std::istringstream iss(line);
		std::string field;
		std::vector<std::string> fields;
		while (std::getline(iss, field, ','))
		{
			fields.push_back(field);
		}

		if (type == "node")
		{
			// Read the nodes
			int node_id = std::stoi(fields[1]); // node ID
			double x = std::stod(fields[2]); // Node coordinate x
			double y = std::stod(fields[3]); // Node coordinate y

			// Add to node Map
			glm::vec2 node_pt = glm::vec2(x, y);
			model_nodes.add_node(node_id, node_pt);
		}
		else if (type == "line")
		{
			int line_id = std::stoi(fields[1]); // line ID
			int start_node_id = std::stoi(fields[2]); // line id start node
			int end_node_id = std::stoi(fields[3]); // line id end node
			int material_id = std::stoi(fields[4]); // materail ID of the line

			// Add to line Map (Note that Nodes needed to be added before the start of line addition !!!!)
			model_lineelements.add_elementline(line_id, &model_nodes.nodeMap[start_node_id], &model_nodes.nodeMap[end_node_id], material_id);
		}
		else if (type == "cnst")
		{
			int cnst_nd_id = std::stoi(fields[1]); // constraint node ID
			int cnst_type = std::stoi(fields[2]); // constraint type 
			double cnst_angle = std::stod(fields[3]); // constraint angle

			// Add to constraint map
			model_constarints.add_constraint(cnst_nd_id, model_nodes.nodeMap[cnst_nd_id].node_pt, cnst_type, cnst_angle);
		}
		else if (type == "load")
		{
			int load_nd_id = std::stoi(fields[1]); // load node ID
			double load_val = std::stod(fields[2]); // load value
			double load_angle = std::stod(fields[3]); // load angle
			double load_start_time = std::stod(fields[4]); // load start time
			double load_end_time = std::stod(fields[5]); // load end time
			double load_loc_x = std::stod(fields[6]); // load loc x
			double load_loc_y = std::stod(fields[7]); // load loc y

			glm::vec2 load_loc = glm::vec2(load_loc_x, load_loc_y);

			// Add to load map
			model_loads.add_load(load_nd_id, load_loc, load_start_time, load_end_time, load_val, load_angle);
		}
		else if (type == "ptms")
		{
			int ptm_nd_id = std::stoi(fields[1]); // load node ID
			double ptm_x = std::stod(fields[2]); // point mass x
			double ptm_y = std::stod(fields[3]); // point mass y

			// Add to point mass map
			model_ptmass.add_pointmass(ptm_nd_id, model_nodes.nodeMap[ptm_nd_id].node_pt, glm::vec2(0), ptm_x, ptm_y, false);
		}
		else if (type == "ilcd")
		{
			int ilcd_nd_id = std::stoi(fields[1]); // initial condition node ID
			double ilcd_displ_x = std::stod(fields[2]); // initial displacement x
			double ilcd_displ_y = std::stod(fields[3]); // initial displacement y
			double ilcd_velo_x = std::stod(fields[4]); // initial velocity x
			double ilcd_velo_y = std::stod(fields[5]); // initial velocity y

			// Add to the initial condition map
			model_inlcond.add_inlcondition(ilcd_nd_id, model_nodes.nodeMap[ilcd_nd_id].node_pt, ilcd_displ_x, ilcd_displ_y, ilcd_velo_x, ilcd_velo_y);
		}
		else if (type == "mtrl")
		{
			// Material data
			material_data inpt_material;
			inpt_material.material_id = std::stoi(fields[1]); // Get the material id
			inpt_material.material_name = fields[2]; // Get the material name
			inpt_material.youngs_mod = std::stod(fields[3]); // Get the material youngs modulus
			inpt_material.mat_density = std::stod(fields[4]); // Get the material density 
			inpt_material.cs_area = std::stod(fields[5]); // Get the material cross section area

			// Add to materail list
			mat_data[inpt_material.material_id] = inpt_material;
		}

		// Iterate line
		j++;
	}


	// Data loaded create the geometry
	if (model_nodes.node_count < 1 || model_lineelements.elementline_count < 1)
	{
		// No elements added
		return;
	}

	//Add the materail list
	mat_window->material_list = mat_data;

	// Re-instantitize geom_store object using the nodeMap and lineMap
	create_geometry(model_nodes, model_lineelements, model_constarints, model_loads, model_ptmass, model_inlcond);
}

void geom_store::write_rawdata(std::ofstream& output_file)
{
	// Write all the nodes
	for (auto& node : model_nodes.nodeMap)
	{
		// Print the node details
		node_store nd_val = node.second;

		output_file << "node, "
			<< nd_val.node_id << ", "
			<< nd_val.node_pt.x << ", "
			<< nd_val.node_pt.y << std::endl;
	}

	// Write all the lines
	for (auto& line : model_lineelements.elementlineMap)
	{
		// Print the line details
		elementline_store ln_val = line.second;

		output_file << "line, "
			<< ln_val.line_id << ", "
			<< ln_val.startNode->node_id << ", "
			<< ln_val.endNode->node_id << ", "
			<< ln_val.material_id << std::endl;
	}

	// Write all the constraints
	for (auto& cnst : model_constarints.constraintMap)
	{
		// Print the constraint details
		constraint_data cn_val = cnst.second;

		output_file << "cnst, "
			<< cn_val.node_id << ", "
			<< cn_val.constraint_type << ", "
			<< cn_val.constraint_angle << std::endl;
	}

	// Write all the loads
	for (auto& ld : model_loads.loadMap)
	{
		// Print the load details
		load_data ld_val = ld.second;

		output_file << "load, "
			<< ld_val.node_id << ", "  // load node ID
			<< ld_val.load_value << ", " // load value
			<< ld_val.load_angle << ", " // load angle
			<< ld_val.load_start_time << ", " // load start time
			<< ld_val.load_end_time << ", " // load end time
			<< ld_val.load_loc.x << ", " // load loc x
			<< ld_val.load_loc.y << std::endl; // load loc y
	}

	// Write all the point mass
	for (auto& ptmass : model_ptmass.ptmassMap)
	{
		// Print the Point Mass details
		nodepointmass_data ptm = ptmass.second;

		output_file << "ptms, "
			<< ptm.node_id << ", "
			<< ptm.ptmass_x << ", "
			<< ptm.ptmass_y << std::endl;
	}

	// Write all the initial condition data
	for (auto& ilcd_m : model_inlcond.inlcondMap)
	{
		// Print the initial condition details
		nodeinl_condition_data ilcd = ilcd_m.second;

		output_file << "ilcd, "
			<< ilcd.node_id << ", "
			<< ilcd.inl_displacement_x << ", "
			<< ilcd.inl_displacement_y << ", "
			<< ilcd.inl_velocity_x << ", "
			<< ilcd.inl_velocity_y << std::endl;
	}

	// Write all the material property
	for (auto& mat : mat_window->material_list)
	{
		material_data mat_d = mat.second;
		output_file << "mtrl, "
			<< mat_d.material_id << ", "
			<< mat_d.material_name << ", "
			<< mat_d.youngs_mod << ", "
			<< mat_d.mat_density << ", "
			<< mat_d.cs_area << std::endl;
	}
}


void geom_store::update_WindowDimension(const int& window_width, const int& window_height)
{
	// Update the window dimension
	this->geom_param.window_width = window_width;
	this->geom_param.window_height = window_height;

	if (is_geometry_set == true)
	{
		// Update the model matrix
		update_model_matrix();
		// !! Zoom to fit operation during window resize is handled in mouse event class !!
	}
}

void geom_store::update_model_matrix()
{
	// Set the model matrix for the model shader
	// Find the scale of the model (with 0.9 being the maximum used)
	int max_dim = geom_param.window_width > geom_param.window_height ? geom_param.window_width : geom_param.window_height;

	double normalized_screen_width = 1.8f * (static_cast<double>(geom_param.window_width) / static_cast<double>(max_dim));
	double normalized_screen_height = 1.8f * (static_cast<double>(geom_param.window_height) / static_cast<double>(max_dim));


	geom_param.geom_scale = std::min(normalized_screen_width / geom_param.geom_bound.x,
		normalized_screen_height / geom_param.geom_bound.y);

	// Translation
	glm::vec3 geom_translation = glm::vec3(-1.0f * (geom_param.max_b.x + geom_param.min_b.x) * 0.5f * geom_param.geom_scale,
		-1.0f * (geom_param.max_b.y + geom_param.min_b.y) * 0.5f * geom_param.geom_scale,
		0.0f);

	glm::mat4 g_transl = glm::translate(glm::mat4(1.0f), geom_translation);

	geom_param.modelMatrix = g_transl * glm::scale(glm::mat4(1.0f), glm::vec3(static_cast<float>(geom_param.geom_scale)));

	// Update the model matrix
	model_nodes.update_geometry_matrices(true, false, false, false, false);
	model_lineelements.update_geometry_matrices(true, false, false, false, false);
	model_constarints.update_geometry_matrices(true, false, false, false, false);
	model_loads.update_geometry_matrices(true, false, false, false, false);
	model_ptmass.update_geometry_matrices(true, false, false, false, false);
	model_inlcond.update_geometry_matrices(true, false, false, false, false);

	// Update the modal analysis result matrix
	modal_result_lineelements.update_geometry_matrices(true, false, false, false, false);
	modal_result_nodes.update_geometry_matrices(true, false, false, false, false);

	// Update the pulse analysis result matrix
	pulse_result_lineelements.update_geometry_matrices(true, false, false, false, false);
	pulse_result_nodes.update_geometry_matrices(true, false, false, false, false);
}

void geom_store::update_model_zoomfit()
{
	if (is_geometry_set == false)
		return;

	// Set the pan translation matrix
	geom_param.panTranslation = glm::mat4(1.0f);

	// Set the zoom scale
	geom_param.zoom_scale = 1.0f;

	// Update the zoom scale and pan translation
	model_nodes.update_geometry_matrices(false, true, true, false, false);
	model_lineelements.update_geometry_matrices(false, true, true, false, false);
	model_constarints.update_geometry_matrices(false, true, true, false, false);
	model_loads.update_geometry_matrices(false, true, true, false, false);
	model_ptmass.update_geometry_matrices(false, true, true, false, false);
	model_inlcond.update_geometry_matrices(false, true, true, false, false);

	// Update the modal analysis result matrix
	modal_result_lineelements.update_geometry_matrices(false, true, true, false, false);
	modal_result_nodes.update_geometry_matrices(false, true, true, false, false);

	// Update the pulse analysis result matrix
	pulse_result_lineelements.update_geometry_matrices(false, true, true, false, false);
	pulse_result_nodes.update_geometry_matrices(false, true, true, false, false);
}

void geom_store::update_model_pan(glm::vec2& transl)
{
	if (is_geometry_set == false)
		return;

	// Pan the geometry
	geom_param.panTranslation = glm::mat4(1.0f);

	geom_param.panTranslation[0][3] = -1.0f * transl.x;
	geom_param.panTranslation[1][3] = transl.y;

	// Update the pan translation
	model_nodes.update_geometry_matrices(false, true, false, false, false);
	model_lineelements.update_geometry_matrices(false, true, false, false, false);
	model_constarints.update_geometry_matrices(false, true, false, false, false);
	model_loads.update_geometry_matrices(false, true, false, false, false);
	model_ptmass.update_geometry_matrices(false, true, false, false, false);
	model_inlcond.update_geometry_matrices(false, true, false, false, false);

	// Update the modal analysis result matrix
	modal_result_lineelements.update_geometry_matrices(false, true, false, false, false);
	modal_result_nodes.update_geometry_matrices(false, true, false, false, false);

	// Update the pulse analysis result matrix
	pulse_result_lineelements.update_geometry_matrices(false, true, false, false, false);
	pulse_result_nodes.update_geometry_matrices(false, true, false, false, false);
}

void geom_store::update_model_zoom(double& z_scale)
{
	if (is_geometry_set == false)
		return;

	// Zoom the geometry
	geom_param.zoom_scale = z_scale;

	// Update the Zoom
	model_nodes.update_geometry_matrices(false, false, true, false, false);
	model_lineelements.update_geometry_matrices(false, false, true, false, false);
	model_constarints.update_geometry_matrices(false, false, true, false, false);
	model_loads.update_geometry_matrices(false, false, true, false, false);
	model_ptmass.update_geometry_matrices(false, false, true, false, false);
	model_inlcond.update_geometry_matrices(false, false, true, false, false);

	// Update the modal analysis result matrix
	modal_result_lineelements.update_geometry_matrices(false, false, true, false, false);
	modal_result_nodes.update_geometry_matrices(false, false, true, false, false);

	// Update the pulse analysis result matrix
	pulse_result_lineelements.update_geometry_matrices(false, false, true, false, false);
	pulse_result_nodes.update_geometry_matrices(false, false, true, false, false);
}

void geom_store::update_model_transperency(bool is_transparent)
{
	if (is_geometry_set == false)
		return;

	if (is_transparent == true)
	{
		// Set the transparency value
		geom_param.geom_transparency = 0.2f;
	}
	else
	{
		// remove transparency
		geom_param.geom_transparency = 1.0f;
	}

	// Update the model transparency
	model_nodes.update_geometry_matrices(false, false, false, true, false);
	model_lineelements.update_geometry_matrices(false, false, false, true, false);
	model_constarints.update_geometry_matrices(false, false, false, true, false);
	model_loads.update_geometry_matrices(false, false, false, true, false);
	model_ptmass.update_geometry_matrices(false, false, false, true, false);
	model_inlcond.update_geometry_matrices(false, false, false, true, false);

	// Update the modal analysis result matrix
	// modal_result_lineelements.update_geometry_matrices(false, false, false, true, false);
	// modal_result_nodes.update_geometry_matrices(false, false, false, true, false);
}

void geom_store::set_nodal_constraint(glm::vec2 mouse_click_loc, int& constraint_type, double& constraint_angle, bool is_add)
{
	// geometry is set so check whether node is hit
	int node_hit_id = -1;

	if (is_geometry_set == true)
	{
		// Check whether the node is hit or not
		node_hit_id = model_nodes.is_node_hit(mouse_click_loc);;

		if (node_hit_id != -1)
		{
			// Node is hit
			if (is_add == true)
			{
				// Add constraints
				model_constarints.add_constraint(node_hit_id, model_nodes.nodeMap[node_hit_id].node_pt, constraint_type, constraint_angle);
				model_constarints.set_buffer();
			}
			else
			{
				// remove constraint
				model_constarints.delete_constraint(node_hit_id);
				model_constarints.set_buffer();
			}
		}
	}
}

void geom_store::set_member_load(glm::vec2 mouse_click_loc, double& load_start_time, double& load_end_time,
	double& load_value, double& load_angle, bool is_add)
{
	int node_hit_id = -1;

	if (is_geometry_set == true)
	{
		// Check whether the node is hit or not
		node_hit_id = model_nodes.is_node_hit(mouse_click_loc);;

		if (node_hit_id != -1)
		{
			// node is hit
			if (is_add == true)
			{
				// Get the location of the load (node)
				glm::vec2 load_loc = model_nodes.nodeMap[node_hit_id].node_pt;

				// Add Load
				model_loads.add_load(node_hit_id, load_loc, load_start_time, load_end_time, load_value, load_angle);
				model_loads.set_buffer();
			}
			else
			{
				// remove all the loads on the member
				model_loads.delete_load(node_hit_id);
				model_loads.set_buffer();
			}
		}
	}
}

void geom_store::set_elementline_material(glm::vec2 mouse_click_loc)
{
	// Set the element line material
	int line_hit_id = -1;

	if (is_geometry_set == true)
	{
		// Check whether the line is hit or not
		line_hit_id = model_lineelements.is_line_hit(mouse_click_loc);

		if (line_hit_id != -1)
		{
			// Line Material
			const int& selected_material_option = mat_window->selected_material_option;
			model_lineelements.elementlineMap[line_hit_id].material_id = mat_window->material_list[selected_material_option].material_id;

			// Update the material ID label
			model_lineelements.update_material_id_labels();
		}
	}
}

void geom_store::set_nodal_pointmass(glm::vec2 mouse_click_loc, double& pt_mass_x, double& pt_mass_y, bool is_add)
{
	// Set the nodal point mass
	int node_hit_id = -1;

	if (is_geometry_set == true)
	{
		// Check whether the node is hit or not
		node_hit_id = model_nodes.is_node_hit(mouse_click_loc);;
		if (node_hit_id != -1)
		{
			// Node is hit
			if (is_add == true)
			{
				// Add Point mass
				model_ptmass.add_pointmass(node_hit_id, model_nodes.nodeMap[node_hit_id].node_pt, glm::vec2(0), pt_mass_x, pt_mass_y, false);
				model_ptmass.set_buffer();
			}
			else
			{
				// remove Point mass
				model_ptmass.delete_pointmass(node_hit_id);
				model_ptmass.set_buffer();
			}
		}
	}
}

void geom_store::set_nodal_initialcondition(glm::vec2 mouse_click_loc, double& inl_displ_x, double& inl_displ_y,
	double& inl_velo_x, double& inl_velo_y, bool is_add)
{
	// Set the nodal initial condition
	int node_hit_id = -1;

	if (is_geometry_set == true)
	{
		// Check whether the node is hit or not
		node_hit_id = model_nodes.is_node_hit(mouse_click_loc);;
		if (node_hit_id != -1)
		{
			// Node is hit
			if (is_add == true)
			{
				// Add initial condition
				model_inlcond.add_inlcondition(node_hit_id, model_nodes.nodeMap[node_hit_id].node_pt,inl_displ_x,inl_displ_y,inl_velo_x,inl_velo_y);
				model_inlcond.set_buffer();
			}
			else
			{
				// remove initial condition
				model_inlcond.delete_inlcondition(node_hit_id);
				model_inlcond.set_buffer();
			}
		}
	}
}

void geom_store::paint_geometry()
{
	if (is_geometry_set == false)
		return;

	// Clean the back buffer and assign the new color to it
	glClear(GL_COLOR_BUFFER_BIT);

	// Paint the model
	paint_model();

	// Modal Analysis
	paint_modal_analysis();

	// Pulse Response Analysis
	paint_pulse_analysis();
}

void geom_store::paint_model()
{
	//____________________________________________________________
	// Postprocessing is in progress

	if (sol_modal_window->is_show_window == true && sol_modal_window->show_undeformed_model == false)
	{
		return;
	}

	//____________________________________________________________

	// Paint the model
	model_constarints.paint_constraints();
	model_lineelements.paint_elementlines();
	model_nodes.paint_model_nodes();
	model_loads.paint_loads();
	model_ptmass.paint_pointmass();

	if (op_window->is_show_nodenumber == true)
	{
		// Show model node number
		model_nodes.paint_label_node_ids();
	}

	if (op_window->is_show_nodecoord == true)
	{
		// Show model node coordinate
		model_nodes.paint_label_node_coords();
	}

	if (op_window->is_show_linenumber == true)
	{
		// Show line ID label
		model_lineelements.paint_label_line_ids();
	}

	if (op_window->is_show_linelength == true)
	{
		// Show line length label
		model_lineelements.paint_label_line_lengths();
	}

	if (op_window->is_show_loadvalue == true)
	{
		// Show load value label
		model_loads.paint_load_labels();
		model_ptmass.paint_pointmass_label();

		// Show the initial condition label
		model_inlcond.paint_inlcondition_label();
	}

	if (mat_window->is_show_window == true)
	{
		// Show the materials of line member
		if (mat_window->execute_delete_materialid != -1)
		{
			// Delete material
			update_delete_material(mat_window->execute_delete_materialid);
			mat_window->execute_delete_materialid = -1;
		}
		// Show the material ID
		model_lineelements.paint_lines_material_id();
	}
}

void geom_store::paint_modal_analysis()
{
	// Check closing sequence for modal analysis window
	if (sol_modal_window->execute_close == true)
	{
		// Execute the close sequence
		if (is_modal_analysis_complete == true)
		{
			// Modal analysis is complete remove the transparency for the model
			update_model_transperency(false);
		}

		sol_modal_window->execute_close = false;
	}

	// Check whether the modal analysis solver window is open or not
	if (sol_modal_window->is_show_window == false)
	{
		return;
	}

	// Paint the modal analysis result
	if (is_modal_analysis_complete == true)
	{
		// Change the buffer depending on the selected mode
		if (sol_modal_window->is_mode_selection_changed == true)
		{
			// Update the buffers
			// Modal Line buffer
			modal_result_lineelements.set_buffer(sol_modal_window->selected_modal_option);

			// Modal Node buffer
			modal_result_nodes.set_buffer(sol_modal_window->selected_modal_option);

			sol_modal_window->is_mode_selection_changed = false;
		}

		// Update the deflection scale
		geom_param.normalized_defl_scale = std::abs(sol_modal_window->normailzed_defomation_scale);
		geom_param.defl_scale = sol_modal_window->deformation_scale;

		// Update the deflection scale
		modal_result_lineelements.update_geometry_matrices(false, false, false, false, true);
		modal_result_nodes.update_geometry_matrices(false, false, false, false, true);
		// ______________________________________________________________________________________

		// Paint the modal lines
		modal_result_lineelements.paint_modal_elementlines();

		// Paint the modal nodes
		modal_result_nodes.paint_modal_nodes();

		// Paint result text
		if (sol_modal_window->show_result_text_values == true)
		{
			// Paint the modal result vector
			modal_result_nodes.paint_label_mode_vectors();
		}
	}


	if (sol_modal_window->execute_open == true)
	{
		// Execute the open sequence
		if (is_modal_analysis_complete == true)
		{
			// update the modal window list box
			sol_modal_window->mode_result_str = modal_results.mode_result_str;

			// Set the buffer
			sol_modal_window->is_mode_selection_changed = true;

			// Modal analysis is already complete so set the transparency for the model
			update_model_transperency(true);
		}
		sol_modal_window->execute_open = false;
	}

	if (sol_modal_window->execute_modal_analysis == true)
	{
		// Execute the Modal Analysis
		md_solver.modal_analysis_start(model_nodes,
			model_lineelements,
			model_constarints,
			model_ptmass,
			mat_window->material_list,
			sol_modal_window->is_include_consistent_mass_matrix,
			modal_results,
			modal_result_nodes,
			modal_result_lineelements,
			is_modal_analysis_complete);

		// reset the frequency response and pulse response solution
		is_freq_analysis_complete = false;
		is_pulse_analysis_complete = false;

		// Check whether the modal analysis is complete or not
		if (is_modal_analysis_complete == true)
		{
			// update the modal window list box
			sol_modal_window->mode_result_str = modal_results.mode_result_str;

			// Set the buffer
			sol_modal_window->is_mode_selection_changed = true;

			// Modal analysis is already complete so set the transparency for the model
			update_model_transperency(true);
		}

		sol_modal_window->execute_modal_analysis = false;
	}

}

void geom_store::paint_pulse_analysis()
{
	// Check closing sequence for Pulse response analysis window
	if (sol_pulse_window->execute_close == true)
	{
		// Execute the close sequence
		if (is_pulse_analysis_complete == true)
		{
			// Pulse response is complete (but clear the results anyway beacuse results will be loaded at open)
			sol_pulse_window->pulse_response_analysis_complete = false;

			// Pulse response analysis is complete
			update_model_transperency(false);
		}

		sol_pulse_window->execute_close = false;
	}

	// Check whether the modal analysis solver window is open or not
	if (sol_pulse_window->is_show_window == false)
	{
		return;
	}

	// Paint the pulse analysis result
	if (is_pulse_analysis_complete == true)
	{
		// Update the deflection scale
		geom_param.normalized_defl_scale = 1.0f;
		geom_param.defl_scale = sol_pulse_window->deformation_scale_max;

		// Update the deflection scale
		pulse_result_lineelements.update_geometry_matrices(false, false, false, false, true);
		pulse_result_nodes.update_geometry_matrices(false, false, false, false, true);
		// ______________________________________________________________________________________

		// Paint the pulse lines
		pulse_result_lineelements.paint_pulse_elementlines(sol_pulse_window->time_step);

		// Paint the pulse nodes
		pulse_result_nodes.paint_pulse_nodes(sol_pulse_window->time_step);
	}


	if (sol_pulse_window->execute_open == true)
	{
		// Execute the open sequence
		if (is_modal_analysis_complete == false)
		{
			// Exit the window (when modal analysis is not complete)
			sol_pulse_window->is_show_window = false;
		}
		else
		{
			// Modal analysis Results
			sol_pulse_window->number_of_modes = static_cast<int>(modal_results.eigen_values.size());
			sol_pulse_window->modal_first_frequency = std::sqrt(modal_results.eigen_values.at(0)) / (2.0 * m_pi); // std::sqrt(modal_results.eigen_values[i]) / (2.0 * m_pi);
			sol_pulse_window->modal_end_frequency = std::sqrt(modal_results.eigen_values.at(sol_pulse_window->number_of_modes - 1)) / (2.0 * m_pi);

			// Modal analysis is complete (check whether frequency response analysis is complete or not)
			if (is_pulse_analysis_complete == true)
			{
				// Set the pulse response analysis result
				sol_pulse_window->pulse_response_analysis_complete = true;
				sol_pulse_window->time_interval_atrun = pulse_response_result.time_interval;
				sol_pulse_window->time_step_count = pulse_response_result.time_step_count;

				// Reset the buffers for pulse result nodes and lines
				pulse_result_lineelements.set_buffer();
				pulse_result_nodes.set_buffer();

				// Pulse response analysis is complete
				update_model_transperency(true);
			}

		}
		sol_pulse_window->execute_open = false;
	}

	if (sol_pulse_window->execute_pulse_analysis == true)
	{
		// Execute the Frequency response Analysis
		pulse_analysis_solver pulse_solver;
		pulse_solver.pulse_analysis_start(model_nodes,
			model_lineelements,
			model_constarints,
			model_loads,
			model_ptmass,
			model_inlcond,
			mat_window->material_list,
			sol_modal_window->is_include_consistent_mass_matrix,
			md_solver,
			modal_results,
			sol_pulse_window->total_simulation_time,
			sol_pulse_window->time_interval,
			sol_pulse_window->damping_ratio,
			pulse_response_result,
			pulse_result_nodes,
			pulse_result_lineelements,
			is_pulse_analysis_complete);

		// Check whether the modal analysis is complete or not
		if (is_pulse_analysis_complete == true)
		{
			// Set the pulse response analysis result
			sol_pulse_window->pulse_response_analysis_complete = true;
			sol_pulse_window->time_interval_atrun = pulse_response_result.time_interval;
			sol_pulse_window->time_step_count = pulse_response_result.time_step_count;

			// Reset the buffers for pulse result nodes and lines
			pulse_result_lineelements.set_buffer();
			pulse_result_nodes.set_buffer();

			// Pulse response analysis is complete
			update_model_transperency(true);
		}
		sol_pulse_window->execute_pulse_analysis = false;
	}
}

void geom_store::create_geometry(nodes_list_store& model_nodes, elementline_list_store& model_lineelements,
	nodeconstraint_list_store& model_constarints, nodeload_list_store& model_loads, nodepointmass_list_store& model_ptmass,
	nodeinlcond_list_store& model_inlcond)
{
	// Reinitialize the model geometry
	is_geometry_set = false;
	is_modal_analysis_complete = false;

	this->model_nodes.init(&geom_param);
	this->model_lineelements.init(&geom_param);
	this->model_constarints.init(&geom_param);
	this->model_loads.init(&geom_param);
	this->model_ptmass.init(&geom_param);
	this->model_inlcond.init(&geom_param);

	//________________________________________________

	// Add to model nodes
	for (auto& nd : model_nodes.nodeMap)
	{
		// create a temporary node
		node_store temp_node;
		temp_node = nd.second;

		// Add to the node list
		this->model_nodes.add_node(temp_node.node_id, temp_node.node_pt);
	}

	// Add to model lines
	for (auto& ln : model_lineelements.elementlineMap)
	{
		// create a temporary line element
		elementline_store temp_line;
		temp_line = ln.second;

		// Add to the element list
		this->model_lineelements.add_elementline(temp_line.line_id, &this->model_nodes.nodeMap[temp_line.startNode->node_id],
			&this->model_nodes.nodeMap[temp_line.endNode->node_id], temp_line.material_id);
	}

	// Add to model constraints
	for (auto& cnst : model_constarints.constraintMap)
	{
		// create a temporary constraint
		constraint_data temp_cnst;
		temp_cnst = cnst.second;

		// Add to the constraint list
		this->model_constarints.add_constraint(temp_cnst.node_id, temp_cnst.constraint_loc,
			temp_cnst.constraint_type, temp_cnst.constraint_angle);
	}

	// Add to model loads
	for (auto& load : model_loads.loadMap)
	{
		// create a temporary load
		load_data temp_load;
		temp_load = load.second;

		// Add to the load list
		this->model_loads.add_load(temp_load.node_id, temp_load.load_loc, temp_load.load_start_time,
			temp_load.load_end_time, temp_load.load_value, temp_load.load_angle);
	}

	// Add to model point loads
	for (auto& ptmass : model_ptmass.ptmassMap)
	{
		// Create a temporart point mass
		nodepointmass_data temp_ptmass;
		temp_ptmass = ptmass.second;

		// Add to the point mass list
		this->model_ptmass.add_pointmass(temp_ptmass.node_id, temp_ptmass.ptmass_loc, temp_ptmass.ptmass_defl,
			temp_ptmass.ptmass_x, temp_ptmass.ptmass_y, temp_ptmass.is_offset);
	}

	// Add to model initial conditions
	for (auto& inlcond_m : model_inlcond.inlcondMap)
	{
		// Create a temporary initial condition data
		nodeinl_condition_data temp_inlcond;
		temp_inlcond = inlcond_m.second;


		// Add to the initial condition list
		this->model_inlcond.add_inlcondition(temp_inlcond.node_id, temp_inlcond.inlcond_loc, temp_inlcond.inl_displacement_x,
			temp_inlcond.inl_displacement_y, temp_inlcond.inl_velocity_x, temp_inlcond.inl_velocity_y);
	}


	// Geometry is loaded
	is_geometry_set = true;

	// Set the boundary of the geometry
	std::pair<glm::vec2, glm::vec2> result = findMinMaxXY(model_nodes.nodeMap);
	this->geom_param.min_b = result.first;
	this->geom_param.max_b = result.second;
	this->geom_param.geom_bound = geom_param.max_b - geom_param.min_b;

	// Set the center of the geometry
	this->geom_param.center = findGeometricCenter(model_nodes.nodeMap);

	// Set the geometry
	update_model_matrix();
	update_model_zoomfit();

	// Set the geometry buffers
	this->model_nodes.set_buffer();
	this->model_lineelements.set_buffer();
	this->model_constarints.set_buffer();
	this->model_loads.set_buffer();
	this->model_ptmass.set_buffer();
	this->model_inlcond.set_buffer();

	// Initialize the modal analysis result nodes and lines
	modal_results.clear_data();
	modal_result_nodes.init(&geom_param);
	modal_result_lineelements.init(&geom_param);

	// Clear the modal results in the window
	sol_modal_window->init();
	sol_pulse_window->init();
}

std::pair<glm::vec2, glm::vec2> geom_store::findMinMaxXY(const std::unordered_map<int, node_store>& model_nodes)
{
	// Initialize min and max values to first node in map
	glm::vec2 firstNode = model_nodes.begin()->second.node_pt;
	glm::vec2 minXY = glm::vec2(firstNode.x, firstNode.y);
	glm::vec2 maxXY = minXY;

	// Loop through all nodes in map and update min and max values
	for (auto it = model_nodes.begin(); it != model_nodes.end(); ++it)
	{
		const auto& node = it->second.node_pt;
		if (node.x < minXY.x)
		{
			minXY.x = node.x;
		}
		if (node.y < minXY.y)
		{
			minXY.y = node.y;
		}
		if (node.x > maxXY.x)
		{
			maxXY.x = node.x;
		}
		if (node.y > maxXY.y)
		{
			maxXY.y = node.y;
		}
	}

	// Return pair of min and max values
	return { minXY, maxXY };
}

glm::vec2 geom_store::findGeometricCenter(const std::unordered_map<int, node_store>& model_nodes)
{
	// Function returns the geometric center of the nodes
		// Initialize the sum with zero
	glm::vec2 sum(0);

	// Sum the points
	for (auto it = model_nodes.begin(); it != model_nodes.end(); ++it)
	{
		sum += it->second.node_pt;
	}
	return sum / static_cast<float>(model_nodes.size());
}

void geom_store::update_delete_material(int& del_material_id)
{
	// Update delete material
	bool is_del_material_found = false;

	// Delete the material
	for (int i = 0; i < model_lineelements.elementlineMap.size(); i++)
	{
		if (model_lineelements.elementlineMap[i].material_id == del_material_id)
		{
			// Delete material is removed and the material ID of that element to 0
			model_lineelements.elementlineMap[i].material_id = 0;
			is_del_material_found = true;
		}
	}

	// Update the material ID label
	if (is_del_material_found == true)
	{
		model_lineelements.update_material_id_labels();
	}
}
