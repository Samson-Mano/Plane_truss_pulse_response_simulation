#include "modal_analysis_result_store.h"

modal_analysis_result_store::modal_analysis_result_store()
{
	// Empty constructor
}

modal_analysis_result_store::~modal_analysis_result_store()
{
	// Empty destructor
}

void modal_analysis_result_store::clear_data()
{
	// Clear the eigen values and eigen vectors
	number_of_modes = 0;
	nodeid_map.clear();
	eigen_values.clear();
	eigen_vectors.clear();
	eigen_vectors_reduced.clear();
	mode_result_str.clear();
}

void modal_analysis_result_store::add_node_map(std::unordered_map<int, int>& nodeid_map)
{
	// Add the node map for the results
	this->nodeid_map = nodeid_map;
}

void modal_analysis_result_store::add_eigen_data(int mode_number, double eigen_value, std::vector<double> eigen_vectors, std::vector<double> eigen_vectors_reduced)
{
	// insert the eigen value
	eigen_values.insert({ mode_number,eigen_value });

	// insert the eigen vectors of this particular eigen value
	this->eigen_vectors.insert({ mode_number,eigen_vectors });

	this->eigen_vectors_reduced.insert({ mode_number, eigen_vectors_reduced });

	// Iterate the number of modes
	number_of_modes++;
}
