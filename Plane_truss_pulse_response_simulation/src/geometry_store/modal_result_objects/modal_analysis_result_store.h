#pragma once
#include <iostream>
#include <unordered_map>
#include <vector>

class modal_analysis_result_store
{
public:
	int number_of_modes = 0;
	std::unordered_map<int, int> nodeid_map;
	std::unordered_map<int, double> eigen_values;
	std::unordered_map<int, std::vector<double>> eigen_vectors;
	std::unordered_map<int, std::vector<double>> eigen_vectors_reduced;
	std::vector<std::string> mode_result_str;

	modal_analysis_result_store();
	~modal_analysis_result_store();
	void clear_data();
	void add_node_map(std::unordered_map<int, int>& nodeid_map);
	void add_eigen_data(int mode_number, double eigen_value, std::vector<double> eigen_vectors, std::vector<double> eigen_vectors_reduced);

private:
};
