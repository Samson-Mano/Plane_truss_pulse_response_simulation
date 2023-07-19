#pragma once
#include <iostream>
#include <unordered_map>
#include <vector>


class pulse_analysis_result_store
{
public:
	int time_step_count = 0;
	double time_interval = 0.0;
	double total_simulation_time = 0.0;

	pulse_analysis_result_store();
	~pulse_analysis_result_store();
	void set_analysis_setting(const int& time_step_count, const double& time_interval, const double& total_simulation_time);
	void clear_results();

private:

};
