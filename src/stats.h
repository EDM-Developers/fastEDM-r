#include <vector>

double median(std::vector<double> u);
std::vector<int> rank(const std::vector<double>& v_temp);
std::vector<double> remove_value(const std::vector<double>& vec, double target);
double correlation(const std::vector<double>& y1, const std::vector<double>& y2);
double mean_absolute_error(const std::vector<double>& y1, const std::vector<double>& y2);
double default_missing_distance(const std::vector<double>& x);
double default_dt_weight(const std::vector<double>& t, const std::vector<double>& x);