#include <vector>

double median(std::vector<double> u);
std::vector<int> rank(const std::vector<double>& v_temp);
template<typename T>
std::vector<T> remove_value(const std::vector<T>& vec, T target);
double correlation(const std::vector<double>& y1, const std::vector<double>& y2);
double mean_absolute_error(const std::vector<double>& y1, const std::vector<double>& y2);
double standard_deviation(const std::vector<double>& vec);
double default_missing_distance(const std::vector<double>& x);
