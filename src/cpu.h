#include <thread>
#include <vector>

size_t num_logical_cores();
size_t num_physical_cores();
void distribute_threads(std::vector<std::thread>& threads);