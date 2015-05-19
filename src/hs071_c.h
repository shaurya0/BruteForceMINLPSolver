#include "coin/IpStdCInterface.h"
#include "sys/types.h"
#include <stdint.h>

// Custom functions
Number solve(Number *fixed_binary_values);
Number *create_binary_input_combinations();
size_t get_num_combinations();
size_t get_num_binary_variables();
size_t get_num_vars();
Index* get_binary_indices();

