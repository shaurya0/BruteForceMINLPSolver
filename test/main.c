#include "hs071_c.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>


void brute_force_minlp_solve(Number *obj_value, Number *optimal_binary_input, size_t NUM_BINARY_VARIABLES)
{
    Index *binary_indices             = get_binary_indices();
    const size_t NUM_VARIABLES        = get_num_vars();
    const size_t NUM_COMBINATIONS     = get_num_combinations();

    Number *fixed_binary_values = (Number*)malloc(NUM_BINARY_VARIABLES*sizeof(Number));
    Number *obj_values          = (Number*)malloc(NUM_COMBINATIONS*sizeof(Number));
    Number *constraint_inputs   = create_binary_input_combinations();

    Number minimum_value = 2e19;
    Index min_index      = 0;


    for (int i = 0, k = 0; i < NUM_COMBINATIONS*NUM_BINARY_VARIABLES; i+=NUM_BINARY_VARIABLES, k++)
    {

        memcpy(fixed_binary_values, &constraint_inputs[i], NUM_BINARY_VARIABLES*sizeof(Number));
        obj_values[k] = solve(fixed_binary_values);
        if (obj_values[k] < minimum_value)
        {
            minimum_value = obj_values[k];
            min_index = k;
        }
    }
    *obj_value = minimum_value;
    memcpy(optimal_binary_input, &constraint_inputs[min_index*NUM_BINARY_VARIABLES], NUM_BINARY_VARIABLES*sizeof(Number));

    free(constraint_inputs);
    free(obj_values);
    free(fixed_binary_values);
    free(binary_indices);
}


int main(int argc, char const *argv[])
{
    const size_t NUM_BINARY_VARIABLES = get_num_binary_variables();
    Number *minimum_value             = (Number*)malloc(sizeof(Number));
    Number *optimal_binary_input      = (Number*)malloc(NUM_BINARY_VARIABLES*sizeof(Number));

    brute_force_minlp_solve(minimum_value, optimal_binary_input, NUM_BINARY_VARIABLES);
    free(minimum_value);
    free(optimal_binary_input);

    return 0;
}
