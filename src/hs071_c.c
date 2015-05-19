#include "hs071_c.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* This is an example how user_data can be used. */
struct MyUserData
{
	Number g_offset[2]; /* This is an offset for the constraints.  */
};


static const size_t NUM_VARIABLES = 4;
static const size_t NUM_BINARY_VARIABLES = 2;
static const size_t NUM_CONSTRAINTS = 2;
static const Index BINARY_INDICES[NUM_BINARY_VARIABLES] = {0,1};
static Number g_L[NUM_CONSTRAINTS] = {0.0, 0.0};
static Number g_U[NUM_CONSTRAINTS] = {2e19, 40};
/* Number of nonzeros in the Jacobian of the constraints */
static const Index nele_jac = 8;
/* Number of nonzeros in the Hessian of the Lagrangian (lower or
	 upper triangual part only) */
static const Index nele_hess = 10;
/* indexing style for matrices */
static const Index index_style = 0; /* C-style; start counting of rows and column
					 indices at 0 */


/* Function Implementations */
Bool eval_f(Index n, Number* x, Bool new_x,
	Number* obj_value, UserDataPtr user_data)
{
	assert(n == 4);

	*obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

	return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
	Number* grad_f, UserDataPtr user_data)
{
	assert(n == 4);

	grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad_f[1] = x[0] * x[3];
	grad_f[2] = x[0] * x[3] + 1;
	grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

	return TRUE;
}


Bool eval_g(Index n, Number* x, Bool new_x,
	Index m, Number* g, UserDataPtr user_data)
{
	struct MyUserData* my_data = (struct MyUserData*)user_data;

	assert(n == 4);
	assert(m == 2);

	g[0] = x[0] * x[1] * x[2] * x[3] + my_data->g_offset[0];
	g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + my_data->g_offset[1];

	return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
	Index m, Index nele_jac,
	Index *iRow, Index *jCol, Number *values,
	UserDataPtr user_data)
{
	if (values == NULL) {
	/* return the structure of the jacobian */

	/* this particular jacobian is dense */
		iRow[0] = 0;
		jCol[0] = 0;
		iRow[1] = 0;
		jCol[1] = 1;
		iRow[2] = 0;
		jCol[2] = 2;
		iRow[3] = 0;
		jCol[3] = 3;
		iRow[4] = 1;
		jCol[4] = 0;
		iRow[5] = 1;
		jCol[5] = 1;
		iRow[6] = 1;
		jCol[6] = 2;
		iRow[7] = 1;
		jCol[7] = 3;
	}
	else {
	/* return the values of the jacobian of the constraints */

	values[0] = x[1]*x[2]*x[3]; /* 0,0 */
	values[1] = x[0]*x[2]*x[3]; /* 0,1 */
	values[2] = x[0]*x[1]*x[3]; /* 0,2 */
	values[3] = x[0]*x[1]*x[2]; /* 0,3 */

	values[4] = 2*x[0];         /* 1,0 */
	values[5] = 2*x[1];         /* 1,1 */
	values[6] = 2*x[2];         /* 1,2 */
	values[7] = 2*x[3];         /* 1,3 */
	}

	return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
	Index m, Number *lambda, Bool new_lambda,
	Index nele_hess, Index *iRow, Index *jCol,
	Number *values, UserDataPtr user_data)
{
Index idx = 0; /* nonzero element counter */
Index row = 0; /* row counter for loop */
Index col = 0; /* col counter for loop */
	if (values == NULL) {
	/* return the structure. This is a symmetric matrix, fill the lower left
	 * triangle only. */

	/* the hessian for this problem is actually dense */
		idx=0;
		for (row = 0; row < 4; row++) {
			for (col = 0; col <= row; col++) {
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
			}
		}

		assert(idx == nele_hess);
	}
	else {
	/* return the values. This is a symmetric matrix, fill the lower left
	 * triangle only */

	/* fill the objective portion */
	values[0] = obj_factor * (2*x[3]);               /* 0,0 */

	values[1] = obj_factor * (x[3]);                 /* 1,0 */
	values[2] = 0;                                   /* 1,1 */

	values[3] = obj_factor * (x[3]);                 /* 2,0 */
	values[4] = 0;                                   /* 2,1 */
	values[5] = 0;                                   /* 2,2 */

	values[6] = obj_factor * (2*x[0] + x[1] + x[2]); /* 3,0 */
	values[7] = obj_factor * (x[0]);                 /* 3,1 */
	values[8] = obj_factor * (x[0]);                 /* 3,2 */
	values[9] = 0;                                   /* 3,3 */


	/* add the portion for the first constraint */
	values[1] += lambda[0] * (x[2] * x[3]);          /* 1,0 */

	values[3] += lambda[0] * (x[1] * x[3]);          /* 2,0 */
	values[4] += lambda[0] * (x[0] * x[3]);          /* 2,1 */

	values[6] += lambda[0] * (x[1] * x[2]);          /* 3,0 */
	values[7] += lambda[0] * (x[0] * x[2]);          /* 3,1 */
	values[8] += lambda[0] * (x[0] * x[1]);          /* 3,2 */

	/* add the portion for the second constraint */
	values[0] += lambda[1] * 2;                      /* 0,0 */

	values[2] += lambda[1] * 2;                      /* 1,1 */

	values[5] += lambda[1] * 2;                      /* 2,2 */

	values[9] += lambda[1] * 2;                      /* 3,3 */
	}

	return TRUE;
}

Bool intermediate_cb(Index alg_mod, Index iter_count, Number obj_value,
	Number inf_pr, Number inf_du, Number mu, Number d_norm,
	Number regularization_size, Number alpha_du,
	Number alpha_pr, Index ls_trials, UserDataPtr user_data)
{
	printf("Testing intermediate callback in iteration %d\n", iter_count);
	if (inf_pr < 1e-4) return FALSE;

	return TRUE;
}

Number *create_binary_input_combinations()
{
	const size_t NUM_COMBINATIONS = get_num_combinations();
    Number *constraint_inputs = (Number*) malloc(NUM_COMBINATIONS * NUM_BINARY_VARIABLES * sizeof(Number));

    for (int i = 0; i < NUM_COMBINATIONS; ++i)
    {
        for (int j = 0; j < NUM_BINARY_VARIABLES; ++j)
        {
        	size_t k = i*NUM_BINARY_VARIABLES + j;
            constraint_inputs[k] = (Number)((i >> j) & (0x1));
        }

    }
    return constraint_inputs;
}


Number solve(Number *fixed_binary_values)
{
	// Number* x_L = NULL;                  /* lower bounds on x */
	// Number* x_U = NULL;                  /* upper bounds on x */

	IpoptProblem nlp = NULL;             /* IpoptProblem */
	enum ApplicationReturnStatus status; /* Solve return code */
	Number* x = NULL;                    /* starting point and solution vector */
	Number* mult_g = NULL;               /* constraint multipliers at the solution */
	Number* mult_x_L = NULL;             /* lower bound multipliers at the solution */
	Number* mult_x_U = NULL;             /* lower bound multipliers at the solution */
	Number* x_L = NULL;
	Number* x_U = NULL;
	Number obj;                          /* objective value */


	/* our user data for the function evalutions. */
	struct MyUserData user_data;

	// /* set the number of variables and allocate space for the bounds */
	// x_L = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);
	// x_U = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);
	// /* allocate space for the initial point and set the values */
	x = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);
	x[0] = 1.0;
	x[1] = 5.0;
	x[2] = 5.0;
	x[3] = 1.0;
	/* allocate space to store the bound multipliers at the solution */
	mult_g = (Number*)malloc(sizeof(Number)*NUM_CONSTRAINTS);
	mult_x_L = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);
	mult_x_U = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);

	x_L = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);
	x_U = (Number*)malloc(sizeof(Number)*NUM_VARIABLES);


	/* set the values for the variable bounds */
	for (int i=0; i<NUM_VARIABLES; i++)
	{
		x_L[i] = 1.0;
		x_U[i] = 5.0;
	}
	for (int i = 0; i < NUM_BINARY_VARIABLES; ++i)
	{
		Index j = BINARY_INDICES[i];
		x_L[j] = x_U[j] = fixed_binary_values[i];
	}

	/* create the IpoptProblem */
	nlp = CreateIpoptProblem(NUM_VARIABLES, x_L, x_U, NUM_CONSTRAINTS, g_L, g_U, nele_jac, nele_hess, index_style, &eval_f, &eval_g, &eval_grad_f, &eval_jac_g, &eval_h);
	status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);
	FreeIpoptProblem(nlp);
	free(x);
	free(mult_g);
	free(mult_x_L);
	free(mult_x_U);
	if(status != Solve_Succeeded)
	{
		uint64_t bits = 0x7ff0000000000000;
		obj = *(Number *)&bits;
	}
	return obj;
}

size_t get_num_combinations()
{
	return  (size_t)pow(2, NUM_BINARY_VARIABLES);
}

size_t get_num_binary_variables()
{
	return NUM_BINARY_VARIABLES;
}

size_t get_num_vars()
{
	return NUM_VARIABLES;
}

Index* get_binary_indices()
{
	static const size_t SIZE = NUM_BINARY_VARIABLES*sizeof(Index);
	Index *binary_indices = (Index*)malloc( SIZE );
	memcpy(binary_indices, BINARY_INDICES, SIZE);
	return binary_indices;
}
