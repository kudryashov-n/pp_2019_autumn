// Copyright 2019 Kudryashov Nikita
#ifndef MODULES_TASK_3_KUDRYASHOV_N_FOX_ALGORITHM_FOX_ALGORITHM_H_
#define MODULES_TASK_3_KUDRYASHOV_N_FOX_ALGORITHM_FOX_ALGORITHM_H_

#include <mpi.h>

double* fox_mult(double* a, unsigned int a_size, double* b, unsigned int b_size);
double* subtask_matr_mult(double* a_block, unsigned int a_block_size,
    double* b_block, unsigned int b_block_size);  // Matrix multiplication.
bool is_square(int size);

#endif  // MODULES_TASK_3_KUDRYASHOV_N_FOX_ALGORITHM_FOX_ALGORITHM_H_
