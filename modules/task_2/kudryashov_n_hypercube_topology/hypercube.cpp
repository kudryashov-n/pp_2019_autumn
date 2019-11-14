// Copyright 2019 Kudryashov Nikita
#include <mpi.h>
#include "../../../modules/task_2/kudryashov_n_hypercube_topology/hypercube.h"

MPI_Comm createHcube(unsigned int dim) {
    if (dim == 0) {
        throw "Dimensions of hypercube can't be equal to zero.";
    }
    int* dimSize = new int[dim];
    int* periods = new int[dim];
    for (unsigned int i = 0; i < dim; i++) {
        dimSize[i] = 2;
        periods[i] = 0;
    }

    MPI_Comm hypercubeComm;
    MPI_Cart_create(MPI_COMM_WORLD, dim, dimSize, periods, 1, &hypercubeComm);

    delete[] dimSize;
    delete[] periods;

    return hypercubeComm;
}

