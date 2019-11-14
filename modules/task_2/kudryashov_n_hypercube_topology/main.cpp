// Copyright 2019 Kudryashov Nikita
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include "./hypercube.h"

TEST(HCubeTopology, can_create_hcube_dim_1) {
    MPI_Comm Hcube;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ASSERT_NO_THROW(Hcube = createHcube(1));
}

TEST(HCubeTopology, can_create_hcube_dim_2) {
    MPI_Comm Hcube;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ASSERT_NO_THROW(Hcube = createHcube(2));
}

TEST(HCubeTopology, can_create_hcube_dim_3) {
    MPI_Comm Hcube;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ASSERT_NO_THROW(Hcube = createHcube(3));
}

TEST(HCubeTopology, can_create_hcube_dim_4) {
    MPI_Comm Hcube;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ASSERT_NO_THROW(Hcube = createHcube(4));
}

TEST(HCubeTopology, throw_with_zero_dim) {
    MPI_Comm Hcube;
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ASSERT_ANY_THROW(Hcube = createHcube(0));
}

TEST(HCubeTopology, topology_check) {
    MPI_Comm Hcube;
    int status;

    Hcube = createHcube(4);

    MPI_Topo_test(Hcube, &status);
    ASSERT_EQ(status, MPI_CART);
}

TEST(HCubeTopology, dims_4_param_dims_check) {
    MPI_Comm Hcube;
    int rank, dims;

    Hcube = createHcube(4);

    MPI_Comm_rank(Hcube, &rank);

    if (rank == 0) {
        MPI_Cartdim_get(Hcube, &dims);
        ASSERT_EQ(dims, 4);
    }
}

TEST(HCubeTopology, dims_4_param_periods_check) {
    MPI_Comm Hcube;
    int rank, dims;

    Hcube = createHcube(4);

    MPI_Comm_rank(Hcube, &rank);

    if (rank == 0) {
        MPI_Cartdim_get(Hcube, &dims);

        int* dim = new int[dims];
        int* period = new int[dims];
        int* coords = new int[dims];
        bool check = true;

        MPI_Cart_get(Hcube, dims, dim, period, coords);

        for (int i = 0; i < dims; i++) {
            if (period[i] != 0) {
                check = false;
                break;
            }
        }
    }
}

TEST(HCubeTopology, dims_4_param_coords_check) {
    MPI_Comm Hcube;
    int rank, dims;

    Hcube = createHcube(4);

    MPI_Comm_rank(Hcube, &rank);

    if (rank == 0) {
        MPI_Cartdim_get(Hcube, &dims);

        int* dim = new int[dims];
        int* period = new int[dims];
        int* coords = new int[dims];
        int* rank_coords = new int[dims];
        bool check = true;

        MPI_Cart_get(Hcube, dims, dim, period, coords);
        MPI_Cart_coords(Hcube, 0, dims, rank_coords);

        for (int i = 0; i < dims; i++) {
            if (rank_coords[i] != coords[i]) {
                check = false;
                break;
            }
        }
        ASSERT_EQ(check, true);
    }
}

TEST(HCubeTopology, dims_4_param_dim_array_check) {
    MPI_Comm Hcube;
    int rank, dims;

    Hcube = createHcube(4);

    MPI_Comm_rank(Hcube, &rank);

    if (rank == 0) {
        MPI_Cartdim_get(Hcube, &dims);

        int* dim = new int[dims];
        int* period = new int[dims];
        int* coords = new int[dims];
        bool check = true;

        MPI_Cart_get(Hcube, dims, dim, period, coords);

        for (int i = 0; i < dims; i++) {
            if (dim[i] != 2) {
                check = false;
                break;
            }
        }
        ASSERT_EQ(check, true);
    }
}

TEST(HCubeTopology, can_send_msg) {
    MPI_Comm Hcube;
    int rank;

    Hcube = createHcube(4);

    MPI_Comm_rank(Hcube, &rank);

    if (rank == 0) {
        int msg = 15;

        MPI_Send(&msg, 1, MPI_INT, 5, 1, Hcube);
    }
    if (rank == 5) {
        int msg, ret;
        MPI_Status status;

        ret = MPI_Recv(&msg, 1, MPI_INT, 0, 1, Hcube, &status);

        ASSERT_TRUE(ret == MPI_SUCCESS && msg == 15);
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
