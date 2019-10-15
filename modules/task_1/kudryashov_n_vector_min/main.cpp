// Copyright 2019 Kudryashov Nikita
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./vector_min.h"

TEST(SequentialVectorMin, low_positive) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> low(4);
    int res;
    for (int i = 0; i < 4; i++) {
        low[i] = i;
    }

    res = getVectorMinSequential(low);

    if (rank == 0) {
        ASSERT_EQ(res, 0);
    }
}

TEST(SequentialVectorMin, all_positive_except_one_negative) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> one_pos(100);
    int res;
    for (int i = 100; i < 199; i++) {
        one_pos[i - 100] = i * i;
    }
    one_pos[99] = - 1;
    res = getVectorMinSequential(one_pos);

    if (rank == 0) {
        ASSERT_EQ(res, -1);
    }
}

TEST(ParallelVectorMin, low_amount) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> low(5);

    for (int i = 0; i < 5; i++) {
        low[i] = i;
    }

    int paral_res = getVectorMinParallel(low);
    int seq_res = getVectorMinSequential(low);

    if (rank == 0) {
        ASSERT_EQ(paral_res, seq_res);
    }
}

TEST(ParallelVectorMin, high_amount_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> high;

    high = getRandomVector(100);

    int paral_res = getVectorMinParallel(high);
    int seq_res = getVectorMinSequential(high);

    if (rank == 0) {
        ASSERT_EQ(paral_res, seq_res);
    }
}

TEST(ParallelVectorMin, high_amount_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> high;

    high = getRandomVector(500);

    int paral_res = getVectorMinParallel(high);
    int seq_res = getVectorMinSequential(high);

    if (rank == 0) {
        ASSERT_EQ(paral_res, seq_res);
    }
}

TEST(ParallelVectorMin, high_amount_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> high;

    high = getRandomVector(1000);

    int paral_res = getVectorMinParallel(high);
    int seq_res = getVectorMinSequential(high);

    if (rank == 0) {
        ASSERT_EQ(paral_res, seq_res);
    }
}

TEST(ParallelVectorMin, all_positive_decreasing) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> decr(100);

    for (int i = 0; i < 100; i++) {
        decr[i] = 100 - i;
    }

    int paral_res = getVectorMinParallel(decr);
    int seq_res = getVectorMinSequential(decr);

    if (rank == 0) {
        ASSERT_EQ(seq_res, paral_res);
    }
}

TEST(ParallelVectorMin, all_negative_increasing) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> incr(100);

    for (int i = 0; i < 100; i++) {
        incr[i] = -100 + i;
    }

    int paral_res = getVectorMinParallel(incr);
    int seq_res = getVectorMinSequential(incr);

    if (rank == 0) {
        ASSERT_EQ(seq_res, paral_res);
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
