// Copyright 2019 Kudryashov Nikita

#include "../../../modules/task_3/kudryashov_n_fox_algorithm/fox_algorithm.h"
#include <cmath>

bool is_square(int size) {
    int i = 1;
    while (i <= size) {
        if (i * i == size) {
            return true;
        } else {
            i++;
        }
    }
    return false;
}

double* subtask_matr_mult(double* a_block, unsigned int a_block_size, double* b_block, unsigned int b_block_size) {
    double* c = new double[a_block_size * a_block_size];
    for (unsigned int i = 0; i < a_block_size; i++) {
        for (unsigned int j = 0; j < a_block_size; j++) {
            c[i * a_block_size + j] = 0.0;
            for (unsigned int k = 0; k < a_block_size; k++) {
                c[i * a_block_size + j] += a_block[i * a_block_size + k] * b_block[k * a_block_size + j];
            }
        }
    }
    return c;
}

double* fox_mult(double* a, unsigned int a_size, double* b, unsigned int b_size) {
    if (a_size != b_size) {
        throw "Matrix must have equal size.";
    }

    bool matr_size_changed = false;
    int save_size;
    int size, rank, size_inside;
    MPI_Status status;
    MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size == 1) {
        double* c;
        return c = subtask_matr_mult(a, a_size, b, b_size);
    }

    if (!is_square(size) && static_cast<unsigned int>(size) < a_size * a_size) {
        unsigned int i = 1;
        while (i*i < static_cast<unsigned int>(size)) {
            i++;
        }
        i--;
        size_inside = i * i;
    } else {
        if (static_cast<unsigned int>(size) < a_size * a_size) {
            size_inside = size;
        } else {
            size_inside = a_size * a_size;
        }
    }

    save_size = a_size;
    double* temp_a = new double[save_size * save_size];
    double* temp_b = new double[save_size * save_size];

    if (a_size % ((unsigned int)sqrt(size_inside)) != 0) {
    // Creating matrix of new size and surrounding by zeros'.
        matr_size_changed = true;
        a_size = a_size + (unsigned int)sqrt(size_inside) - a_size % ((unsigned int)sqrt(size_inside));
        // Saving A.
        for (int i = 0; i < save_size; i++) {
            for (int j = 0; j < save_size; j++) {
                temp_a[i*save_size+j] = a[i*save_size+j];
            }
        }
        delete[] a;
        a = new double[a_size * a_size];
        for (int i = 0; i < save_size; i++) {
            for (int j = 0; j < save_size; j++) {
                a[i*a_size+j] = temp_a[i*save_size+j];
            }
        }
        // Left: not full.
        for (int i = save_size; i < static_cast<int>(a_size); i++) {
            for (int j = 0; j < save_size; j++) {
                a[i*a_size+j] = 0.0;
            }
        }
        // Right: full.
        for (unsigned int i = 0; i < a_size; i++) {
            for (int j = save_size; j < static_cast<int>(a_size); j++) {
                a[i*a_size+j] = 0.0;
            }
        }
        // Finished creating new A.

        // Saving B.
        for (int i = 0; i < save_size; i++) {
            for (int j = 0; j < save_size; j++) {
                temp_b[i*save_size+j] = b[i*save_size+j];
            }
        }
        delete[] b;
        b = new double[a_size * a_size];
        for (int i = 0; i < save_size; i++) {
            for (int j = 0; j < save_size; j++) {
                b[i*a_size+j] = temp_b[i*save_size+j];
            }
        }
        // Left: not full.
        for (int i = save_size; i < static_cast<int>(a_size); i++) {
            for (int j = 0; j < save_size; j++) {
                b[i*a_size+j] = 0.0;
            }
        }
        // Right: full.
        for (int i = 0; i < static_cast<int>(a_size); i++) {
            for (int j = save_size; j < static_cast<int>(a_size); j++) {
                b[i*a_size+j] = 0.0;
            }
        }
            // Finished creating new B.
    }

    int block_size = a_size / static_cast<int>(sqrt(size_inside));

    int color, key = rank;
    MPI_Comm inside;
    if (rank < size_inside) {
        color = 1;
    } else {
        color = 0;
    }
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &inside);
    MPI_Comm_rank(inside, &rank);

    if (color == 1) {
        double* a_block = new double[block_size * block_size];
        double* b_block = new double[block_size * block_size];
        double* c;

        // Initializing process
        int blocks_per_row = a_size / block_size;
        int blocks_per_column = a_size / block_size;
        for (int i = 0; i < block_size; i++) {
            for (int j = 0; j < block_size; j++) {
                a_block[block_size * i + j] = a[rank / blocks_per_row * block_size * a_size +
                    ((rank + rank / blocks_per_row) % blocks_per_row) * block_size + i * a_size + j];
                b_block[block_size * i + j] = b[((rank +(rank % blocks_per_column) *
                    blocks_per_column) % (blocks_per_column * blocks_per_column))/ blocks_per_row *
                    block_size * a_size + rank % blocks_per_row * block_size + i * a_size + j];
                // b_block form operation is the same as a_block form operation, but matrix b is transposed by blocks
            }
        }
        double* res = new double[a_size * a_size];
        for (unsigned int i = 0; i < a_size; i++) {
            for (unsigned int j = 0; j < a_size; j++) {
                res[i * a_size + j] = 0.0;
            }
        }
        // Initialization over

        // Subtask mult and move
        for (int k = 0; k < blocks_per_row; k++) {
            // Subtask mult and collecting the result
            MPI_Barrier(inside);
            c = subtask_matr_mult(a_block, block_size, b_block, block_size);
            if (rank == 0) {
                double* recv_arr = new double[block_size * block_size];
                for (int i = 0; i < block_size; i++) {
                    for (int j = 0; j < block_size; j++) {
                        res[i * a_size + j] += c[i * block_size +j];
                    }
                }

                int coord_i;
                int coord_j;
                for (int i = 1; i < size_inside; i++) {
                    coord_i = i / blocks_per_row;
                    coord_j = i % blocks_per_column;
                    MPI_Recv(recv_arr, block_size * block_size, MPI_DOUBLE, i, 0, inside, &status);
                    for (int d1 = 0; d1 < block_size; d1++) {
                        for (int d2 = 0; d2 < block_size; d2++) {
                            res[block_size * a_size * coord_i + block_size*coord_j + d1*a_size + d2]
                                += recv_arr[d1 * block_size + d2];
                        }
                    }
                }
            } else {
                MPI_Send(c, block_size * block_size, MPI_DOUBLE, 0, 0, inside);
            }

            // Move
            // A blocks move
            int additive;

            if (rank - 1  < blocks_per_row * (rank / blocks_per_row)) {  // Getting lesser than first in a row
                additive = blocks_per_row;
            } else {
                additive = 0;
            }

            MPI_Barrier(inside);
            MPI_Isend(a_block, block_size * block_size, MPI_DOUBLE, rank - 1 + additive, 0, inside, &request);
            MPI_Barrier(inside);

            if (rank + 1  > blocks_per_row * (rank / blocks_per_row)
                    + blocks_per_row - 1) {  // Getting bigger than last in a row
                additive = blocks_per_row;
            } else {
                additive = 0;
            }

            MPI_Barrier(inside);
            MPI_Irecv(a_block, block_size * block_size, MPI_DOUBLE, rank + 1 - additive, 0, inside, &request);
            MPI_Barrier(inside);

            // B blocks move
            int dest;

            // Mod for negative number
            if (rank - blocks_per_column >= 0) {
                dest = (rank - blocks_per_column) % (blocks_per_row * blocks_per_column);
            } else {
                dest = blocks_per_column * blocks_per_column + (rank - blocks_per_column);
            }

            MPI_Barrier(inside);
            MPI_Isend(b_block, block_size * block_size, MPI_DOUBLE, dest, 0, inside, &request);
            MPI_Barrier(inside);

            MPI_Barrier(inside);
            MPI_Irecv(b_block, block_size * block_size, MPI_DOUBLE,
                (rank + blocks_per_column) % (blocks_per_row * blocks_per_column), 0, inside, &request);
            MPI_Barrier(inside);
        }

        delete[] a_block;
        delete[] b_block;
        MPI_Barrier(inside);
        MPI_Comm_free(&inside);

        if (rank == 0) {
            if (matr_size_changed) {
                double* res_new = new double[save_size * save_size];
                // Restoring result.
                for (int i = 0; i < save_size; i++) {
                    for (int j = 0; j < save_size; j++) {
                        res_new[i*save_size+j] = res[i*a_size+j];
                    }
                }
                delete[] res;

                // Restoring A.
                delete[] a;
                a = new double[save_size * save_size];
                for (int i = 0; i < save_size; i++) {
                    for (int j = 0; j < save_size; j++) {
                        a[i*save_size+j] = temp_a[i*save_size+j];
                    }
                }

                // Restoring B.
                delete[] b;
                b = new double[save_size * save_size];
                for (int i = 0; i < save_size; i++) {
                    for (int j = 0; j < save_size; j++) {
                        b[i*save_size+j] = temp_b[i*save_size+j];
                    }
                }

                delete[] temp_a;
                delete[] temp_b;
                return res_new;
            } else {
                delete[] temp_a;
                delete[] temp_b;
                return res;
            }
        } else  {
            delete[] temp_a;
            delete[] temp_b;
            delete[] res;
            return a;
        }
    }
    delete[] temp_a;
    delete[] temp_b;
    MPI_Barrier(inside);
    MPI_Comm_free(&inside);
    return a;
}
