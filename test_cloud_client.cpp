//
// Created by Ruochen WANG on 16/4/2020.
//

#include "test_cloud_client.h"
#include "secret_sharing.h"
#include "utils.h"
#include "decision_tree.h"
#include "network.h"
#include <chrono>
#include <cmath>
#include <iostream>

extern gmp_randclass gmp_prn;

const int param_nd[5][2] = {{9,  8},
                            {13, 3},
                            {13, 13},
                            {15, 4},
                            {57, 17}};

auto start = std::chrono::steady_clock::now(), end = std::chrono::steady_clock::now();
#define CLOCK_START {start = std::chrono::steady_clock::now();}
#define CLOCK_END {end = std::chrono::steady_clock::now();}
#define ELAPSED std::chrono::duration<double, std::nano>(end - start).count()

void test_by_phases(std::vector<int> phases, int num_trial) {
    for (auto phase : phases) {
        if (phase == 1) {
            phase1(num_trial);
        } else if (phase == 2) {
            phase2(num_trial);
        } else if (phase == 3) {
            phase3(num_trial);
        } else if (phase == 4) {
            phase4(num_trial);
        } else {
            std::cout << "Phase No. out of range." << std::endl;
        }
    }
}

void phase1(int num_trial) {
    // node selection
    // TODO
    int n, m;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        m = pow(2, param_nd[i][1]) - 1;

        ss_tuple_mz sel_mat(m, n);

        sel_mat.encrypt();

        ss_tuple_mz feature(n, 1);
        feature.encrypt();

        tri_tuple_mz tri(m, n, 1);
        tri.encrypt();

        ss_tuple_mz buf_E(m, n);
        ss_tuple_mz buf_f(n, 1);

        ss_tuple_mz xsigma(m, 1);

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            secure_multiplication(sel_mat.share, feature.share,
                                  tri.share,
                                  buf_E, buf_f,
                                  xsigma.share,
                                  0);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }

        printf("secure node selection (n=%d, d=%d, m=%d): %f ns\n", n, param_nd[i][1], m,
               time_total / num_trial / 2); // count only one party
    }
}

void phase2(int num_trial) {
    // node evaluation
//    TODO
    int n, m;
// secure node evaluation
    triplet_b tri_b;
    triplet_z tri_z;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        m = pow(2, param_nd[i][1]) - 1;

        // ugly staff for preparation
        mpz_class (*x)[2] = new mpz_class[m][2], (*y)[2] = new mpz_class[m][2];
        for (int j = 0; j < m; ++j) {
            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), x[j]);
            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), y[j]);
        }

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            for (int k = 0; k < m; ++k)
                secure_node_evaluation(x[k], y[k], tri_z, tri_b);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }

        printf("secure node evaluation (n=%d, d=%d, m=%d): %f ns\n", n, param_nd[i][1], m,
               time_total / num_trial / 2); // count only one party

        // for(int j=0; j<m; ++j) {
        // 	delete[] x[j];
        // 	delete[] y[j];
        // }
        delete[] x;
        delete[] y;
    }


}

void phase3(int num_trial) {
    // secure class generation via path cost
//    TODO
    int n, d;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        d = param_nd[i][1];
        int num_edge = pow(2, d + 1) - 1, num_leaf = pow(2, d);

        std::vector<mpz_class> edges(num_edge, 0), leaf_value(num_leaf, 0), interm_rlt(num_leaf - 1, 0), path_cost(
                num_leaf, 0);
        for (int j = 0; j < num_edge; ++j)
            edges[j] = gmp_prn.get_z_range(CONFIG_P);
        for (int j = 0; j < num_leaf; ++j)
            leaf_value[j] = gmp_prn.get_z_range(CONFIG_P);

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            secure_class_generation_path_cost(edges, leaf_value, interm_rlt, path_cost, d);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }

        printf("secure class generation via path cost (n=%d, d=%d): %f ns\n", n, d,
               time_total / num_trial); // one party
    }
}

void phase4(int num_trial) {
    // secure class generation via polynomial
    int n, d;
//    TODO
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        d = param_nd[i][1];
        int num_edge = pow(2, d + 1) - 1, num_leaf = pow(2, d);

        mpz_class (*edges)[2] = new mpz_class[num_edge][2], (*leaf_value)[2] = new mpz_class[num_leaf][2], (*interm_rlt)[2] = new mpz_class[
        num_leaf - 1][2], (*path_mul)[2] = new mpz_class[num_leaf][2];
        for (int j = 0; j < num_edge; ++j) {
            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), edges[j]);
        }

        for (int j = 0; j < num_leaf; ++j) {
            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), leaf_value[j]);
            if (j != (num_leaf - 1)) {
                interm_rlt[j][0] = interm_rlt[j][1] = 0;
            }
            path_mul[j][0] = path_mul[j][1] = 0;
        }

        triplet_z tri;

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            secure_class_generation_polynomial(edges, leaf_value, interm_rlt, path_mul, tri, d);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }

        delete[] edges;
        delete[] leaf_value;
        delete[] interm_rlt;
        delete[] path_mul;

        printf("secure class generation via polynomial (n=%d, d=%d): %f ns\n", n, d, time_total / num_trial / 2);
    }
}