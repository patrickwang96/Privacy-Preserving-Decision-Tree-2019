//
// Created by Ruochen WANG on 16/4/2020.
//

#include "test_cloud_client.h"
#include "secret_sharing.h"
#include "utils.h"
#include "secret_sharing_efficient_tools.h"
#include "client_decision_tree.h"
#include "network.h"
#include <chrono>
#include <cmath>
#include <iostream>

extern gmp_randclass gmp_prn;

auto start = std::chrono::steady_clock::now(), end = std::chrono::steady_clock::now();
#define CLOCK_START {start = std::chrono::steady_clock::now();}
#define CLOCK_END {end = std::chrono::steady_clock::now();}
#define ELAPSED std::chrono::duration<double, std::nano>(end - start).count()

void test_by_phases(std::vector<int> phases, int num_trial) {
    auto *net = new NetAdapter(1);
    for (auto phase : phases) {
        if (phase == 1) {
            phase1(num_trial, net);
        } else if (phase == 2) {
            phase2(num_trial, net);
        } else if (phase == 3) {
            phase25(num_trial, net);
        } else if (phase == 4) {
            phase3(num_trial, net);
        } else if (phase == 5) {
            phase4(num_trial, net);
        } else {
            std::cout << "Phase No. out of range." << std::endl;
        }
    }
    net->close();
}

void phase1(int num_trial, NetAdapter *net) {
    // node selection
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
        uint64_t s_bytes = net->get_send_bytes();
        uint64_t r_bytes = net->get_rev_bytes();

        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            client_secure_multiplication(sel_mat.share[1], feature.share[1],
                                         tri.share[1],
                                         buf_E, buf_f,
                                         xsigma.share[1],
                                         0, net);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }
        s_bytes = net->get_send_bytes() - s_bytes;
        r_bytes = net->get_rev_bytes() - r_bytes;

        double send_mb = 1.0 * s_bytes / num_trial / 1024 / 1024, recv_mb = 1.0 * r_bytes / num_trial / 1024 / 1024;


        printf("secure node selection (n=%d, d=%d, m=%d): %f ns, s bytes: %f, r bytes: %f\n", n, param_nd[i][1], m,
               time_total / num_trial / 2, send_mb, recv_mb); // count only one party
    }
}

void phase2(int num_trial, NetAdapter *net) {
    // node evaluation
    int n, m;
// secure node evaluation
    triplet_b tri_b;
    triplet_z tri_z;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        m = pow(2, param_nd[i][1]) - 1;

        // ugly staff for preparation
//        mpz_class (*x)[2] = new mpz_class[m][2], (*y)[2] = new mpz_class[m][2];
//        for (int j = 0; j < m; ++j) {
//            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), x[j]);
//            ss_encrypt(gmp_prn.get_z_bits(CONFIG_L), y[j]);
//        }
        std::vector<mpz_class> x(m);
        std::vector<mpz_class> y(m);
        for (int j = 0; j < m; j++) {
            x[i] = gmp_prn.get_z_bits(CONFIG_L);
            y[i] = gmp_prn.get_z_bits(CONFIG_L);
        }

        uint64_t s_bytes = net->get_send_bytes();
        uint64_t r_bytes = net->get_rev_bytes();

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
//            for (int k = 0; k < m; ++k)
            client_secure_node_evaluation(x, y, tri_z, tri_b, net);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }
        s_bytes = net->get_send_bytes() - s_bytes;
        r_bytes = net->get_rev_bytes() - r_bytes;

        double send_mb = 1.0 * s_bytes / num_trial / 1024 / 1024, recv_mb = 1.0 * r_bytes / num_trial / 1024 / 1024;

        printf("secure node evaluation (n=%d, d=%d, m=%d): %f ns, s bytes: %f, r bytes: %f\n", n, param_nd[i][1], m,
               time_total / num_trial / 2, send_mb, recv_mb); // count only one party

    }


}


void phase25(int num_trial, NetAdapter *net) {
    // node evaluation
    int n, m;
// secure node evaluation
    triplet_b tri_b;
    triplet_z tri_z;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        m = pow(2, param_nd[i][1]) - 1;

        // ugly staff for preparation
        std::vector<int> a(m);
        std::vector<mpz_class> result(m);
        for (int j = 0; j < m; j++) {
            result[i] = gmp_prn.get_z_bits(1);
            a[i] = result[i].get_ui();
            result[i] = gmp_prn.get_z_bits(CONFIG_L);
        }

        uint64_t s_bytes = net->get_send_bytes();
        uint64_t r_bytes = net->get_rev_bytes();

        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
//            for (int k = 0; k < m; ++k)
            client_node_eval_phase2(a, result, tri_b, net);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }
        s_bytes = net->get_send_bytes() - s_bytes;
        r_bytes = net->get_rev_bytes() - r_bytes;

        double send_mb = 1.0 * s_bytes / num_trial / 1024 / 1024, recv_mb = 1.0 * r_bytes / num_trial / 1024 / 1024;

        printf("secure node evaluation phase 2 (n=%d, d=%d, m=%d): %f ns, bytes: %f\n", n, param_nd[i][1], m,
               time_total / num_trial / 2, send_mb + recv_mb); // count only one party

    }

}


void phase3(int num_trial, NetAdapter *net) {
    // secure class generation via path cost
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

        uint64_t s_bytes = net->get_send_bytes();
        uint64_t r_bytes = net->get_rev_bytes();
        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            client_secure_class_generation_path_cost(edges, leaf_value, interm_rlt, path_cost, d, net);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }
        s_bytes = net->get_send_bytes() - s_bytes;
        r_bytes = net->get_rev_bytes() - r_bytes;

        double send_mb = 1.0 * s_bytes / num_trial / 1024 / 1024, recv_mb = 1.0 * r_bytes / num_trial / 1024 / 1024;

        printf("secure class generation via path cost (n=%d, d=%d): %f ns, s bytes: %f, r bytes: %f\n", n, d,
               time_total / num_trial, send_mb, recv_mb); // one party
    }
}

void phase4(int num_trial, NetAdapter *net) {
    // secure class generation via polynomial
    int n, d;
    for (int i = 0; i < 5; ++i) {
        n = param_nd[i][0];
        d = param_nd[i][1];
        int num_edge = pow(2, d + 1) - 1, num_leaf = pow(2, d);

        std::vector<mpz_class> edges(num_edge);
        std::vector<mpz_class> leaf_value(num_leaf);
        std::vector<mpz_class> interm_rlt(num_leaf - 1, 0);
        std::vector<mpz_class> path_cost(num_leaf, 0);

        for (int j = 0; j < num_edge; ++j) {
            edges[i] = gmp_prn.get_z_bits(CONFIG_L);
        }

        for (int j = 0; j < num_leaf; ++j) {
            leaf_value[i] = gmp_prn.get_z_bits(CONFIG_L);
        }

        triplet_b tri;

        uint64_t s_bytes = net->get_send_bytes();
        uint64_t r_bytes = net->get_rev_bytes();
        double time_total = 0;
        for (int j = 0; j < num_trial; ++j) {
            CLOCK_START
            fig6(leaf_value, path_cost, d, net);
            CLOCK_END
            time_total += ELAPSED;

            cache_flusher();
        }
        s_bytes = net->get_send_bytes() - s_bytes;
        r_bytes = net->get_rev_bytes() - r_bytes;

        double send_mb = 1.0 * s_bytes / num_trial / 1024 / 1024, recv_mb = 1.0 * r_bytes / num_trial / 1024 / 1024;

        printf("secure class generation via polynomial (n=%d, d=%d): %f ns, s bytes: %f, r bytes: %f\n", n, d,
               time_total / num_trial / 2, send_mb, recv_mb);
    }
}