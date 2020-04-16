//
// Created by Ruochen WANG on 16/4/2020.
//

#include "server_decision_tree.h"
#include "secret_sharing_efficient_tools.h"

extern gmp_randclass gmp_prn;

void server_secure_multiplication(const matrix_z shareA, const matrix_z shareB,
                                  const triplet_mz shareTri,
        // intermediate buffers
                                  ss_tuple_mz &U, ss_tuple_mz &V,
        // output
                                  matrix_z shareAB,
        // whether to conduct piecewise multiplication
                                  int piecewise, NetAdapter *net) {
    // 1)
    U.share[0] = shareA - shareTri.X;
    V.share[0] = shareB - shareTri.Y;

    // 2)
    ss_decrypt_server(U.plain, U.share[0], net);
    ss_decrypt_server(V.plain, V.share[0], net);

    // 3)

    // for server (0)
    shareAB = U.plain * shareTri.Y;
    shareAB += shareTri.X * V.plain;
    shareAB += shareTri.Z;


    mod_2exp(shareAB, CONFIG_L);

}

void server_secure_node_evaluation(std::vector<mpz_class> &x, std::vector<mpz_class> &y, const triplet_z &tri_z,
                                   const triplet_b &tri_b,
                                   NetAdapter *net) {
    int m = x.size();

    // 1)
    auto a = new mpz_class[m];

    for (int i = 0; i < m; i++) {
        a[i] = y[i] - x[i];
        mod_2exp(a[i], CONFIG_L);
    }

    // 2)
    auto p = new int[m * CONFIG_L];
    memset(p, 0, m * CONFIG_L * sizeof(int));
    auto q = new int[m * CONFIG_L];

    for (int i = 0; i < CONFIG_L; i++) {
        for (int j = 0; j < m; j++) {
            q[i * m + j] = extract_bit(a[j], i);
        }
    }

    // 3)
    auto c = new int[m * CONFIG_L];
    // compute only the 1st bit
    secure_mul_server_batch_compressed(p, q, c, m, tri_b, net);

    // 4)
    auto d = new int[m * CONFIG_L];
    secure_mul_server_batch_compressed(p, q, d, m, tri_b, net);
    for (int i = 0; i < m * CONFIG_L; i++) d[i] ^= 1;

    auto e = new int[m * CONFIG_L];


    for (int i = 1; i < CONFIG_L - 1; i++) {
        secure_mul_server_batch_compressed(q + i * m, c + (i - 1) * m, e + i * m, m, tri_b, net);
        for (int j = 0; j < m; j++) e[i * m + j] ^= 1; // plus one
        secure_mul_server_batch_compressed(e + i * m, d + i * m, c + i * m, m, tri_b, net);
        for (int j = 0; j < m; j++) c[i * m + j] ^= 1; // plus one
    }

    // 5)
    auto a_2 = new int[m];
    for (int i = 0; i < m; i++) a_2[i] = q[m * (CONFIG_L - 1) + i] ^ c[m * (CONFIG_L - 2) + i];

    delete[] a;
    delete[] p;
    delete[] q;
    delete[] c;
    delete[] d;
    delete[] e;

    // 6)
//    std::vector<mpz_class> t1(m, 0);
//    std::vector<mpz_class> t2(m);
    auto t1 = new uint64_t[m];
    auto t2 = new uint64_t[m];

    for (int i = 0; i < m; i++) {
        t1[i] = 0;
        t2[i] = a_2[i];
    }

    // 7)

    auto result_a = new mpz_class[m];
    auto product = new uint64_t[m];

    secure_mul_server_batch(t1, t2, product, m, tri_b, net);
    for (int i = 0; i < m; i++) {
        result_a[i] = (long) (t1[i] + t2[i] - 2 * product[i]);
    }

    delete[] t1;
    delete[] t2;
    delete[] product;
    delete[] result_a; // this is the result;
}

void server_secure_class_generation_path_cost(const std::vector<mpz_class> &edges, std::vector<mpz_class> &leaf_value,
                                              std::vector<mpz_class> &interm_rlt, std::vector<mpz_class> &path_cost,
                                              int depth, NetAdapter *net) {
    int first_leaf = pow(2, depth) - 1;
    int num_leaf = pow(2, depth);

    // non-leaf nodes
    int id = 0;
    for (id = 1; id < first_leaf; ++id)
        interm_rlt[id] += interm_rlt[(id - 1) / 2] + edges[id];

    // leaf nodes
    for (int i = 0; i < num_leaf; ++i, ++id)
        path_cost[i] = interm_rlt[(id - 1) / 2] + edges[id];

    /* masking */
    for (int i = 0; i < num_leaf; ++i) {
        leaf_value[i] += gmp_prn.get_z_range(CONFIG_P) * path_cost[i];
        mod_prime(leaf_value[i], CONFIG_P);

        path_cost[i] *= gmp_prn.get_z_range(CONFIG_P);
        mod_prime(path_cost[i], CONFIG_P);
    }

    // permutation
    std::random_shuffle(path_cost.begin(), path_cost.end());
    std::random_shuffle(leaf_value.begin(), leaf_value.end());


}

void
server_secure_class_generation_polynomial(std::vector<mpz_class> &edges, std::vector<mpz_class> &leaf_value,
                                          std::vector<mpz_class> &interm_rlt,
                                          std::vector<mpz_class> &path_mul, const triplet_b &tri, int depth,
                                          NetAdapter *net) {

    // complete traversal of the binary tree
    //             0
    //       1            2          depth=1
    //     3   4       5     6       depth=2
    //    7 8 9 10   11 12 13 14     depth=3
    //   ...      ...         ...    depth=4 (leaves)

    //   egde_id goes with node_id

    int first_leaf = pow(2, depth) - 1;
    int num_leaf = pow(2, depth);

    // non-leaf nodes
    int id = 0;
    interm_rlt[0] = edges[0];
    auto tmp1 = new uint64_t[first_leaf];
    auto tmp2 = new uint64_t[first_leaf];
    auto buffer = new uint64_t[first_leaf];

    id = 1;
    for (int d = 1; d < depth; d++) {
        int layer_count = (1 << d);
        for (int l = 0; l < layer_count; l++) {
            tmp1[l] = mpz_to_u64(interm_rlt[(id - 1) / 2]);
            tmp2[l] = mpz_to_u64(edges[id]);
            id++;
        }
        secure_mul_server_batch(tmp1, tmp2, buffer, layer_count, tri, net);
        for (int l = id - layer_count; l < id; l++)
            interm_rlt[l] = (long) buffer[l - (id - layer_count)];

    }
    delete[] tmp1;
    delete[] tmp2;
    delete[] buffer;


    // leaf nodes
    auto iterm_rlt_uint = new uint64_t[num_leaf];
    auto edges_uint = new uint64_t[num_leaf];
    auto path_mul_uint = new uint64_t[num_leaf];
    id = first_leaf;
    for (int i = 0; i < num_leaf; i++, id++) {
        iterm_rlt_uint[i] = mpz_to_u64(interm_rlt[(id - 1) / 2]);
        edges_uint[i] = mpz_to_u64(edges[id]);
    }
    secure_mul_server_batch(iterm_rlt_uint, edges_uint, path_mul_uint, num_leaf, tri, net);
    for (int i = 0; i < num_leaf; i++) path_mul[i] = (long) path_mul_uint[i];

    uint64_t final = 0; // final value;
    auto leaf_value_uint = new uint64_t[num_leaf];
    for (int i = 0; i < num_leaf; i++) leaf_value_uint[i] = mpz_to_u64(leaf_value[i]);
    secure_mul_server_batch(path_mul_uint, leaf_value_uint, leaf_value_uint, num_leaf, tri, net);
    for (int i = 0; i < num_leaf; i++) final += leaf_value_uint[i];

    delete[] iterm_rlt_uint;
    delete[] edges_uint;
    delete[] path_mul_uint;
}
