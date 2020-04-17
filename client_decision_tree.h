//
// Created by Ruochen WANG on 16/4/2020.
//

#ifndef PRIVACY_PRESERVING_DECISION_TREE_2019_CLIENT_DECISION_TREE_H
#define PRIVACY_PRESERVING_DECISION_TREE_2019_CLIENT_DECISION_TREE_H

#include "types.h"
#include "network.h"
#include "utils.h"
#include <vector>

void client_secure_multiplication(const matrix_z shareA, const matrix_z shareB,
                                  const triplet_mz shareTri,
        // intermediate buffers
                                  ss_tuple_mz &U, ss_tuple_mz &V,
        // output
                                  matrix_z shareAB,
        // whether to conduct piecewise multiplication
                                  int piecewise, NetAdapter *net);

void client_secure_node_evaluation(std::vector<mpz_class> &x, std::vector<mpz_class> &y, const triplet_z &tri_z,
                                   const triplet_b &tri_b,
                                   NetAdapter *net);

void client_secure_class_generation_path_cost(const std::vector<mpz_class> &edges, std::vector<mpz_class> &leaf_value,
                                              std::vector<mpz_class> &interm_rlt, std::vector<mpz_class> &path_cost,
                                              int depth, NetAdapter *net);

void fig6(std::vector<mpz_class> &leaf_value, std::vector<mpz_class> &path_cost, int depth, NetAdapter *net);

void
client_secure_class_generation_polynomial(std::vector<mpz_class> &edges, std::vector<mpz_class> &leaf_value,
                                          std::vector<mpz_class> &interm_rlt,
                                          std::vector<mpz_class> &path_mul, const triplet_b &tri, int depth,
                                          NetAdapter *net);

#endif //PRIVACY_PRESERVING_DECISION_TREE_2019_CLIENT_DECISION_TREE_H
