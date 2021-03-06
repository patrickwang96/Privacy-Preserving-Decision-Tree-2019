//
// Created by Ruochen WANG on 16/4/2020.
//

#ifndef PRIVACY_PRESERVING_DECISION_TREE_2019_TEST_CLOUD_SERVER_H
#define PRIVACY_PRESERVING_DECISION_TREE_2019_TEST_CLOUD_SERVER_H

#include <vector>
#include "network.h"

void test_by_phases(std::vector<int> phases, int num_trial);

void phase1(int num_trial, NetAdapter *net);

void phase2(int num_trial, NetAdapter *net);

void phase25(int num_trial, NetAdapter *net);

void phase3(int num_trial, NetAdapter *net);

void phase4(int num_trial, NetAdapter *net);

#endif //PRIVACY_PRESERVING_DECISION_TREE_2019_TEST_CLOUD_SERVER_H
