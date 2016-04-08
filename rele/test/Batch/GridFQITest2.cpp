/*
 * rele,
 *
 *
 * Copyright (C) 2015  Davide Tateo & Matteo Pirotta
 * Versione 1.0
 *
 * This file is part of rele.
 *
 * rele is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rele is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rele.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rele/core/Core.h"
#include "rele/core/BatchCore.h"
#include "rele/core/PolicyEvalAgent.h"
#include "rele/policy/q_policy/e_Greedy.h"
#include "rele/algorithms/batch/td/W-FQI.h"
#include "rele/algorithms/td/Q-Learning.h"
#include "rele/approximators/features/DenseFeatures.h"
#include "rele/core/FiniteMDP.h"
#include "rele/generators/GridWorldGenerator.h"
#include "rele/approximators/regressors/nn/FFNeuralNetwork.h"
#include "rele/approximators/regressors/nn/FFNeuralNetworkEnsemble.h"
#include "rele/approximators/basis/IdentityBasis.h"
#include "rele/approximators/regressors/trees/KDTree.h"
#include "rele/approximators/regressors/trees/ExtraTree.h"
#include "rele/approximators/regressors/trees/ExtraTreeEnsemble.h"

#include <iostream>

using namespace std;
using namespace ReLe;
using namespace arma;


// This simple test is used to verify the correctness of the FQI implementation
int main(int argc, char *argv[])
{
    GridWorldGenerator generator;
    generator.load(argv[1]);

    double gamma = 0.9;
    FiniteMDP&& mdp = generator.getMDP(gamma);

    unsigned int nActions = mdp.getSettings().actionsNumber;
    unsigned int nStates = mdp.getSettings().statesNumber;

    BasisFunctions bfs;
    bfs = IdentityBasis::generate(2);

    DenseFeatures phi(bfs);

    //FFNeuralNetwork QRegressor(phi, 50, 1);
    //QRegressor.getHyperParameters().lambda = 0.0005;
    //QRegressor.getHyperParameters().maxIterations = 10;

    arma::vec defaultValue = {0};
    EmptyTreeNode<arma::vec> defaultNode(defaultValue);
    KDTree<arma::vec, arma::vec> QRegressor(phi, defaultNode, 1, 1);

    FiniteW_FQI batchAgent(QRegressor, nStates, nActions, 1e-8);

    auto&& core = buildBatchCore(mdp, batchAgent);

    core.getSettings().nEpisodes = 100;
    core.getSettings().agentLogger = new BatchAgentPrintLogger<FiniteAction, FiniteState>();
    core.getSettings().episodeLength = 100;
    core.getSettings().maxBatchIterations = 100;

    arma::mat Q(nStates, nActions, arma::fill::zeros);
    e_Greedy policy;
    policy.setQ(&Q);

    core.run(policy);

    return 0;
}
