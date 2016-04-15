/*
 * BicycleTest.cpp
 *
 *  Created on: 15 apr 2016
 *      Author: samuele
 */

#include "rele/environments/Bicycle.h"
#include "rele/core/Core.h"
#include "rele/algorithms/td/LinearSARSA.h"
#include "rele/algorithms/td/DenseSARSA.h"
#include "rele/approximators/basis/PolynomialFunction.h"
#include "rele/approximators/features/DenseFeatures.h"
#include "rele/approximators/basis/GaussianRbf.h"
#include "rele/approximators/basis/ConditionBasedFunction.h"
#include "rele/policy/q_policy/e_Greedy.h"
#include "rele/utils/FileManager.h"

using namespace std;
using namespace ReLe;

int main(int argc, char *argv[])
{
    unsigned int episodes = 10000;
    Bicycle mdp;

    BasisFunctions bVector = PolynomialFunction::generate(1, mdp.getSettings().statesNumber + 1);
    BasisFunctions basis = AndConditionBasisFunction::generate(bVector, 100, mdp.getSettings().actionsNumber);

    DenseFeatures phi(basis);

    e_GreedyApproximate policy;
    policy.setEpsilon(0.05);
    ConstantLearningRateDense alpha(0.1);
    LinearGradientSARSA agent(phi, policy, alpha);
    agent.setLambda(0.8);

    FileManager fm("bc", "linearSarsa");
    fm.createDir();
    fm.cleanDir();
    auto&& core = buildCore(mdp, agent);
    core.getSettings().loggerStrategy = new WriteStrategy<FiniteAction, DenseState>(fm.addPath("bc.log"));

    for (int i = 0; i < episodes; i++)
    {
        core.getSettings().episodeLength = 10000;
        cout << "Starting episode: " << i << endl;
        core.runEpisode();
    }

}




