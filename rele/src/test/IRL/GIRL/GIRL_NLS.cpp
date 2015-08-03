/*
 * rele,
 *
 *
 * Copyright (C) 2015 Davide Tateo & Matteo Pirotta
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

#include "features/SparseFeatures.h"
#include "features/DenseFeatures.h"
#include "regressors/GaussianMixtureModels.h"
#include "basis/IdentityBasis.h"
#include "basis/GaussianRbf.h"

#include "parametric/differentiable/NormalPolicy.h"

#include "NLS.h"

#include "PolicyEvalAgent.h"
#include "Core.h"
#include "ParametricRewardMDP.h"
#include "algorithms/GIRL.h"
#include "algorithms/PGIRL.h"
#include "policy_search/gradient/onpolicy/GPOMDPAlgorithm.h"

#include "FileManager.h"

using namespace std;
using namespace arma;
using namespace ReLe;

int main(int argc, char *argv[])
{
//  RandomGenerator::seed(45423424);
//  RandomGenerator::seed(8763575);

    IRLGradType atype = IRLGradType::R;
    int nbEpisodes = 10000;

    FileManager fm("nls", "GIRL");
    fm.createDir();
    fm.cleanDir();
    std::cout << std::setprecision(OS_PRECISION);


    NLS mdp;

    //Setup expert policy
    int dim = mdp.getSettings().continuosStateDim;

    BasisFunctions basis = IdentityBasis::generate(dim);
    DenseFeatures phi(basis);

    BasisFunctions stdBasis = IdentityBasis::generate(dim);
    DenseFeatures stdPhi(stdBasis);
    arma::vec stdWeights(stdPhi.rows());
    stdWeights.fill(0.5);

    NormalStateDependantStddevPolicy expertPolicy(phi, stdPhi, stdWeights);

    arma::vec p(2);
    p(0) = 6.5178;
    p(1) = -2.5994;

    expertPolicy.setParameters(p);

    PolicyEvalAgent<DenseAction, DenseState> expert(expertPolicy);

    // Generate expert dataset
    Core<DenseAction, DenseState> expertCore(mdp, expert);
    CollectorStrategy<DenseAction, DenseState> collection;
    expertCore.getSettings().loggerStrategy = &collection;
    expertCore.getSettings().episodeLenght = mdp.getSettings().horizon;
    expertCore.getSettings().testEpisodeN = nbEpisodes;
    expertCore.runTestEpisodes();
    Dataset<DenseAction,DenseState>& data = collection.data;


    // Create parametric reward
    //BasisFunctions basisReward = IdentityBasis::generate(2);
    arma::vec cT = {0, 0};
    arma::vec cF = {10, 10};
    BasisFunction* bfT = new GaussianRbf(cT, 1);
    BasisFunction* bfF = new GaussianRbf(cF, 1);
    BasisFunctions basisReward;
    basisReward.push_back(bfT);
    basisReward.push_back(bfF);

    /*unsigned int nbasis = 3;
    for(unsigned int i = 0; i < nbasis; i++)
    {
        double angle = i*2.0/3.0*M_PI;
        basisReward.push_back(new GaussianRbf({cos(angle), sin(angle)}, 1));
    }*/

    DenseFeatures phiReward(basisReward);

    //GaussianRegressor rewardRegressor(phiReward);
    LinearApproximator rewardRegressor(phiReward);
    /*PlaneGIRL<DenseAction,DenseState> irlAlg(data, expertPolicy, basisReward,
            mdp.getSettings().gamma, atype);*/
    GIRL<DenseAction,DenseState> irlAlg(data, expertPolicy, rewardRegressor,
                                        mdp.getSettings().gamma, atype, true);

    //Info print
    std::cout << "Basis size: " << phiReward.rows();
    std::cout << " | Params: " << expertPolicy.getParameters().t() << std::endl;
    std::cout << "Features Expectation " << data.computefeatureExpectation(phiReward, mdp.getSettings().gamma).t();

    //Run GIRL
    irlAlg.run();
    arma::vec weights = irlAlg.getWeights();


    rewardRegressor.setParameters(weights);
    cout << "weights: " << weights.t();


    //Try to recover the initial policy
    int episodesPerPolicy = 1;
    int policyPerUpdate = 100;
    int updates = 400;
    int episodes = episodesPerPolicy*policyPerUpdate*updates;

    NormalStateDependantStddevPolicy imitatorPolicy(phi, stdPhi, stdWeights);
    AdaptiveStep stepRule(0.01);
    int nparams = phi.rows();
    arma::vec mean(nparams, fill::zeros);
    mean[0] = -0.42;
    mean[1] =  0.42;

    imitatorPolicy.setParameters(mean);
    GPOMDPAlgorithm<DenseAction, DenseState> imitator(imitatorPolicy, policyPerUpdate,
            mdp.getSettings().horizon, stepRule, GPOMDPAlgorithm<DenseAction, DenseState>::BaseLineType::MULTI);

    ParametricRewardMDP<DenseAction, DenseState> prMdp(mdp, rewardRegressor);
    Core<DenseAction, DenseState> imitatorCore(prMdp, imitator);
    EmptyStrategy<DenseAction, DenseState> emptyStrategy;
    imitatorCore.getSettings().loggerStrategy = &emptyStrategy;
    imitatorCore.getSettings().episodeLenght = mdp.getSettings().horizon;
    imitatorCore.getSettings().episodeN = episodes;
    imitatorCore.getSettings().testEpisodeN = nbEpisodes;
    imitatorCore.runEpisodes();

    cout << "Learned Parameters: " << imitatorPolicy.getParameters().t() << endl;
    cout << arma::as_scalar(imitatorCore.runBatchTest()) << endl;

    //Evaluate policy against the real mdp
    Core<DenseAction, DenseState> evaluationCore(mdp, imitator);
    CollectorStrategy<DenseAction, DenseState> collector2;
    evaluationCore.getSettings().loggerStrategy = &collector2;
    evaluationCore.getSettings().episodeLenght = mdp.getSettings().horizon;
    evaluationCore.getSettings().episodeN = episodes;
    evaluationCore.getSettings().testEpisodeN = nbEpisodes;

    evaluationCore.runTestEpisodes();

    Dataset<DenseAction,DenseState>& data2 = collector2.data;

    double gamma = mdp.getSettings().gamma;
    cout << "Features Expectation ratio: " << (data2.computefeatureExpectation(phiReward, gamma)/data.computefeatureExpectation(phiReward, gamma)).t();
    cout << "reward: " << arma::as_scalar(evaluationCore.runBatchTest());

    //calculate full grid function
    int samplesPositive = 100;
    arma::vec valuesG(samplesPositive);
    arma::vec valuesJ(samplesPositive);

#ifdef PRINT_GRADIENT
    for(int i = 0; i < samplesPositive; i++)
    {
        double step = 0.01;
        arma::vec wm(2);
        wm(0) = i*step;
        wm(1) = 1.0 - wm(0);
        rewardRegressor.setParameters(wm);
        arma::mat gGrad(2, 2);
        arma::vec grad = irlAlg.ReinforceBaseGradient(gGrad);
        //gradients.row(i+samplesPositive) = gGrad.t()*grad;
        double Je = irlAlg.computeJ();
        valuesG(i) = as_scalar(0.5*grad.t()*grad);
        valuesJ(i) = 0.5*std::pow(Je, 2);
        cout << i << endl;
    }

    valuesG.save("/tmp/ReLe/normG.txt", arma::raw_ascii);
    valuesJ.save("/tmp/ReLe/normJ.txt", arma::raw_ascii);
    //gradients.save("/tmp/ReLe/gradient.txt", arma::raw_ascii);
#endif

    return 0;
}
