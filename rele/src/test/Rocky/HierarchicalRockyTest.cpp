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

#include "Rocky.h"

#include "Core.h"

#include "policy_search/gradient/hierarchical/HierarchicalGPOMDP.h"

#include "parametric/differentiable/GibbsPolicy.h"
#include "features/DenseFeatures.h"
#include "basis/IdentityBasis.h"

#include "FileManager.h"
#include "ConsoleManager.h"

#include "RockyOptions.h"

#include <iostream>

using namespace std;
using namespace ReLe;
using namespace arma;


int main(int argc, char *argv[])
{
    FileManager fm("Rocky", "HPG");
    fm.createDir();
    fm.cleanDir();

    Rocky rocky;


    //-- Features
    BasisFunctions basis = IdentityBasis::generate(rocky.getSettings().continuosStateDim);
    DenseFeatures phi(basis);


    //-- options
    vector<FiniteAction> actions;
    for(int i = 0; i < 4; ++i)
        actions.push_back(FiniteAction(i));


    Eat eat;
    Home home;
    Feed feed;
    Escape escape;

    vector<Option<DenseAction, DenseState>*> options;
    options.push_back(&eat);
    options.push_back(&home);
    options.push_back(&feed);
    options.push_back(&escape);

    ParametricGibbsPolicy<DenseState> rootPolicyOption(actions, phi, 5);
    DifferentiableOption<DenseAction, DenseState> rootOption(rootPolicyOption, options);

    //-- agent
    int nbepperpol = 25, nbstep = 20000;
    AdaptiveStep stepRule(0.01);
    HierarchicalGPOMDPAlgorithm<DenseAction, DenseState> agent(rootOption, nbepperpol, nbstep, stepRule);

    Core<DenseAction, DenseState> core(rocky, agent);
    //--


    int episodes = 3000;
    core.getSettings().episodeLenght = 10000;
    core.getSettings().loggerStrategy = new WriteStrategy<DenseAction, DenseState>(fm.addPath("Rocky.log"),
            WriteStrategy<DenseAction, DenseState>::outType::AGENT);


    ConsoleManager console(episodes, 1);
    console.printInfo("starting learning");
    for (int i = 0; i < episodes; i++)
    {
        console.printProgress(i);
        core.runEpisode();
    }

    delete core.getSettings().loggerStrategy;

    console.printInfo("Starting evaluation episode");
    core.getSettings().loggerStrategy = new WriteStrategy<DenseAction, DenseState>(fm.addPath("Rocky.log"),
            WriteStrategy<DenseAction, DenseState>::outType::TRANS);
    core.runTestEpisode();

    delete core.getSettings().loggerStrategy;

    return 0;

}