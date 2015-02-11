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

#include "policy_search/REPS.h"

using namespace arma;

namespace ReLe
{

TabularREPS::TabularREPS()
{
	x = 0;
	u = 0;
	eta = 1;

	//default parameters
	N = 1;
	eps = 0.5;

	//sample iteration counter
	currentIteration = 0;
}

void TabularREPS::initEpisode(const FiniteState& state, FiniteAction& action)
{
	x = state.getStateN();
	u = policy(x);

	action.setActionN(u);

	resetSamples();
}

void TabularREPS::sampleAction(const FiniteState& state, FiniteAction& action)
{
	x = state.getStateN();
	u = policy(x);

	action.setActionN(u);
}

void TabularREPS::step(const Reward& reward, const FiniteState& nextState,
			FiniteAction& action)
{

	size_t xn = nextState.getStateN();
	unsigned int un = policy(xn);
	double r = reward[0];

	updateSamples(xn, r);

	if (currentIteration >= N)
	{
		currentIteration = 0;
		updatePolicy();
		resetSamples();
	}

	//update action and state
	x = xn;
	u = un;

	//set next action
	action.setActionN(u);
}

void TabularREPS::endEpisode(const Reward& reward)
{
	updatePolicy();
	printStatistics();
}

void TabularREPS::endEpisode()
{
	printStatistics();
}

TabularREPS::~TabularREPS()
{

}

void TabularREPS::updatePolicy()
{

}

void TabularREPS::updateSamples(size_t xn, double r)
{
	ndelta(x, u) += r + theta[xn] - theta[x];
	nlambda(x, u) += 0; //TODO FIXME
	d(x, u) += 1;

	currentIteration++;
}

void TabularREPS::resetSamples()
{
	ndelta.zeros(task.finiteStateDim, task.finiteActionDim);
	nlambda.zeros(task.finiteStateDim, task.finiteActionDim);
	d.zeros(task.finiteStateDim, task.finiteActionDim);

	currentIteration = 0;
}

void TabularREPS::init()
{
	//Init policy and parameters
	policy.init(task.finiteStateDim, task.finiteActionDim);
	theta = vec(task.finiteStateDim, fill::ones);
	eta = 1;

	//setup optimization algorithm
	optimizator = nlopt::opt(nlopt::algorithm::LD_LBFGS, theta.n_elem + 1);
}

void TabularREPS::printStatistics()
{
	cout << endl << endl << "### Tabular REPS ###";
	cout << endl << endl << "Using " << policy.getPolicyName() << " policy"
				<< endl << endl;

	cout << "--- Parameters ---" << endl << endl;
	cout << "N: " << N << endl;
	cout << "eps: " << eps << endl;

	cout << endl << endl << "--- Learning results ---" << endl << endl;
	cout << "- Policy" << endl;
	cout << policy.printPolicy();
}

}
