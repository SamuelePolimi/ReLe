/*
 * rele,
 *
 *
 * Copyright (C) 2015 Matteo Pirotta & Davide Tateo
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

#ifndef INCLUDE_ALGORITHMS_TD_LINEARSARSA_H_
#define INCLUDE_ALGORITHMS_TD_LINEARSARSA_H_

#include "td/TD.h"
#include <armadillo>

namespace ReLe
{

class LinearGradientSARSA: public LinearTD
{
public:
    LinearGradientSARSA(LinearApproximator &la);
    virtual void initEpisode(const DenseState& state, FiniteAction& action);
    virtual void sampleAction(const DenseState& state, FiniteAction& action);
    virtual void step(const Reward& reward, const DenseState& nextState,
                      FiniteAction& action);
    virtual void endEpisode(const Reward& reward);
    virtual void endEpisode();

    virtual ~LinearGradientSARSA();

protected:
    virtual void printStatistics();

private:
    arma::vec eligibility;
    double lambda;

};

}

#endif /* INCLUDE_ALGORITHMS_TD_LINEARSARSA_H_ */