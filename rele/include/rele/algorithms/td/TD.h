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

#ifndef INCLUDE_ALGORITHMS_TD_TD_H_
#define INCLUDE_ALGORITHMS_TD_TD_H_

#include "Agent.h"
#include <armadillo>
#include "LinearApproximator.h"

namespace ReLe
{

class FiniteTD: public Agent<FiniteAction, FiniteState>
{
public:
    FiniteTD();

    void setAlpha(double alpha)
    {
        this->alpha = alpha;
    }

    void setEpsilon(double eps)
    {
        this->eps = eps;
    }

protected:
    virtual void init();
    unsigned int policy(std::size_t x);
    void printStatistics();

protected:
    //Action-value function
    arma::mat Q;

    //current an previous actions and states
    size_t x;
    unsigned int u;

    //algorithm parameters
    double alpha;
    double eps;

};

class LinearTD : public Agent<FiniteAction, DenseState>
{
public:
    LinearTD(LinearApproximator &la);

    void setAlpha(double alpha)
    {
        this->alpha = alpha;
    }

    void setEpsilon(double eps)
    {
        this->eps = eps;
    }

    void setLinearApproximator(LinearApproximator &la)
    {
        this->Q = la;
    }

protected:
    unsigned int policy(DenseState state);
    void printStatistics();

protected:
    //Linear action-value function
    LinearApproximator Q;

    //current an previous actions and states
    DenseState x;
    unsigned int u;

    //algorithm parameters
    double alpha;
    double eps;
};

}//end namespace

#endif /* INCLUDE_ALGORITHMS_TD_TD_H_ */
