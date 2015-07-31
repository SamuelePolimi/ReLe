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

#ifndef INCLUDE_RELE_ENVIRONMENTS_GAUSSIANREWARDMDP_H_
#define INCLUDE_RELE_ENVIRONMENTS_GAUSSIANREWARDMDP_H_

#include "ContinuousMDP.h"

namespace ReLe
{

class GaussianRewardMDP: public ContinuousMDP
{

public:
    GaussianRewardMDP(unsigned int dimension, double mu = 0.0, double sigma = 1.0,
                      double gamma = 0.9, unsigned int horizon = 50);
    GaussianRewardMDP(arma::mat& A, arma::mat& B, arma::vec& mu, arma::mat& sigma,
                      double gamma = 0.9, unsigned int horizon = 50);

    virtual void step(const DenseAction& action, DenseState& nextState,
                      Reward& reward);
    virtual void getInitialState(DenseState& state);

private:
    void initialize(unsigned int dimensions, double mu_s, double sigma_s);

private:
    arma::mat A, B;
    arma::vec mu;
    arma::mat sigma;

};

}



#endif /* INCLUDE_RELE_ENVIRONMENTS_GAUSSIANREWARDMDP_H_ */
