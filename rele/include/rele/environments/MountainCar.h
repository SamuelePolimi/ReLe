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

#ifndef MOUNTAINCAR_H_
#define MOUNTAINCAR_H_

#include "rele/core/DenseMDP.h"

namespace ReLe
{

/*!
 * This class implements the Mountain Car environment.
 * In this problem a car is placed on a valley between two hills
 * and has to reach the top of one of them.
 * Unfortunately, it is not able to do so by only accelerating and,
 * therefore, it has to accelerate in the opposite direction to use
 * the slope of the other hill to acquire inertia that may allow it
 * to reach the goal.
 * For further information see: https://books.google.it/books?id=2S64-ZZ1fREC&pg=PA479&lpg=PA479&dq=mountain+car+optimization&source=bl&ots=BA_U_Ghnn5&sig=nvBPwQFvIxrE_PEC4TpIFX0FolY&hl=it&sa=X&ved=0ahUKEwij4oys2rTLAhUBoBQKHeQkDzcQ6AEILzAC#v=onepage&q=mountain%20car%20optimization&f=false
 */
class MountainCar: public DenseMDP
{
public:
    enum StateLabel
    {
        velocity = 0, position = 1
    };

    /*!
     * Type of configuration obtained from different experiments
     * in respective articles.
     */
    enum ConfigurationsLabel
    {
        Sutton, Klein, Random
    };

public:

    /*!
     * Constructor.
     * \param label configuration type
     */
    MountainCar(ConfigurationsLabel label = Sutton);

    /*!
     * \see Environment::step
     */
    virtual void step(const FiniteAction& action, DenseState& nextState,
                      Reward& reward) override;

    /*!
     * \see Environment::getInitialState
     */
    virtual void getInitialState(DenseState& state) override;

    //! Configuration type.
    ConfigurationsLabel s0type;
};

}

#endif /* MOUNTAINCAR_H_ */
