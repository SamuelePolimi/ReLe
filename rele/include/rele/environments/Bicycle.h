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

#ifndef BICYCLE_H_
#define BICYCLE_H_

#include "rele/core/DenseMDP.h"
#include <armadillo>

using namespace std;
using namespace arma;

namespace ReLe
{

class Bicycle: public DenseMDP
{
public:
    enum StateLabel
    {
    	omega = 0, omega_dot = 1, omega_ddot = 2, theta = 3, theta_dot = 4
    };

    enum ConfigurationsLabel
    {
        Sutton, Klein, Random
    };

public:

    Bicycle(ConfigurationsLabel label = Sutton);
    virtual void step(const FiniteAction& action, DenseState& nextState,
                      Reward& reward) override;
    virtual void getInitialState(DenseState& state) override;

    ConfigurationsLabel s0type;
protected:
    bool isAtGoal();

    double vector_angle(mat u, mat v);

    double const reward_fall=-1.o;
    double const reward_fall = -1.0;
    double const reward_goal = 0.01;
    double const goal_rsqrd = 1000.0;
    	//# Square of the radius around the goal (10m)^2
    double navigate = false;
    double const reward_shaping = 0.001;
    double goal_loc_x = 1000;
    double goal_loc_y = 0;
    // Units in Meters and Kilograms
    double const c = 0.66;       	// Horizontal dist between bottom of front wheel and center of mass
    double const d_cm = 0.30;     	// Vertical dist between center of mass and the cyclist
    double const h = 0.94;       	//  Height of the center of mass over the ground
    double const l = 1.11;       	// Dist between front tire and back tire at point on ground
    double const M_c = 15.0;		// Mass of bicycle
    double const M_d = 1.7;      	// Mass of tire
    double const M_p = 60;	        // Mass of cyclist
    double const r = 0.34;		    // Radius of tire
    double const v = 10.0 / 3.6;	// Velocity of bicycle (converted from km/h to m/s)


    // Useful precomputations
    double const M = M_p + M_c;
    double const Inertia_bc = (13.0/3.0) * M_c * h*h + M_p * (h + d_cm)*(h + d_cm);
    double const Inertia_dv = M_d * r*r;
    double const Inertia_dl = 0.5 * M_d * r*r;
    double const sigma_dot = v / r;

    // Simulation Constants
    double const gravity = 9.8;
    double const delta_time = 0.02;
    double const sim_steps = 10;

    double const noise = 0.04;

    double x_f = 0, y_f = 0, x_b = 0, y_b = 0, psi = 0;

    double** state_range = {{-M_PI * 12.0/180.0, M_PI * 12.0/180.0},
                            {-M_PI * 2.0/180.0, M_PI * 2.0/180.0},
                            {-M_PI, M_PI},
                            {-M_PI * 80.0/180.0, M_PI * 80.0/180.0},
                            {-M_PI * 2.0/180.0, M_PI * 2.0/180.0}};

};

}

#endif /* MOUNTAINCAR_H_ */
