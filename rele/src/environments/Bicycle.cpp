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

#include "rele/environments/Bicycle.h"
#include "rele/utils/RandomGenerator.h"

using namespace std;

namespace ReLe
{

Bicycle::Bicycle(ConfigurationsLabel label) :
    //Sutton's article
    //    DenseMDP(2, 3, 1, false, true)
    //Klein's articles
    DenseMDP(5, 9, 1, false, true, 0.95, 1000), s0type(label)
{
		//state_range = state_range.t();
	}

void Bicycle::step(const FiniteAction& action,
                       DenseState& nextState, Reward& reward)
{

    double T = 2. * ((action.getActionN()/3) /*- 1*/); 				//Torque on handle bars
    double d = 0.02 * ((action.getActionN() % 3) /*- 1*/); 			// Displacement of center of mass (in meters)

    if(noise > 0){
          d += (rand()-0.5)*noise;		 							// Noise between [-0.02, 0.02] meters
    }

    for(int step= 0;step<10;step++){
    	double r_f, r_b, r_CM;
        if (currentState[theta] == 0) // Infinite radius tends to not be handled well
        	r_f = r_b = r_CM = 1.e8;
        else{
        	r_f = l / abs(sin(currentState[theta]));
            r_b = l / abs(tan(currentState[theta]));
            r_CM = sqrt(pow((l - c),2) + (l*l / pow(tan(currentState[theta]),2)));
        }
        double varphi = currentState[omega] + atan(d / h);

        currentState[omega_ddot] = h * M * gravity * sin(varphi);

        currentState[omega_ddot] -= cos(varphi) * (Inertia_dv * sigma_dot * currentState[theta_dot] + copysign(1.0,currentState[theta])*v*v*(M_d * r *(1.0/r_f + 1.0/r_b) + M*h/r_CM));
        currentState[omega_ddot] /= Inertia_bc;

        double theta_ddot = (T - Inertia_dv * sigma_dot * currentState[omega_dot]) / Inertia_dl;

        double df = (delta_time / 10.0); //10.0 = simStep
        currentState[omega_dot] += df * currentState[omega_ddot];
        currentState[omega] += df * currentState[omega_dot];
        currentState[theta_dot] += df * theta_ddot;
        currentState[theta] += df * currentState[theta_dot];

        //Handle bar limits (80 deg.)

        if(currentState[theta]<state_range(3,0))
            currentState[theta] = state_range(3,0);
        if(currentState[theta]>state_range(3,1))
            currentState[theta] = state_range(3,1);

        // Update position (x,y) of tires
        double front_term = psi + currentState[theta] + copysign(1.0,psi + currentState[theta])* asin(v * df / (2.*r_f));
        double back_term = psi + copysign(1.0,psi)*asin(v * df / (2.*r_b));
        x_f += -sin(front_term);
        y_f += cos(front_term);
        x_b += -sin(back_term);
        y_b += cos(back_term);

        // Handle Roundoff errors, to keep the length of the bicycle constant
        double dist = sqrt((x_f-x_b)*(x_f-x_b)+ (y_f-y_b)*(y_f-y_b));
        if (abs(dist - l) > 0.01){
            x_b += (x_b - x_f) * (l - dist)/dist;
            y_b += (y_b - y_f) * (l - dist)/dist;
        }
        // Update psi
        if(x_f==x_b && y_f-y_b < 0)
        	psi = M_PI;
        else if(y_f - y_b > 0)
        	psi = atan((x_b - x_f)/(y_f - y_b));
        else
        	psi = copysign(1.0,x_b - x_f)*(M_PI/2.0) - atan((y_f - y_b)/(x_b-x_f));


    }

    double val = currentState[omega];
    double delta = state_range(0,1);
    if(abs(val) > state_range(0,1)){ // Bicycle fell over
                //return -1.0, True
    	reward = {-1.0};
    	currentState.setAbsorbing();
    	//resetState();
    }else if(isAtGoal()){
    	reward = {reward_goal};
		currentState.setAbsorbing();
    	//resetState();
    }else if (!navigate)
        reward = {reward_shaping};
    else{
    	vec g_l ={goal_loc_x,goal_loc_y};
    	vec x = {x_f-x_b, y_f-y_b};
    	double goal_angle = vector_angle(g_l, x) * M_PI / 180.0;
        reward = {(4.0 - goal_angle*goal_angle) * reward_shaping};

    }
    //reward = {reward[0]*10};
    nextState = currentState;
}


bool Bicycle::isAtGoal(){
        // Anywhere in the goal radius
	    double dist = pow((x_f - goal_loc_x),2) + pow((y_f - goal_loc_y),2);
        if(navigate)
        	return (sqrt(max(0.0,(dist - goal_rsqrd))) < 0.00001);
        else
        	return false;
}

double Bicycle::vector_angle(vec u, vec v){
    return (acos(dot(u,  v)/(norm(u)*norm(v))))*180.0/M_PI;
}

void Bicycle::resetState(){
	x_f = 0;
	y_f = 0;
	x_b = 0;
	y_b = l;
	psi =  atan((y_f-x_f)/(y_b - x_b));
}

void Bicycle::getInitialState(DenseState& state){

	currentState[omega] = 0;
	currentState[omega_dot] = 0;
	currentState[omega_ddot] = 0;
	currentState[theta] = 0;
	currentState[theta_dot] = 0;
	resetState();
    currentState.setAbsorbing(false);
    state=currentState;
}
}

