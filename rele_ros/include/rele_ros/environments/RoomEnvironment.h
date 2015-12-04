/*
 * rele_ros,
 *
 *
 * Copyright (C) 2015 Davide Tateo
 * Versione 1.0
 *
 * This file is part of rele_ros.
 *
 * rele_ros is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rele_ros is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rele_ros.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_RELE_ROS_environmentS_ROOMENVIRONMENT_H_
#define INCLUDE_RELE_ROS_environmentS_ROOMENVIRONMENT_H_

#include "../core/SimulatedEnvironment.h"

namespace ReLe_ROS
{

class SimulatedRoomEnvironment : public SimulatedEnvironment<ReLe::DenseAction, ReLe::DenseState>
{

public:
    SimulatedRoomEnvironment(double controlFrequency);
    virtual ~SimulatedRoomEnvironment();

protected:
    virtual void publishAction(const ReLe::DenseAction& action);
    virtual void setState(ReLe::DenseState& state);
    virtual void setReward(const ReLe::DenseAction& action,
                           const ReLe::DenseState& state, ReLe::Reward& reward);

private:
    void writeSettings();

private:
    int leftMotorHandle;
    int rightMotorHandle;

    int positionHandle;

    int objectiveHandle;

    ros::Publisher motorSpeedPub;


    arma::vec objective;

};

}

#endif /* INCLUDE_RELE_ROS_environmentS_ROOMENVIRONMENT_H_ */