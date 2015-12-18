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

/*
 * Written by: Carlo D'Eramo
 */

#include "MAB/InternetAds.h"
#include "MAB/Roulette.h"
#include "td/DoubleQ-Learning.h"
#include "nonparametric/SequentialPolicy.h"
#include "Core.h"

#include <iostream>

using namespace std;
using namespace ReLe;


/*
 * MAB test with InternetAds or Roulette. Sequential policy is used
 * with the purpose to execute each possible action sequentially until
 * the end of the episode.
 */

enum EnvironmentLabel
{
    iAds, R
};

int main(int argc, char *argv[])
{
    EnvironmentLabel e = EnvironmentLabel::R;

    if(e == iAds)
    {
        InternetAds mab(10, 0.95, InternetAds::Second);

        unsigned int episodeLength = 10;

        SequentialPolicy policy(mab.getSettings().finiteActionDim, episodeLength);
        //Q_Learning agent(policy);
        DoubleQ_Learning agent(policy);

        auto&& core = buildCore(mab, agent);
        core.getSettings().episodeLength = episodeLength;
        for(unsigned int i = 0; i < 20000; i++)
        {
            cout << endl << "### Starting episode " << i << " ###" << endl;
            core.runEpisode();
        }
    }
    else
    {
        Roulette mab(0.95);

        unsigned int episodeLength = mab.getSettings().finiteActionDim;

        SequentialPolicy policy(mab.getSettings().finiteActionDim, episodeLength);
        Q_Learning agent(policy);
        //DoubleQ_Learning agent(policy);

        auto&& core = buildCore(mab, agent);
        core.getSettings().episodeLength = episodeLength;
        for(unsigned int i = 1; i <= 2000; i++)
        {
            agent.setAlpha(1.0 / i);
            cout << endl << "### Starting episode " << i << " ###" << endl;
            core.runEpisode();
        }
    }
}
