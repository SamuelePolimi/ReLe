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

#ifndef GIRL_H_
#define GIRL_H_

#include "IRLAlgorithm.h"
#include "Policy.h"
#include "Transition.h"
#include "ArmadilloExtensions.h"
#include <nlopt.hpp>
#include <cassert>

#include "policy_search/step_rules/StepRules.h"

//#define DISPARITY
//#define LOG_OBJ
//#define SQUARE_OBJ

namespace ReLe
{

enum IRLGradType {R, RB, G, GB, ENAC, NATR, NATRB, NATG, NATGB};

template<class ActionC, class StateC>
class GIRL : public IRLAlgorithm<ActionC, StateC>
{
public:

    GIRL(Dataset<ActionC,StateC>& dataset,
         DifferentiablePolicy<ActionC,StateC>& policy,
         ParametricRegressor& rewardf,
         double gamma, IRLGradType aType,
         bool useSimplexConstraints = true)
        : policy(policy), data(dataset), rewardf(rewardf),
          gamma(gamma), atype(aType),
          useSimplexConstraints(useSimplexConstraints), isRewardLinear(false)
    {
        nbFunEvals = 0;
        maxSteps = data.getEpisodeMaxLenght();
        // initially all features are active
        active_feat.set_size(rewardf.getParametersSize());
        std::iota (std::begin(active_feat), std::end(active_feat), 0);
    }

    GIRL(Dataset<ActionC,StateC>& dataset,
         DifferentiablePolicy<ActionC,StateC>& policy,
         LinearApproximator& rewardf,
         double gamma, IRLGradType aType)
        : GIRL(dataset, policy, rewardf, gamma, aType, true)
    {
        isRewardLinear = true;
    }

    virtual ~GIRL() { }

    //======================================================================
    // RUNNERS
    //----------------------------------------------------------------------
    virtual void run()
    {
        run(arma::vec(), 0);
    }

    virtual void run(arma::vec starting,
                     unsigned int maxFunEvals)
    {
        int dpr = rewardf.getParametersSize();
        assert(dpr > 0);

        if (maxFunEvals == 0)
            maxFunEvals = std::min(30*dpr, 600);

        nbFunEvals = 0;

        // initialize active features set
        preprocessing();

        //compute effective parameters dimension
        int effective_dim = active_feat.n_elem;

        // handle the case of only one active reward feature
        if (effective_dim == 1)
        {
            arma::vec x(dpr, arma::fill::zeros);
            x.elem(active_feat).ones();
            rewardf.setParameters(x);
            return;
        }

        //setup optimization algorithm
        nlopt::opt optimizator;

        if(useSimplexConstraints)
        {
            // simplex constraint reduces the parameter by one element
            --effective_dim;

            // optimizator = nlopt::opt(nlopt::algorithm::LN_COBYLA, effective_dim);
            optimizator = nlopt::opt(nlopt::algorithm::LD_MMA, effective_dim);

            std::vector<double> lowerBounds(effective_dim, 0.0);
            std::vector<double> upperBounds(effective_dim, 1.0);
            optimizator.set_lower_bounds(lowerBounds);
            optimizator.set_upper_bounds(upperBounds);

            // equality constraint
            // optimizator.add_equality_constraint(GIRL::OneSumConstraint, NULL, 1e-3);

            // inequality constraint
            // x >= 0 && sum x <= 1
            std::vector<double> tols(effective_dim + 1, 1e-5);
            optimizator.add_inequality_mconstraint(GIRL::InequalitySimplexConstraints, NULL, tols);
        }
        else
        {
            optimizator = nlopt::opt(nlopt::algorithm::LD_SLSQP, effective_dim);
        }

        optimizator.set_min_objective(GIRL::wrapper, this);
        optimizator.set_xtol_rel(1e-8);
        optimizator.set_ftol_rel(1e-8);
        optimizator.set_ftol_abs(1e-8);
        optimizator.set_maxeval(maxFunEvals);

        // define initial point
        if (starting.n_elem == 0)
        {
            starting.ones(effective_dim);
            starting /= arma::sum(starting);
        }
        else
        {
            assert(effective_dim <= starting.n_elem);
        }
        std::vector<double> parameters(effective_dim);
        for (int i = 0; i < effective_dim; ++i)
            parameters[i] = starting[i];

        //optimize dual function
        double minf;
        if (optimizator.optimize(parameters, minf) < 0)
        {
            std::cout << "nlopt failed!" << std::endl;
        }
        else
        {
            std::cout << "found minimum = " << minf << std::endl;

            // reconstruct parameters
            int dim = active_feat.n_elem;
            int n   = parameters.size();

            arma::vec x(dpr, arma::fill::zeros);

            if (n == dim-1)
            {
                // simplex scenario
                double sumx = 0.0;
                for (int i = 0; i < n; ++i)
                {
                    x(active_feat(i)) = parameters[i];
                    sumx += parameters[i];
                }
                x(active_feat(n)) = 1 - sumx;
            }
            else
            {
                // full features
                for (int i = 0; i < dim; ++i)
                {
                    x(active_feat(i)) = parameters[i];
                }
            }

            rewardf.setParameters(x);
        }
    }


    //======================================================================
    // OPTIMIZATION: WRAPPER AND OBJECTIVE FUNCTION
    //----------------------------------------------------------------------
    static double wrapper(unsigned int n, const double* x, double* grad, void* o)
    {
        arma::vec df;
        arma::vec parV(const_cast<double*>(x), n, true);
        double value = static_cast<GIRL*>(o)->objFunction(parV, df);

        //Save gradient
        if(grad)
        {
            for (int i = 0; i < df.n_elem; ++i)
            {
                grad[i] = df[i];
            }
        }

        //Print gradient and value
        //printOptimizationInfo(value, n, x, grad);

        return value;
    }

    double objFunction(const arma::vec& x, arma::vec& df)
    {

        ++nbFunEvals;


        // reconstruct parameters
        int dpr = rewardf.getParametersSize();
        int n = x.n_elem;
        arma::vec parV(dpr, arma::fill::zeros);
        int dim = active_feat.n_elem;
        if (n == dim-1)
        {
            // simplex scenario
            parV(active_feat(arma::span(0,dim-2))) = x;
            parV(active_feat(dim-1)) = 1.0 - sum(x);
//            double sumx = 0.0;
//            for (int i = 0; i < n; ++i)
//            {
//                parV(active_feat(i)) = x[i];
//                sumx += x[i];
//            }
//            parV(active_feat(n)) = 1 - sumx;
        }
        else
        {
            // full features
            parV(active_feat) = x;
//            for (int i = 0; i < dim; ++i)
//            {
//                parV(active_feat(i)) = x[i];
//            }
        }


        // dispatch the right call
        arma::vec gradient;
        arma::mat dGradient;
        rewardf.setParameters(parV);
        if (atype == IRLGradType::R)
        {
            //            std::cout << "GIRL REINFORCE" << std::endl;
            gradient = ReinforceGradient(dGradient);
        }
        else if (atype == IRLGradType::RB)
        {
            //            std::cout << "GIRL REINFORCE BASE" << std::endl;
            gradient = ReinforceBaseGradient(dGradient);
        }
        else if (atype == IRLGradType::G)
        {
            //            std::cout << "GIRL GPOMDP" << std::endl;
            gradient = GpomdpGradient(dGradient);
        }
        else if (atype == IRLGradType::GB)
        {
            //            std::cout << "GIRL GPOMDP BASE" << std::endl;
            gradient = GpomdpBaseGradient(dGradient);
        }
        else if (atype == IRLGradType::ENAC)
        {
            gradient = ENACGradient(dGradient);
        }
        else if ((atype == IRLGradType::NATR) || (atype == IRLGradType::NATRB) ||
                 (atype == IRLGradType::NATG) || (atype == IRLGradType::NATGB))
        {
            gradient = NaturalGradient(dGradient);
        }
        else
        {
            std::cerr << "GIRL ERROR" << std::endl;
            abort();
        }

        // select only active elements in the derivative of the policy
        // gradient
        if (n == dim-1)
        {
            dGradient = dGradient.cols(active_feat(arma::span(0,dim-2)));
        }
        else if (dpr != dim)
        {
            dGradient = dGradient.cols(active_feat);
        }


        double g2 = arma::as_scalar(gradient.t()*gradient);

#if defined LOG_OBJ || defined SQUARE_OBJ
        arma::vec dJ;
        double J = computeJ(dJ);
        double J4 = std::pow(J, 4);
#elif defined DISPARITY
        arma::vec dD;
        double D = computeDisparity(dD);
        double D2 = D*D;
        std::cout << dD << std::endl;
#endif


#ifdef LOG_OBJ
        double f = std::log(g2) - std::log(J4);
#elif defined SQUARE_OBJ
        double f = g2/J4;
#elif defined DISPARITY
        double f = g2/D2;
#else
        double f = g2;
#endif

        arma::vec dg2 = 2.0*dGradient.t() * gradient;
#if defined LOG_OBJ || defined SQUARE_OBJ
        arma::vec dJ4 = 4.0*dJ*std::pow(J, 3);

#ifdef LOG_OBJ
        df = dg2/g2 - dJ4/J4;
#elif defined SQUARE_OBJ
        df = (dg2*J4 - dJ4*g2) / (J4*J4);
#endif
#elif defined DISPARITY
        arma::vec dD2 = 2.0*dD*D;
        df = (dg2*D2 - dD2*g2) / (D2*D2);
#else
        df = dg2;
#endif

        std::cout << "g2: " << g2 << std::endl;
        std::cout << "f: " << f << std::endl;
        std::cout << "dwdj: " << dGradient;
        std::cout << "df: " << df.t();
        std::cout << "x:  " << x.t();
        std::cout << "x_full:  " << parV.t();
        std::cout << "-----------------------------------------" << std::endl;

        // abort();
        return f;

    }

    //======================================================================
    // METRICS THAT CAN BE INTEGRATED IN THE GIRL OPTIMIZATION
    //----------------------------------------------------------------------
    double computeJ(arma::vec& dR)
    {
        double J = 0;
        dR.zeros(rewardf.getParametersSize());
        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();
            double df = 1.0;

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                auto tr = data[i][t];
                J += df*arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dR += df*rewardf.diff(vectorize(tr.x, tr.u, tr.xn));
                df *= gamma;
            }
        }

        return J/static_cast<double>(nbEpisodes);

    }

    double computeDisparity(arma::vec& dD)
    {
        double D = 0;
        dD.zeros(rewardf.getParametersSize());

        int nEpisodes = data.size();
        for (int i = 0; i < nEpisodes; ++i)
        {
            double R = 0;
            double R2 = 0;
            arma::vec dR(rewardf.getParametersSize(), arma::fill::zeros);
            arma::vec dR2(rewardf.getParametersSize(), arma::fill::zeros);

            double Gamma = 0;

            //Compute episode J
            int nbSteps = data[i].size();
            double df = 1.0;

            for (int t = 0; t < nbSteps; ++t)
            {
                Gamma += df;
                auto tr = data[i][t];
                double r = arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                arma::vec dr = rewardf.diff(vectorize(tr.x, tr.u, tr.xn));
                R += df*r;
                dR += df*dr;
                R2 += df*r*r;
                dR2 += df*dr*r;

                df *= gamma;
            }

            R /= Gamma;
            R2 /= Gamma;

            D += R2 - R*R;
            dD += 2.0*(dR2 - R*dR)/Gamma;

        }

        double N = nEpisodes;

        dD /= N;

        return D/N;
    }

    //======================================================================
    // PRE-PROCESSING
    //----------------------------------------------------------------------

    void preprocessing()
    {
        // get reward parameters dimension
        int dpr = rewardf.getParametersSize();

        //initially all features are active
        active_feat.set_size(dpr);
        std::iota (std::begin(active_feat), std::end(active_feat), 0);

        //if the reward is linear perform preprocessing
        if (isRewardLinear)
        {
            // performs preprocessing in order to remove the features
            // that are constant and the one that are almost never
            // under the given samples
            arma::uvec const_ft;
            arma::vec mu = preproc_linear_reward(const_ft);

            // normalize features
            mu = arma::normalise(mu);

            //find non-zero features
            arma::uvec q = arma::find( abs(mu) > 1e-6);

            //sort indexes
            q = arma::sort(q);
            const_ft = arma::sort(const_ft);

            //compute set difference in order to obtain active features
            auto it = std::set_difference(q.begin(), q.end(), const_ft.begin(), const_ft.end(), active_feat.begin());
            active_feat.resize(it-active_feat.begin());

            std::cout << "--- LINEAR REWARD: PRE-PROCESSING ---" << std::endl;
            std::cout << "Feature expectation\n mu: " << mu.t();
            std::cout << "Constant features\n cf: " << const_ft.t();
            std::cout << "based on mu test, the following features are preserved\n q: " << q.t();
            std::cout << "Finally the active features are\n q - cf: " << active_feat.t();
            std::cout << "--------------------------------------------" << std::endl;

            // force simplex constraint with linear reward parametrizations
            useSimplexConstraints = true;
        }
    }

    /**
     * @brief Compute the feature expectation and identifies the constant features
     * The function computes the features expectation over trajecteries that can be
     * used to remove the features that are never or rarelly visited under the given
     * samples.
     * Moreover, it identifies the features that are constant. We consider a feature
     * constant when its range (max-min) over an episode is less then a threshold.
     * Clearly, this condition must hold for every episode.
     * @param const_features vector storing the indexis of the constant features
     * @param tol threshold used to test the range
     * @return the feature expectation
     */
    arma::vec preproc_linear_reward(arma::uvec& const_features, double tol = 1e-4)
    {
        int nEpisodes = data.size();
        unsigned int dpr = rewardf.getParametersSize();
        unsigned int dp  = policy.getParametersSize();
        arma::vec mu(dpr, arma::fill::zeros);

        arma::mat constant_reward(dpr, nEpisodes, arma::fill::zeros);

        for (int ep = 0; ep < nEpisodes; ++ep)
        {

            //Compute episode J
            int nbSteps = data[ep].size();
            double df = 1.0;


            // store immediate reward over trajectory
            arma::mat reward_vec(dp, nbSteps, arma::fill::zeros);

            for (int t = 0; t < nbSteps; ++t)
            {
                auto tr = data[ep][t];

                reward_vec.col(t) = rewardf.diff(vectorize(tr.x, tr.u, tr.xn));

                mu += df * reward_vec.col(t);

                df *= gamma;
            }

            // check reward range over trajectories
            arma::vec R = range(reward_vec, 1);
            for (int p = 0; p < R.n_elem; ++p)
            {
                // check range along each feature
                if (R(p) <= tol)
                {
                    constant_reward(p,ep) = 1;
                }
            }

        }

        const_features = arma::find( arma::sum(constant_reward,1) == nEpisodes );

        mu /= nEpisodes;

        return mu;

    }

    //======================================================================
    // CONSTRAINTS
    //----------------------------------------------------------------------
    static void InequalitySimplexConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
    {
        result[n] = -1.0;
        for (unsigned int i = 0; i < n; ++i)
        {
            // -x_i <= 0
            result[i] = -x[i];
            // sum x_i - 1 <= 0
            result[n] += x[i];
        }
        //compute the gradient: d c_i / d x_j = g(i*n+j)
        if (grad != nullptr)
        {
            for (unsigned int j = 0; j < n; ++j)
            {
                for (unsigned int i = 0; i < n; ++i)
                {
                    if (i == j)
                    {
                        grad[i*n+j] = -1.0;
                    }
                    else
                    {
                        grad[i*n+j] = 0.0;
                    }

                }
                grad[n*n+j] = 1.0;
            }
        }
    }

    static double OneSumConstraint(unsigned int n, const double *x, double *grad, void *data)
    {
        if(grad != nullptr)
        {
            for (unsigned int i = 0; i < n; ++i)
            {
                grad[i] = 1;
            }
        }

        double val = -1.0;
        //std::cout << "x: ";
        for (unsigned int i = 0; i < n; ++i)
        {
            //std::cout << x[i] << " ";
            val += x[i];
        }
        //std::cout << std::endl << "val: " << val << std::endl;
        return val;
    }


    //======================================================================
    // GETTERS and SETTERS
    //----------------------------------------------------------------------
    virtual arma::vec getWeights()
    {
        return rewardf.getParameters();
    }

    virtual Policy<ActionC, StateC>* getPolicy()
    {
        return &policy;
    }

    void setData(Dataset<ActionC,StateC>& dataset)
    {
        data = dataset;
        maxSteps = data.getEpisodeMaxLenght();
    }

    inline void setIsRewardLinear(bool flag)
    {
        isRewardLinear = flag;
    }

    inline bool getIsRewardLinear()
    {
        return isRewardLinear;
    }

    inline void setUseSimplexConstraints(bool flag)
    {
        useSimplexConstraints = flag;
    }

    inline bool getUseSimplexConstraints()
    {
        return useSimplexConstraints;
    }

    unsigned int getFunEvals()
    {
        return nbFunEvals;
    }


    //======================================================================
    // GRADIENTS
    //----------------------------------------------------------------------
    arma::vec ReinforceGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();

        gGradient.zeros(dp,dpr);

        arma::vec sumGradLog(dp), localg;
        arma::vec gradient_J(dp, arma::fill::zeros);
        double Rew;
        arma::mat dRew(1, dpr);


        int totstep = 0;
        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();


            // *** REINFORCE CORE *** //
            sumGradLog.zeros();
            double df = 1.0;
            Rew = 0.0;
            dRew.zeros();
            // ********************** //

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[i][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** REINFORCE CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                sumGradLog += localg;
                //                std::cout << tr.r[0] << " " << tr.r[1] << std::endl;
                Rew += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew += df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();
                // ********************** //


                ++totstep;
                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

            // *** REINFORCE CORE *** //
            for (int p = 0; p < dp; ++p)
            {
                gradient_J[p] += Rew * sumGradLog(p);
                for (int rp = 0; rp < dpr; ++rp)
                {
                    gGradient(p,rp) += sumGradLog(p) * dRew(0,rp);
                }
            }
            // ********************** //

        }
        // compute mean values
        if (gamma == 1.0)
        {
            gradient_J /= totstep;
            gGradient  /= totstep;
        }
        else
        {
            gradient_J /= nbEpisodes;
            gGradient  /= nbEpisodes;
        }

        return gradient_J;
    }

    arma::vec ReinforceBaseGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();
        int nbEpisodes = data.size();

        // performance (J) gradient
        arma::vec gradient_J(dp, arma::fill::zeros);
        // gradient w.r.t. R weights of performance (J) gradient
        gGradient.zeros(dp,dpr);

        // cumulate reward function over episode
        double Rew;
        // cumulate reward derivative over episode
        arma::mat dRew(1, dpr);

        // sum of log-gradientes and gradient of the policy
        arma::vec sumGradLog(dp), localg;
        // baseline denominator is shared between Rfun and Rder
        arma::vec baseline_den(dp, arma::fill::zeros);
        // the sum of the log-grad is shared between Rfun and Rder
        arma::mat sumGradLog_CompEp(dp,nbEpisodes);

        // variables related to reward function
        arma::vec baseline_Rfun_num(dp, arma::fill::zeros);
        arma::vec return_Rfun_ObjEp(nbEpisodes);

        // variables related to reward derivative
        arma::mat baseline_Rder_num(dp, dpr, arma::fill::zeros);
        std::vector<arma::mat> return_Rder_ObjEp(nbEpisodes, arma::mat());

        int totstep = 0;
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();


            // *** REINFORCE CORE *** //
            sumGradLog.zeros();
            double df = 1.0;
            Rew = 0.0;
            dRew.zeros();
            // ********************** //

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[i][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** REINFORCE CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                sumGradLog += localg;
                Rew  += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew += df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();
                // ********************** //

                ++totstep;
                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

            // *** REINFORCE BASE CORE *** //

            // store the basic elements used to compute the gradients

            return_Rfun_ObjEp(i) = Rew;
            return_Rder_ObjEp[i] = dRew;

            //            for (int p = 0; p < dp; ++p)
            //            {
            //                sumGradLog_CompEp(p,i) = sumGradLog(p);
            //            }

            // compute the baselines
            for (int p = 0; p < dp; ++p)
            {
                //store sum of log-gradients
                sumGradLog_CompEp(p,i) = sumGradLog(p);

                // square sum log-grad
                double tmp = sumGradLog(p) * sumGradLog(p);
                // baseline denominator
                baseline_den(p) += tmp;
                // compute numerator for reward fun (baseline)
                baseline_Rfun_num(p) += Rew * tmp;
                // compute numerator for reward der (baseline)
                for (int rp = 0; rp < dpr; ++rp)
                {
                    baseline_Rder_num(p,rp) += tmp * dRew(0,rp);
                }
            }

            // ********************** //

        }

        // *** REINFORCE BASE CORE *** //

        // compute the gradients
        for (int p = 0; p < dp; ++p)
        {

            double baseline_J = 0;
            if (baseline_den(p) != 0)
            {
                baseline_J = baseline_Rfun_num(p) / baseline_den(p);
            }

            for (int ep = 0; ep < nbEpisodes; ++ep)
            {
                gradient_J[p] += (return_Rfun_ObjEp(ep) - baseline_J) * sumGradLog_CompEp(p,ep);

                for (int rp = 0; rp < dpr; ++rp)
                {
                    double basel = baseline_den(p) != 0 ? baseline_Rder_num(p,rp) / baseline_den(p) : 0.0;
                    gGradient(p,rp) += (return_Rder_ObjEp[ep](0,rp) - basel) * sumGradLog_CompEp(p,ep);
                }
            }
        }
        // in MATLAB the above loops are replaced by
        // dJdtheta = dJdtheta + sumdlogPi .* (ones(dlogpi_r, 1) * sum_rewfun - b);
        // drewdJ = drewdJ + repmat(sumdlogPi,1,size(sum_rewder,2)) .* (repmat(sum_rewder,size(b2,1),1) - b2);

        // ********************** //

        // compute mean values
        if (gamma == 1.0)
        {
            gradient_J /= totstep;
            gGradient  /= totstep;
        }
        else
        {
            gradient_J /= nbEpisodes;
            gGradient  /= nbEpisodes;
        }

        return gradient_J;
    }

    arma::vec GpomdpGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();

        gGradient.zeros(dp,dpr);

        arma::vec sumGradLog(dp), localg;
        arma::vec gradient_J(dp, arma::fill::zeros);
        double Rew;
        arma::mat dRew(1, dpr);

        int totstep = 0;
        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();


            // *** GPOMDP CORE *** //
            sumGradLog.zeros();
            double df = 1.0;
            Rew = 0.0;
            dRew.zeros();
            // ********************** //

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[i][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** GPOMDP CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                sumGradLog += localg;

                // compute the reward gradients
                Rew  = df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew = df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();

                for (int p = 0; p < dp; ++p)
                {
                    gradient_J[p] += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn))) * sumGradLog(p);
                    for (int rp = 0; rp < dpr; ++rp)
                    {
                        gGradient(p,rp) += sumGradLog(p) * dRew(0,rp);
                    }
                }
                // ********************** //

                ++totstep;
                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

        }
        // compute mean values
        if (gamma == 1.0)
        {
            gradient_J /= totstep;
            gGradient  /= totstep;
        }
        else
        {
            gradient_J /= nbEpisodes;
            gGradient  /= nbEpisodes;
        }

        return gradient_J;
    }

    arma::vec GpomdpBaseGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();
        int nbEpisodes = data.size();

        gGradient.zeros(dp,dpr);

        arma::vec sumGradLog(dp), localg;
        arma::vec gradient_J(dp, arma::fill::zeros);


        // baseline denominator is shared between rew fun and rew der
        arma::mat baseline_den(dp, maxSteps, arma::fill::zeros);

        // baseline numerator of reward function is a (dp x maxStep) matrix
        arma::mat baseline_J_num(dp, maxSteps, arma::fill::zeros);
        // for each episode it store the discounted immediate reward
        arma::mat reward_J_ObjEpStep(nbEpisodes, maxSteps);

        // variables related to the reward derivative
        // baseline numerator is a vector of maxSteps elements where each element is a (dp x dpr) matrix
        std::vector<arma::mat> baseline_R_num(maxSteps, arma::mat(dp,dpr, arma::fill::zeros));
        //for each episode, for each step, it store the discounted reward derivative (1 x dpr)
        std::vector<std::vector<arma::mat>> reward_R_ObjEpStep(nbEpisodes,
                                         std::vector<arma::mat>(maxSteps,arma::mat())); //(1 x dpr)

        // store for each episode the sum of log-grad at a given time step
        arma::cube sumGradLog_CompEpStep(dp,nbEpisodes, maxSteps);

        // store the length of each episode
        arma::vec  maxsteps_Ep(nbEpisodes);

        int totstep = 0;
        for (int ep = 0; ep < nbEpisodes; ++ep)
        {
            //core setup
            int nbSteps = data[ep].size();


            // *** GPOMDP CORE *** //
            sumGradLog.zeros();
            double df = 1.0;
            // ********************** //

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[ep][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** GPOMDP CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                sumGradLog += localg;

                // store the basic elements used to compute the gradients
                double creward = df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                arma::mat cdreward = df *rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t(); //(1 x dpr)
                reward_J_ObjEpStep(ep,t) =  creward;
                reward_R_ObjEpStep[ep][t] = cdreward;


                for (int p = 0; p < dp; ++p)
                {
                    sumGradLog_CompEpStep(p,ep,t) = sumGradLog(p);
                }

                // compute the baselines
                for (int p = 0; p < dp; ++p)
                {
                    double tmp = sumGradLog(p) * sumGradLog(p);

                    baseline_J_num(p,t) += creward * tmp;

                    for (int rp = 0; rp < dpr; ++rp)
                    {
                        baseline_R_num[t](p,rp) += cdreward(0,rp) * tmp;
                    }

                    baseline_den(p,t) += tmp;
                }

                //for (int p = 0; p < dp; ++p)
                //{
                //    baseline_den(p,t) += sumGradLog(p) * sumGradLog(p);
                //}
                // ********************** //

                ++totstep;
                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

            // store the actual length of the current episode (<= maxsteps)
            maxsteps_Ep(ep) = nbSteps;

        }

        // *** GPOMDP BASE CORE *** //

        //for (int i = 0; i < baseline_R_num.size(); ++i)
        //{
        //    std::cout << baseline_R_num[i] / arma::repmat(baseline_den.col(i),1,baseline_R_num[i].n_cols) << std::endl;
        //}

        // compute the gradients
        for (int p = 0; p < dp; ++p)
        {
            for (int ep = 0; ep < nbEpisodes; ++ep)
            {
                for (int t = 0; t < maxsteps_Ep(ep); ++t)
                {

                    double baseline_J = 0;
                    if (baseline_den(p,t) != 0)
                    {
                        baseline_J = baseline_J_num(p,t) / baseline_den(p,t);
                    }

                    gradient_J[p] += (reward_J_ObjEpStep(ep,t) - baseline_J) * sumGradLog_CompEpStep(p,ep,t);

                    arma::mat& tmp = reward_R_ObjEpStep[ep][t];
                    for (int rp = 0; rp < dpr; ++rp)
                    {
                        double basel = baseline_den(p) != 0 ? baseline_R_num[t](p,rp) / baseline_den(p,t) : 0.0;
                        gGradient(p,rp) += (tmp(0,rp) - basel) * sumGradLog_CompEpStep(p,ep,t);
                    }
                }
            }
        }
        // ************************ //

        // compute mean values
        if (gamma == 1.0)
        {
            gradient_J /= totstep;
            gGradient  /= totstep;
        }
        else
        {
            gradient_J /= nbEpisodes;
            gGradient  /= nbEpisodes;
        }

        return gradient_J;
    }

    arma::vec ENACGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();

        gGradient.zeros(dp,dpr);

        arma::vec localg;
        double Rew;
        arma::vec g(dp+1, arma::fill::zeros), phi(dp+1);
        arma::mat fisher(dp+1,dp+1, arma::fill::zeros);

        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();


            // *** eNAC CORE *** //
            double df = 1.0;
            Rew = 0.0;
            phi.zeros();
            //    #ifdef AUGMENTED
            phi(dp) = 1.0;
            //    #endif
            // ********************** //

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[i][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** eNAC CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                double creward = arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                Rew += df * creward;

                //Construct basis functions
                for (unsigned int i = 0; i < dp; ++i)
                    phi[i] += df * localg[i];
                // ********************** //

                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

            fisher += phi * phi.t();
            g += Rew * phi;

        }


        arma::vec nat_grad;
        int rnk = arma::rank(fisher);
        //        std::cout << rnk << " " << fisher << std::endl;
        if (rnk == fisher.n_rows)
        {
            nat_grad = arma::solve(fisher, g);
        }
        else
        {
            std::cerr << "WARNING: Fisher Matrix is lower rank (rank = " << rnk << ")!!! Should be " << fisher.n_rows << std::endl;

            arma::mat H = arma::pinv(fisher);
            nat_grad = H * g;
        }

        return nat_grad.rows(0,dp-1);
    }

    arma::vec NaturalGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();
        arma::vec localg;
        arma::mat fisher(dp,dp, arma::fill::zeros);

        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            //core setup
            int nbSteps = data[i].size();

            //iterate the episode
            for (int t = 0; t < nbSteps; ++t)
            {
                Transition<ActionC, StateC>& tr = data[i][t];
                //            std::cout << tr.x << " " << tr.u << " " << tr.xn << " " << tr.r[0] << std::endl;

                // *** eNAC CORE *** //
                localg = policy.difflog(tr.x, tr.u);
                fisher += localg * localg.t();
                // ********************** //

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

        }
        fisher /= nbEpisodes;

        arma::vec gradient;
        if (atype == IRLGradType::NATR)
        {
            gradient = ReinforceGradient(gGradient);
        }
        else if (atype == IRLGradType::NATRB)
        {
            gradient = ReinforceBaseGradient(gGradient);
        }
        else if (atype == IRLGradType::NATG)
        {
            gradient = GpomdpGradient(gGradient);
        }
        else if (atype == IRLGradType::NATGB)
        {
            gradient = GpomdpBaseGradient(gGradient);
        }


        gGradient.zeros(dp,dpr);

        arma::vec nat_grad;
        int rnk = arma::rank(fisher);
        //        std::cout << rnk << " " << fisher << std::endl;
        if (rnk == fisher.n_rows)
        {
            nat_grad = arma::solve(fisher, gradient);
        }
        else
        {
            std::cerr << "WARNING: Fisher Matrix is lower rank (rank = " << rnk << ")!!! Should be " << fisher.n_rows << std::endl;

            arma::mat H = arma::pinv(fisher);
            nat_grad = H * gradient;
        }

        return nat_grad;
    }

protected:
    void printOptimizationInfo(double value, unsigned int n, const double* x,
                               double* grad)
    {
        std::cout << "v= " << value << " ";
        std::cout << "x= ";
        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
        if (grad)
        {
            std::cout << "g= ";
            for (int i = 0; i < n; i++)
            {
                std::cout << grad[i] << " ";
            }
            std::cout << std::endl;
        }
    }

protected:
    Dataset<ActionC,StateC>& data;
    DifferentiablePolicy<ActionC,StateC>& policy;
    ParametricRegressor& rewardf;
    double gamma;
    unsigned int maxSteps;
    IRLGradType atype;
    unsigned int nbFunEvals;
    bool useSimplexConstraints, isRewardLinear;
    arma::uvec active_feat;
};


} //end namespace


#endif /* GIRL_H_ */
