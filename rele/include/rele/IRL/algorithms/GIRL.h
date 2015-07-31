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
#include <nlopt.hpp>
#include <cassert>

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
         double gamma, IRLGradType aType, bool useSimplexConstraints = true)
        : policy(policy), data(dataset), rewardf(rewardf),
          gamma(gamma), maxSteps(0), atype(aType), useSimplexConstraints(useSimplexConstraints)
    {
        nbFunEvals = 0;
    }

    virtual ~GIRL() { }

    virtual void run()
    {
        run(arma::vec(), 0);
    }

    virtual void run(arma::vec starting,
                     unsigned int maxFunEvals)
    {
        int dpr = rewardf.getParametersSize();
        assert(dpr > 0);

        if (starting.n_elem == 0)
        {
            starting.ones(dpr);
            starting /= arma::sum(starting);
            //starting.zeros(dpr);
        }
        else
        {
            assert(dpr == starting.n_elem);
        }

        if (maxFunEvals == 0)
            maxFunEvals = std::min(30*dpr, 600);

        nbFunEvals = 0;

        maxSteps = 0;
        int nbEpisodes = data.size();
        for (int i = 0; i < nbEpisodes; ++i)
        {
            int nbSteps = data[i].size();
            if (maxSteps < nbSteps)
                maxSteps = nbSteps;
        }


        //setup optimization algorithm
        nlopt::opt optimizator;

        if(useSimplexConstraints)
        {
            optimizator = nlopt::opt(nlopt::algorithm::LN_COBYLA, dpr);

            std::vector<double> lowerBounds(dpr, 0.0);
            std::vector<double> upperBounds(dpr, 1.0);
            optimizator.set_lower_bounds(lowerBounds);
            optimizator.set_upper_bounds(upperBounds);
            optimizator.add_equality_constraint(GIRL::OneSumConstraint, NULL, 1e-6);

        }
        else
        {
            optimizator = nlopt::opt(nlopt::algorithm::LD_MMA, dpr);
        }

        optimizator.set_min_objective(GIRL::wrapper, this);
        optimizator.set_xtol_rel(1e-8);
        optimizator.set_ftol_rel(1e-8);
        optimizator.set_ftol_abs(1e-8);
        optimizator.set_maxeval(maxFunEvals);

        //optimize dual function
        std::vector<double> parameters(dpr);
        for (int i = 0; i < dpr; ++i)
            parameters[i] = starting[i];
        double minf;
        if (optimizator.optimize(parameters, minf) < 0)
        {
            std::cout << "nlopt failed!" << std::endl;
        }
        else
        {
            std::cout << "found minimum = " << minf << std::endl;

            arma::vec finalP(dpr);
            for(int i = 0; i < dpr; ++i)
            {
                finalP(i) = parameters[i];
            }
            std::cout << std::endl;

            rewardf.setParameters(finalP);
        }
    }

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
    }

    arma::vec ReinforceGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();

        gGradient.zeros(dp,dpr);

        arma::vec sumGradLog(dp), localg;
        arma::vec gradient_J(dp, arma::fill::zeros);
        double Rew;
        arma::mat dRew(1, dpr);

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
        gradient_J /= nbEpisodes;
        gGradient  /= nbEpisodes;

        return gradient_J;
    }

    arma::vec ReinforceBaseGradient(arma::mat& gGradient)
    {
        int dp  = policy.getParametersSize();
        int dpr = rewardf.getParametersSize();
        int nbEpisodes = data.size();

        gGradient.zeros(dp,dpr);

        arma::vec sumGradLog(dp), localg;
        arma::vec gradient_J(dp, arma::fill::zeros);
        double Rew;
        arma::mat dRew(1, dpr);

        arma::vec baseline_J_num(dp, arma::fill::zeros);
        arma::vec baseline_den(dp, arma::fill::zeros);
        arma::vec return_J_ObjEp(nbEpisodes);
        arma::mat sumGradLog_CompEp(dp,nbEpisodes);

        arma::mat baseline_R_num(dp, dpr, arma::fill::zeros);
        std::vector<arma::mat> return_R_ObjEp(nbEpisodes, arma::mat());

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
                Rew += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew += df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();
                // ********************** //

                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

            // *** REINFORCE BASE CORE *** //

            // store the basic elements used to compute the gradients

            return_J_ObjEp(i) = Rew;
            return_R_ObjEp[i] = dRew;

            for (int p = 0; p < dp; ++p)
            {
                sumGradLog_CompEp(p,i) = sumGradLog(p);
            }

            // compute the baselines
            for (int p = 0; p < dp; ++p)
            {
                baseline_J_num(p) += Rew * sumGradLog(p) * sumGradLog(p);
                baseline_den(p) += sumGradLog(p) * sumGradLog(p);
                for (int rp = 0; rp < dpr; ++rp)
                {
                    baseline_R_num(p,rp) += sumGradLog(p) * sumGradLog(p) * dRew(0,rp);
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
                baseline_J = baseline_J_num(p) / baseline_den(p);
            }

            for (int ep = 0; ep < nbEpisodes; ++ep)
            {
                gradient_J[p] += (return_J_ObjEp(ep) - baseline_J) * sumGradLog_CompEp(p,ep);

                for (int rp = 0; rp < dpr; ++rp)
                {
                    double basel = baseline_den(p) != 0 ? baseline_R_num(p,rp) / baseline_den(p) : 0.0;
                    gGradient(p,rp) += (return_R_ObjEp[ep](0,rp) - basel) * sumGradLog_CompEp(p,ep);
                }
            }
        }

        // ********************** //

        // compute mean values
        gradient_J /= nbEpisodes;
        gGradient  /= nbEpisodes;

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
                Rew += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew += df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();

                // compute the gradients
                Rew += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                dRew += df * rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();

                for (int p = 0; p < dp; ++p)
                {
                    gradient_J[p] += df * arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn))) * sumGradLog(p);
                    for (int rp = 0; rp < dpr; ++rp)
                    {
                        gGradient(p,rp) += sumGradLog(p) * dRew(0,rp);
                    }
                }
                // ********************** //

                df *= gamma;

                if (tr.xn.isAbsorbing())
                {
                    assert(nbSteps == t+1);
                    break;
                }
            }

        }
        // compute mean values
        gradient_J /= nbEpisodes;
        gGradient  /= nbEpisodes;

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
        double Rew;
        arma::mat dRew(1, dpr);


        arma::mat baseline_J_num(dp, maxSteps, arma::fill::zeros);
        std::vector<arma::mat> baseline_R_num(maxSteps, arma::mat(dp,dpr, arma::fill::zeros));
        arma::mat baseline_den(dp, maxSteps, arma::fill::zeros);
        arma::mat reward_J_ObjEpStep(nbEpisodes, maxSteps);
        std::vector<std::vector<arma::mat>> reward_R_ObjEpStep(nbEpisodes,
                                         std::vector<arma::mat>(maxSteps,arma::mat()));
        arma::cube sumGradLog_CompEpStep(dp,nbEpisodes, maxSteps);
        arma::vec  maxsteps_Ep(nbEpisodes);

        for (int ep = 0; ep < nbEpisodes; ++ep)
        {
            //core setup
            int nbSteps = data[ep].size();


            // *** GPOMDP CORE *** //
            sumGradLog.zeros();
            double df = 1.0;
            Rew = 0.0;
            dRew.zeros();
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
                double creward = arma::as_scalar(rewardf(vectorize(tr.x, tr.u, tr.xn)));
                arma::mat cdreward = rewardf.diff(vectorize(tr.x, tr.u, tr.xn)).t();
                Rew += df * creward;
                dRew += df * cdreward;
                reward_J_ObjEpStep(ep,t) = df * creward;
                reward_R_ObjEpStep[ep][t] = df * cdreward;


                for (int p = 0; p < dp; ++p)
                {
                    sumGradLog_CompEpStep(p,ep,t) = sumGradLog(p);
                }

                // compute the baselines
                for (int p = 0; p < dp; ++p)
                {
                    baseline_J_num(p,t) += df * creward * sumGradLog(p) * sumGradLog(p);

                    for (int rp = 0; rp < dpr; ++rp)
                    {
                        baseline_R_num[t](p,rp) += df * cdreward(0,rp) * sumGradLog(p) * sumGradLog(p);
                    }
                }

                for (int p = 0; p < dp; ++p)
                {
                    baseline_den(p,t) += sumGradLog(p) * sumGradLog(p);
                }
                // ********************** //

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
                        double basel = baseline_den(p) != 0 ? baseline_R_num[t](p,rp) / baseline_den(p) : 0.0;
                        gGradient(p,rp) += (tmp(0,rp) - basel) * sumGradLog_CompEpStep(p,ep,t);
                    }
                }
            }
        }
        // ************************ //

        // compute mean values
        gradient_J /= nbEpisodes;
        gGradient  /= nbEpisodes;

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

    double objFunction(unsigned int n, const double* x, double* grad)
    {

        ++nbFunEvals;


        arma::vec gradient;
        arma::mat dGradient;
        arma::vec parV(const_cast<double*>(x), n, true);
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

        //        std::cerr << gradient.t();
        //        std::cerr << dGradient;

        if (grad != nullptr)
        {
            arma::vec g = dGradient.t() * gradient;
            //            std::cerr << g.t();
            for (int i = 0; i < g.n_elem; ++i)
            {
                grad[i] = g[i];
            }
        }

        //double norm22 = arma::norm(gradient,2);
        double f = 0.5 * arma::as_scalar(gradient.t()*gradient);
        //        std::cerr << f << std::endl;
        return f;

    }

    static double OneSumConstraint(unsigned int n, const double *x, double *grad, void *data)
    {
        if(grad != nullptr)
        {
            for (unsigned int i = 0; i < n; ++i)
            {
                grad[0] = 1;
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

    static double f1(unsigned int n, const double *x, double *grad, void *data)
    {
        grad = nullptr;
        double val = -1.0;
        for (unsigned int i = 0; i < n; ++i)
        {
            val += x[i];
            if (grad != nullptr)
                grad[i] = 1;
        }
        return val;
    }
    static double f2(unsigned int n, const double *x, double *grad, void *data)
    {
        double val = 1.0;
        for (unsigned int i = 0; i < n; ++i)
        {
            val += -x[i];
            if (grad != nullptr)
                grad[i] = -1;
        }
        return val;
    }

    static double wrapper(unsigned int n, const double* x, double* grad,
                          void* o)
    {
        double value = reinterpret_cast<GIRL*>(o)->objFunction(n, x, grad);

        std::cout << "v= " << value << " ";

        std::cout << "x= ";
        for(int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }

        std::cout << std::endl;

        if(grad)
        {
            std::cout << "g= ";
            for(int i = 0; i < n; i++)
            {
                std::cout << grad[i] << " ";
            }

            std::cout << std::endl;
        }

        return value;
    }

    unsigned int getFunEvals()
    {
        return nbFunEvals;
    }

protected:
    Dataset<ActionC,StateC>& data;
    DifferentiablePolicy<ActionC,StateC>& policy;
    ParametricRegressor& rewardf;
    double gamma;
    unsigned int maxSteps;
    IRLGradType atype;
    unsigned int nbFunEvals;
    bool useSimplexConstraints;
};

} //end namespace


#endif /* GIRL_H_ */
