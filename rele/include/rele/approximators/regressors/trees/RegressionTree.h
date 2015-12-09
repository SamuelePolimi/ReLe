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

#ifndef INCLUDE_RELE_APPROXIMATORS_REGRESSORS_REGRESSIONTREE_H_
#define INCLUDE_RELE_APPROXIMATORS_REGRESSORS_REGRESSIONTREE_H_

#include "Regressors.h"
#include "nodes/TreeNode.h"
#include "nodes/LeafTreeNode.h"
#include "nodes/InternalTreeNode.h"
#include "nodes/EmptyTreeNode.h"
#include "Features.h"

namespace ReLe
{

template<class InputC, class OutputC>
class RegressionTree: public BatchRegressor_<InputC, OutputC>
{

public:
    RegressionTree(Features_<InputC>& phi,
                   const EmptyTreeNode<OutputC>& emptyNode,
                   unsigned int outputDimensions = 1,
                   unsigned int nMin = 2) :
        BatchRegressor_<InputC, OutputC>(phi, outputDimensions),
        root(nullptr), emptyNode(emptyNode), nMin(nMin), phi(phi)
    {

    }

    virtual OutputC operator() (const InputC& input) override
    {
        if (!root)
            throw std::runtime_error("Empty tree evaluated");

        arma::vec output(this->outputDimension);
        output = root->getValue(phi(input));
        return output;
    }

    /**
     * Set nMin
     * @param nm the minimum number of inputs for splitting
     */
    void setNMin(int nm)
    {
        nMin = nm;
    }

    /**
     * Get nmin
     */
    int getNMin()
    {
        return nMin;
    }

    virtual void trainFeatures(BatchDataFeatures<InputC, OutputC>& featureDataset) override = 0;

    /**
     * Get the root of the tree
     * @return a pointer to the root
     */
    TreeNode<OutputC>* getRoot()
    {
        return root;
    }

    virtual ~RegressionTree()
    {
        cleanTree();
    }

protected:
    void cleanTree()
    {
        if (root && !root->isEmpty())
            delete root;
    }

    void splitDataset(const BatchData<InputC, OutputC>& ds,
                      int cutDir, double cutPoint,
                      arma::uvec& indexesLow,
                      arma::uvec& indexesHigh)
    {
        indexesLow.set_size(ds.size());
        indexesHigh.set_size(ds.size());
        unsigned int lowNumber = 0;
        unsigned int highNumber = 0;

        // split inputs in two subsets
        for (unsigned int i = 0; i < ds.size(); i++)
        {
            auto&& element = phi(ds.getInput(i));
            double tmp = element[cutDir];

            if (tmp < cutPoint)
            {
                indexesLow(lowNumber) = i;
                lowNumber++;
            }
            else
            {
                indexesHigh(highNumber) = i;
                highNumber++;
            }
        }

        indexesLow.resize(lowNumber);
        indexesHigh.resize(highNumber);
    }

    TreeNode<OutputC>* buildLeaf(const BatchData<InputC, OutputC>& ds, LeafType type)
    {
        switch(type)
        {
        case Constant:
            return new LeafTreeNode<InputC, OutputC>(ds);
        case Linear:
            return nullptr; //TODO implement
        case Samples:
            return new SampleLeafTreeNode<InputC, OutputC>(ds.clone());
        default:
            return nullptr;
        }
    }



protected:
    TreeNode<OutputC>* root;
    EmptyTreeNode<OutputC> emptyNode;
    Features_<InputC>& phi;

    unsigned int nMin;  // minimum number of tuples for splitting
};

#define USE_REGRESSION_TREE_MEMBERS               \
	typedef RegressionTree<InputC, OutputC> Base; \
	using Base::phi;                              \
    using Base::root;                             \
    using Base::emptyNode;                        \
    using Base::nMin;

}

#endif /* INCLUDE_RELE_APPROXIMATORS_REGRESSORS_REGRESSIONTREE_H_ */
