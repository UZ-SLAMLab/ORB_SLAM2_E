// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Modified Raul Mur Artal (2014)
// - Stop criterium (solve function)

#include "optimization_algorithm_levenberg.h"

#include <iostream>

#include "../stuff/timeutil.h"

#include "sparse_optimizer.h"
#include "solver.h"
#include "batch_stats.h"


using namespace std;

namespace g2o {

  OptimizationAlgorithmLevenberg::OptimizationAlgorithmLevenberg(Solver* solver) :
    OptimizationAlgorithmWithHessian(solver)
  {
    _currentLambda = -1.;
    _tau = 1e-5;
    _goodStepUpperScale = 2./3.;
    _goodStepLowerScale = 1./3.;
    _userLambdaInit = _properties.makeProperty<Property<double> >("initialLambda", 0.);
    _maxTrialsAfterFailure = _properties.makeProperty<Property<int> >("maxTrialsAfterFailure", 10);
    _ni=2.;
    _levenbergIterations = 0;
    _nBad = 0;
  }

  OptimizationAlgorithmLevenberg::~OptimizationAlgorithmLevenberg()
  {
  }

OptimizationAlgorithm::SolverResult OptimizationAlgorithmLevenberg::solve(int iteration, bool online)
{
    assert(_optimizer && "_optimizer not set");
    assert(_solver->optimizer() == _optimizer && "underlying linear solver operates on different graph");

    if (iteration == 0 && !online)
    { // built up the CCS structure, here due to easy time measure
        bool ok = _solver->buildStructure();
        if (! ok)
        {
            cerr << __PRETTY_FUNCTION__ << ": Failure while building CCS structure" << endl;
            return OptimizationAlgorithm::Fail;
        }
    }

    double t=get_monotonic_time();

    if (bInFEA)
    {
        //vector<vector<float> > vPoints = GetPointCoordinates(pFEA2->vVertices);
        pFEA2->setbfea(true);
        _optimizer->computeActiveErrors();
        pFEA2->setbfea(false);
    }
    else
        _optimizer->computeActiveErrors();

    G2OBatchStatistics* globalStats = G2OBatchStatistics::globalStats();
    if (globalStats)
    {
        globalStats->timeResiduals = get_monotonic_time()-t;
        t=get_monotonic_time();
    }

    double currentChi = _optimizer->activeRobustChi2();
    double tempChi=currentChi;
    double iniChi = currentChi;

    if (bInFEA) if (pFEA2->bDebugMode) cout << "           - iniChi before FEA = " << iniChi << endl;

    _solver->buildSystem();
    if (globalStats)
    {
        globalStats->timeQuadraticForm = get_monotonic_time()-t;
    }

    // core part of the Levenberg algorithm
    if (iteration == 0)
    {
      _currentLambda = computeLambdaInit();
      _ni = 2;
      _nBad = 0;
    }

    double rho=0;
    int& qmax = _levenbergIterations;
    qmax = 0;
    do
    {
        _optimizer->push();
        if (globalStats)
        {
            globalStats->levenbergIterations++;
            t=get_monotonic_time();
        }
        // update the diagonal of the system matrix
        _solver->setLambda(_currentLambda, true);
        bool ok2 = _solver->solve();
        if (globalStats)
        {
            globalStats->timeLinearSolution+=get_monotonic_time()-t;
            t=get_monotonic_time();
        }
        _optimizer->update(_solver->x());
        if (globalStats)
        {
            globalStats->timeUpdate = get_monotonic_time()-t;
        }

        // restore the diagonal
        _solver->restoreDiagonal();

        if (bInFEA)
        {
            pFEA2->setbfea2(true);
            _optimizer->computeActiveErrors();  // Computes reprojection errors and stores them in _error
            pFEA2->setbfea2(false);
        }
        else
            _optimizer->computeActiveErrors();  // Computes reprojection errors and stores them in _error

        tempChi = _optimizer->activeRobustChi2();   // Sums of all the reprojection errors

        if (!ok2)
            tempChi=std::numeric_limits<double>::max();

        if (bInFEA)
        {
            float sE = 0.0;
            float nsE = 0.0;

            vector<vector<float> > vPoints = GetPointCoordinates(pFEA2->vVertices);

            pFEA2->Set_uf(vPoints);
            pFEA2->ComputeDisplacement();
            pFEA2->ComputeForces();

            sE = pFEA2->ComputeStrainEnergy();
            nsE = pFEA2->NormalizeStrainEnergy();

            /*
            for (unsigned int a=0; a<pFEA2->K.size(); a++){
                for (unsigned int b=0; b<pFEA2->K[a].size(); b++){
                    cout << pFEA2->vvf[a][b] << " ";
                }
                cout << endl;
            }
            cout << endl;
            */

            float w_rE = 1.0;
            float w_sE = 5.0;
            if (qmax==0)
            {
                w_rE = 1.0;
                w_sE = 2.0;
                if (pFEA2->bDebugMode) cout << "             currentChi1 / sE / nsE / fFactorFEA = " << currentChi << " / " << sE << " / " << nsE << " / " << fFactorFEA << endl;
                currentChi += nsE;
                if (pFEA2->bDebugMode) cout << "             currentChi2 = " << currentChi << endl;
            }
            else
                if (pFEA2->bDebugMode) cout << "             currentChi = " << currentChi << endl;

            if (pFEA2->bDebugMode) cout << "             wOPT = tempChi = w_rE·rE + w_sE·sE = " << w_rE << "·" << tempChi << " + " << w_sE << "·" << nsE << " = " << (w_rE*tempChi + w_sE*nsE) << endl;
            tempChi = w_rE*tempChi + w_sE*nsE;
        }

        rho = (currentChi-tempChi);         // Difference between reprojection errors before and after the optimization
        double scale = computeScale();
        scale += 1e-3; // make sure it's non-zero :)
        if (bInFEA) if (pFEA2->bDebugMode) cout << "             rho = (currentChi - tempChi)/scale = (" << currentChi << " - " << tempChi << ")/" << scale << " = " << rho/scale << endl;
        rho /= scale;

        if (rho>0 && g2o_isfinite(tempChi)) // last step was good
        {
            double alpha = 1.-pow((2*rho-1),3);
            // crop lambda between minimum and maximum factors
            alpha = (std::min)(alpha, _goodStepUpperScale);
            double scaleFactor = (std::max)(_goodStepLowerScale, alpha);
            _currentLambda *= scaleFactor;
            _ni = 2;
            currentChi=tempChi;
            _optimizer->discardTop();
        }
        else
        {
            _currentLambda*=_ni;
            _ni*=2;
            _optimizer->pop(); // restore the last state before trying to optimize
        }
        qmax++;
    }
    while (rho<0 && qmax < _maxTrialsAfterFailure->value() && ! _optimizer->terminate());

    if (qmax == _maxTrialsAfterFailure->value() || rho==0)
        return Terminate;

    //Stop criteria (Raul)
    if((iniChi-currentChi)*1e3<iniChi)
        _nBad++;
    else
        _nBad=0;

    if(_nBad>=3)
        return Terminate;

    return OK;
}

  double OptimizationAlgorithmLevenberg::computeLambdaInit() const
  {
    if (_userLambdaInit->value() > 0)
      return _userLambdaInit->value();
    double maxDiagonal=0.;
    for (size_t k = 0; k < _optimizer->indexMapping().size(); k++) {
      OptimizableGraph::Vertex* v = _optimizer->indexMapping()[k];
      assert(v);
      int dim = v->dimension();
      for (int j = 0; j < dim; ++j){
        maxDiagonal = std::max(fabs(v->hessian(j,j)),maxDiagonal);
      }
    }
    return _tau*maxDiagonal;
  }

double OptimizationAlgorithmLevenberg::computeScale() const
{
    double scale = 0.;
    for (size_t j=0; j < _solver->vectorSize(); j++)
    {
        scale += _solver->x()[j] * (_currentLambda * _solver->x()[j] + _solver->b()[j]);
        //cout << j << " " << _solver->x()[j] << " " << _solver->b()[j] << endl;
    }
    return scale;
}

  void OptimizationAlgorithmLevenberg::setMaxTrialsAfterFailure(int max_trials)
  {
    _maxTrialsAfterFailure->setValue(max_trials);
  }

  void OptimizationAlgorithmLevenberg::setUserLambdaInit(double lambda)
  {
    _userLambdaInit->setValue(lambda);
  }

  void OptimizationAlgorithmLevenberg::printVerbose(std::ostream& os) const
  {
    os
      << "\t schur= " << _solver->schur()
      << "\t lambda= " << FIXED(_currentLambda)
      << "\t levenbergIter= " << _levenbergIterations;
  }

  void OptimizationAlgorithmLevenberg::setPtrFea(FEA2* pFeaInput)
  {
      pFEA2 = pFeaInput;
  }

  vector<vector<float> > OptimizationAlgorithmLevenberg::GetPointCoordinates(vector<g2o::VertexSBAPointXYZ*> vVertices)
  {
        vector<vector<float> > vPoints;

        for (unsigned int i=0; i<pFEA2->vVertices.size(); i++)
        {
            VertexSBAPointXYZ* v2 = vVertices[i];
            Eigen::Matrix<double, 3, 1> mPoint = v2->estimate();

            vector<float> vPoint;
            vPoint.push_back(mPoint[0]);
            vPoint.push_back(mPoint[1]);
            vPoint.push_back(mPoint[2]);

            vPoints.push_back(vPoint);
        }

        return vPoints;
  }

} // end namespace
