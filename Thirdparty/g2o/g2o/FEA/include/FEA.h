/**
* This file is an addon to ORB-SLAM2.
*
* Its goal is to create a mesh from the 3D points within the computed map.
* Also, it will assembly a rigidity matrix for such mesh, and compute the
* strain deformation energy.
*
* Copyright (C) 2017 Íñigo Cirauqui Viloria <icirauquiviloria at gmail dot com> (University of Zaragoza)
*
* This addon, as ORB-SLAM2, is free software: you can redistribute it and/or
* modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your
* option) any later version.
*
* This addon is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/


/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FEA_H
#define FEA_H

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <thread>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

// G20 LIBS
#include "../../types/types_seven_dof_expmap.h"


using namespace std;

class FEA
{
public:	// FUNCTIONS

	//Constructor
    //FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, bool bSetDebug);
    FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug);

    //Destructor
    ~FEA();

    vector<vector<float> > ComputeKei(vector<vector<float> > vfPts);

    bool MatrixAssembly();

    vector<vector<float> > InvertMatrixEigen(vector<vector<float> > m1);
    vector<vector<float> > MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2);

    void Set_u0();

    void Set_uf(vector<vector<float> > vPoints);

    void ComputeDisplacement();

    void ComputeForces();

    float ComputeStrainEnergy();
    float NormalizeStrainEnergy();

    float GetStrainEnergy();
    float GetNormalizedStrainEnergy();

    void setbfea(bool bSet);
    void setbfea2(bool bSet);
//    void setCurrEdge(int input);

public:	// VARIABLES

    vector<g2o::VertexSBAPointXYZ*> vVertices;

    // Frame
    unsigned int nFrameId;

    // Young Modulus [Pa] & Poisson Coefficient
    unsigned int E;
    float nu;

    // Lamé parameters & Behaviour matrix
    float lambda = 0.0;
    float G = 0.0;
    vector<vector<float> > D;

    // Element depth
    float h = 0.0;

    // Gauss Points
    float fg = 0.5773502692;
    vector<vector<float> > gs;

    //MapPoint coordinates and normals
    vector<vector<float> > vMPsXYZN_t;
    vector<vector<float> > vMPsXYZN_t2;
    vector<vector<float> > vMPsXYZN_u;
    vector<vector<float> > vMPsXYZN_u2;

    unsigned int Ksize;// = 0;
    vector<vector<float> > K;
    float DetK = 0.0;
    vector<vector<float> > K1;
    float DetK1 = 0.0;

    vector<vector<int> > triangles_t;
    vector<vector<int> > triangles_u;
    vector<vector<int> > quads_t;
    vector<vector<int> > quads_u;

    vector<vector<int> > vIdxNewVertices;
    vector<vector<int> > vIdxNewVertices_u;
    vector<vector<int> > vNewPointsBase;
    vector<vector<int> > vNewPointsBase_u;

    // Displacement
    vector<float> u0;  // Starting position of mappoints
    vector<float> uf;  // New position after being modified by the optimization
    vector<float> a;// = vector<float>(Ksize,0.0);     // VectorAssembly(u)

    // Forces
    vector<float> f;// = vector<float>(Ksize,0.0);

    // Strain energy
    float sE;
    float nsE;
    float StartingSE;
    float CurrentSE;
    float TempSE;
    float LastSE = 0.0;

    // Normalization vector
    float fNormFactor;

    bool bInFEA = false;
    bool bInFEA2 = false;
//    int nCurrEdge = 0;

    bool it0 = false;

    // Debug mode
    bool bDebugMode = false;
}; // class FEA

#endif // FEA_H
