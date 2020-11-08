/**
* This file is an addon to ORB-SLAM2.
*
* Its goal is to create a mesh from the 3D points within the computed map.
* Also, it will assembly a rigidity matrix for such mesh, and compute the
* strain deformation energy.
*
* Copyright (C) 2017 Íñigo Cirauqui Viloria <icirauqui at gmail dot com> (University of Zaragoza)
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

#include "../include/FEA.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

#include <iostream>
#include <fstream>

#include <thread>

// Invert matrix
#include <cstdlib>
#include <iomanip>


//FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, bool bSetDebug)
FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg1(input_fg1), fg2(-input_fg1), fg(input_fg1), bDebugMode(bSetDebug)
{
    Ksize = 0.0;
    sE = 0.0;

    // Lamé Parameters
    lambda = ( nu*E ) / ( (1+nu) * (1-2*nu) );
    G = E / ( 2 * (1+nu) );

    float D_temp[6][6] = { (lambda+2*G) ,    lambda    ,    lambda    ,     0.0     ,     0.0     ,     0.0     ,
                              lambda    , (lambda+2*G) ,    lambda    ,     0.0     ,     0.0     ,     0.0     ,
                              lambda    ,    lambda    , (lambda+2*G) ,     0.0     ,     0.0     ,     0.0     ,
                                0.0     ,      0.0     ,      0.0     ,      G      ,     0.0     ,     0.0     ,
                                0.0     ,      0.0     ,      0.0     ,     0.0     ,      G      ,     0.0     ,
                                0.0     ,      0.0     ,      0.0     ,     0.0     ,     0.0     ,      G      };

    for (unsigned int i=0; i<6; i++)
        for (unsigned int j=0; j<6; j++)
            D[i][j] = D_temp[i][j];

    vector<float> g0;     g0.push_back(+fg);  g0.push_back(+fg);  g0.push_back(+fg);     gs.push_back(g0);
    vector<float> g1;     g1.push_back(+fg);  g1.push_back(+fg);  g1.push_back(-fg);     gs.push_back(g1);
    vector<float> g2;     g2.push_back(+fg);  g2.push_back(-fg);  g2.push_back(+fg);     gs.push_back(g2);
    vector<float> g3;     g3.push_back(+fg);  g3.push_back(-fg);  g3.push_back(-fg);     gs.push_back(g3);
    vector<float> g4;     g4.push_back(-fg);  g4.push_back(+fg);  g4.push_back(+fg);     gs.push_back(g4);
    vector<float> g5;     g5.push_back(-fg);  g5.push_back(+fg);  g5.push_back(-fg);     gs.push_back(g5);
    vector<float> g6;     g6.push_back(-fg);  g6.push_back(-fg);  g6.push_back(+fg);     gs.push_back(g6);
    vector<float> g7;     g7.push_back(-fg);  g7.push_back(-fg);  g7.push_back(-fg);     gs.push_back(g7);
}


FEA::~FEA()
{
}



vector<vector<float> > FEA::ComputeKei(vector<vector<float> > vfPts)
{
    vector<vector<float> > vBtDB;
    vector<float> vBtDBi = vector<float>(24,0.0);
        for (unsigned int i=0; i<24; i++)
            vBtDB.push_back(vBtDBi);

    for (unsigned int ops=0; ops<gs.size(); ops++)
    {
        float xi = gs[ops][0];
        float eta = gs[ops][1];
        float zeta = gs[ops][2];

        float dN1dxi = -0.125 * ( (1-eta) * (1-zeta) );
        float dN2dxi = +0.125 * ( (1-eta) * (1-zeta) );
        float dN3dxi = +0.125 * ( (1+eta) * (1-zeta) );
        float dN4dxi = -0.125 * ( (1+eta) * (1-zeta) );
        float dN5dxi = -0.125 * ( (1-eta) * (1+zeta) );
        float dN6dxi = +0.125 * ( (1-eta) * (1+zeta) );
        float dN7dxi = +0.125 * ( (1+eta) * (1+zeta) );
        float dN8dxi = -0.125 * ( (1+eta) * (1+zeta) );

        float dN1deta = -0.125 * ( (1-xi) * (1-zeta) );
        float dN2deta = -0.125 * ( (1+xi) * (1-zeta) );
        float dN3deta = +0.125 * ( (1+xi) * (1-zeta) );
        float dN4deta = +0.125 * ( (1-xi) * (1-zeta) );
        float dN5deta = -0.125 * ( (1-xi) * (1+zeta) );
        float dN6deta = -0.125 * ( (1+xi) * (1+zeta) );
        float dN7deta = +0.125 * ( (1+xi) * (1+zeta) );
        float dN8deta = +0.125 * ( (1-xi) * (1+zeta) );

        float dN1dzeta = -0.125 * ( (1-xi) * (1-eta) );
        float dN2dzeta = -0.125 * ( (1+xi) * (1-eta) );
        float dN3dzeta = -0.125 * ( (1+xi) * (1+eta) );
        float dN4dzeta = -0.125 * ( (1-xi) * (1+eta) );
        float dN5dzeta = +0.125 * ( (1-xi) * (1-eta) );
        float dN6dzeta = +0.125 * ( (1+xi) * (1-eta) );
        float dN7dzeta = +0.125 * ( (1+xi) * (1+eta) );
        float dN8dzeta = +0.125 * ( (1-xi) * (1+eta) );

        /*dxdxi*/   float J_00 = dN1dxi*vfPts[0][0] + dN2dxi*vfPts[1][0] + dN3dxi*vfPts[2][0] + dN4dxi*vfPts[3][0] + dN5dxi*vfPts[4][0] + dN6dxi*vfPts[5][0] + dN7dxi*vfPts[6][0] + dN8dxi*vfPts[7][0];
        /*dydxi*/   float J_01 = dN1dxi*vfPts[0][1] + dN2dxi*vfPts[1][1] + dN3dxi*vfPts[2][1] + dN4dxi*vfPts[3][1] + dN5dxi*vfPts[4][1] + dN6dxi*vfPts[5][1] + dN7dxi*vfPts[6][1] + dN8dxi*vfPts[7][1];
        /*dzdxi*/   float J_02 = dN1dxi*vfPts[0][2] + dN2dxi*vfPts[1][2] + dN3dxi*vfPts[2][2] + dN4dxi*vfPts[3][2] + dN5dxi*vfPts[4][2] + dN6dxi*vfPts[5][2] + dN7dxi*vfPts[6][2] + dN8dxi*vfPts[7][2];
        /*dxdeta*/  float J_10 = dN1deta*vfPts[0][0] + dN2deta*vfPts[1][0] + dN3deta*vfPts[2][0] + dN4deta*vfPts[3][0] + dN5deta*vfPts[4][0] + dN6deta*vfPts[5][0] + dN7deta*vfPts[6][0] + dN8deta*vfPts[7][0];
        /*dydeta*/  float J_11 = dN1deta*vfPts[0][1] + dN2deta*vfPts[1][1] + dN3deta*vfPts[2][1] + dN4deta*vfPts[3][1] + dN5deta*vfPts[4][0] + dN6deta*vfPts[5][0] + dN7deta*vfPts[6][0] + dN8deta*vfPts[7][0];
        /*dzdeta*/  float J_12 = dN1deta*vfPts[0][2] + dN2deta*vfPts[1][2] + dN3deta*vfPts[2][2] + dN4deta*vfPts[3][2] + dN5deta*vfPts[4][0] + dN6deta*vfPts[5][0] + dN7deta*vfPts[6][0] + dN8deta*vfPts[7][0];
        /*dxdzeta*/ float J_20 = dN1dzeta*vfPts[0][0] + dN2dzeta*vfPts[1][0] + dN3dzeta*vfPts[2][0] + dN4dzeta*vfPts[3][0] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0] + dN7dzeta*vfPts[6][0] + dN8dzeta*vfPts[7][0];
        /*dydzeta*/ float J_21 = dN1dzeta*vfPts[0][1] + dN2dzeta*vfPts[1][1] + dN3dzeta*vfPts[2][1] + dN4dzeta*vfPts[3][1] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0] + dN7dzeta*vfPts[6][0] + dN8dzeta*vfPts[7][0];
        /*dzdzeta*/ float J_22 = dN1dzeta*vfPts[0][2] + dN2dzeta*vfPts[1][2] + dN3dzeta*vfPts[2][2] + dN4dzeta*vfPts[3][2] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0] + dN7dzeta*vfPts[6][0] + dN8dzeta*vfPts[7][0];

        float Jac = J_00*J_11*J_22 + J_01*J_12*J_20 + J_10*J_21*J_02 - J_20*J_11*J_02 - J_10*J_01*J_22 - J_21*J_12*J_00;

        //float J_00 = dxdxi;     float J_01 = dydxi;     float J_02 = dzdxi;
        //float J_10 = dxdeta;    float J_11 = dydeta;    float J_12 = dzdeta;
        //float J_20 = dxdzeta;   float J_21 = dydzeta;   float J_22 = dzdzeta;

        float J1_00 = (+1) * ( (J_11*J_22) - (J_21*J_12) ) / Jac;     float J1_01 = (-1) * ( (J_01*J_22) - (J_21*J_02) ) / Jac;     float J1_02 = (-1) * ( (J_01*J_12) - (J_11*J_02) ) / Jac;
        float J1_10 = (-1) * ( (J_10*J_22) - (J_20*J_12) ) / Jac;     float J1_11 = (-1) * ( (J_00*J_22) - (J_20*J_02) ) / Jac;     float J1_12 = (-1) * ( (J_00*J_12) - (J_10*J_02) ) / Jac;
        float J1_20 = (+1) * ( (J_10*J_21) - (J_20*J_11) ) / Jac;     float J1_21 = (-1) * ( (J_00*J_21) - (J_20*J_01) ) / Jac;     float J1_22 = (-1) * ( (J_00*J_11) - (J_10*J_01) ) / Jac;

        //float Jac1 = J1_00*J1_11*J1_22 + J1_01*J1_12*J1_20 + J1_10*J1_21*J1_02 - J1_20*J1_11*J1_02 - J1_10*J1_01*J1_22 - J1_21*J1_12*J1_00;

        float dN1dx = J1_00*dN1dxi + J1_01*dN1deta + J1_02*dN1dzeta;     float dN1dy = J1_10*dN1dxi + J1_11*dN1deta + J1_12*dN1dzeta;     float dN1dz = J1_20*dN1dxi + J1_21*dN1deta + J1_22*dN1dzeta;
        float dN2dx = J1_00*dN2dxi + J1_01*dN2deta + J1_02*dN2dzeta;     float dN2dy = J1_10*dN2dxi + J1_11*dN2deta + J1_12*dN2dzeta;     float dN2dz = J1_20*dN2dxi + J1_21*dN2deta + J1_22*dN2dzeta;
        float dN3dx = J1_00*dN3dxi + J1_01*dN3deta + J1_02*dN3dzeta;     float dN3dy = J1_10*dN3dxi + J1_11*dN3deta + J1_12*dN3dzeta;     float dN3dz = J1_20*dN3dxi + J1_21*dN3deta + J1_22*dN3dzeta;
        float dN4dx = J1_00*dN4dxi + J1_01*dN4deta + J1_02*dN4dzeta;     float dN4dy = J1_10*dN4dxi + J1_11*dN4deta + J1_12*dN4dzeta;     float dN4dz = J1_20*dN4dxi + J1_21*dN4deta + J1_22*dN4dzeta;
        float dN5dx = J1_00*dN5dxi + J1_01*dN5deta + J1_02*dN5dzeta;     float dN5dy = J1_10*dN5dxi + J1_11*dN5deta + J1_12*dN5dzeta;     float dN5dz = J1_20*dN5dxi + J1_21*dN5deta + J1_22*dN5dzeta;
        float dN6dx = J1_00*dN6dxi + J1_01*dN6deta + J1_02*dN6dzeta;     float dN6dy = J1_10*dN6dxi + J1_11*dN6deta + J1_12*dN6dzeta;     float dN6dz = J1_20*dN6dxi + J1_21*dN6deta + J1_22*dN6dzeta;
        float dN7dx = J1_00*dN7dxi + J1_01*dN7deta + J1_02*dN7dzeta;     float dN7dy = J1_10*dN7dxi + J1_11*dN7deta + J1_12*dN7dzeta;     float dN7dz = J1_20*dN7dxi + J1_21*dN7deta + J1_22*dN7dzeta;
        float dN8dx = J1_00*dN8dxi + J1_01*dN8deta + J1_02*dN8dzeta;     float dN8dy = J1_10*dN8dxi + J1_11*dN8deta + J1_12*dN8dzeta;     float dN8dz = J1_20*dN8dxi + J1_21*dN8deta + J1_22*dN8dzeta;

        float B[6][24] = {  dN1dx ,  0.0  ,  0.0  , dN2dx ,  0.0  ,  0.0  , dN3dx ,  0.0  ,  0.0  , dN4dx ,  0.0  ,  0.0  , dN5dx ,  0.0  ,  0.0  , dN6dx ,  0.0  ,  0.0  , dN7dx ,  0.0  ,  0.0  , dN8dx ,  0.0  ,  0.0  ,
                             0.0  , dN1dy ,  0.0  ,  0.0  , dN2dy ,  0.0  ,  0.0  , dN3dy ,  0.0  ,  0.0  , dN4dy ,  0.0  ,  0.0  , dN5dy ,  0.0  ,  0.0  , dN6dy ,  0.0  ,  0.0  , dN7dy ,  0.0  ,  0.0  , dN8dy ,  0.0  ,
                             0.0  ,  0.0  , dN1dz ,  0.0  ,  0.0  , dN2dz ,  0.0  ,  0.0  , dN3dz ,  0.0  ,  0.0  , dN4dz ,  0.0  ,  0.0  , dN5dz ,  0.0  ,  0.0  , dN6dz ,  0.0  ,  0.0  , dN7dz ,  0.0  ,  0.0  , dN8dz ,
                            dN1dy , dN1dx ,  0.0  , dN2dy , dN2dx ,  0.0  , dN3dy , dN3dx ,  0.0  , dN4dy , dN4dx ,  0.0  , dN5dy , dN5dx ,  0.0  , dN6dy , dN6dx ,  0.0  , dN7dy , dN7dx ,  0.0  , dN8dy , dN8dx ,  0.0  ,
                            dN1dz ,  0.0  , dN1dx , dN2dz ,  0.0  , dN2dx , dN3dz ,  0.0  , dN3dx , dN4dz ,  0.0  , dN4dx , dN5dz ,  0.0  , dN5dx , dN6dz ,  0.0  , dN6dx , dN7dz ,  0.0  , dN7dx , dN8dz ,  0.0  , dN8dx ,
                             0.0  , dN1dz , dN1dy ,  0.0  , dN2dz , dN2dy ,  0.0  , dN3dz , dN3dy ,  0.0  , dN4dz , dN4dy ,  0.0  , dN5dz , dN5dy ,  0.0  , dN6dz , dN6dy ,  0.0  , dN7dz , dN7dy ,  0.0  , dN8dz , dN8dy  };
       //                     0       1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      23

        float BtD[24][6] = {0.0};

        BtD[0][0] =  B[0][0]*D[0][0];
        BtD[0][1] =  B[0][0]*D[0][1];
        BtD[0][2] =  B[0][0]*D[0][2];
        BtD[0][3] =  B[3][0]*D[3][3];
        BtD[0][4] =  B[4][0]*D[4][4];
        //BtD[0][5] =  0.0;

        BtD[1][0] =  B[1][1]*D[1][0];
        BtD[1][1] =  B[1][1]*D[1][1];
        BtD[1][2] =  B[1][1]*D[1][2];
        BtD[1][3] =  B[3][1]*D[3][3];
        //BtD[1][4] =  0.0;
        BtD[1][5] =  B[5][1]*D[5][5];

        BtD[2][0] =  B[2][2]*D[2][0];
        BtD[2][1] =  B[2][2]*D[2][1];
        BtD[2][2] =  B[2][2]*D[2][2];
        //BtD[2][3] =  0.0;
        BtD[2][4] =  B[4][2]*D[4][4];
        BtD[2][5] =  B[5][2]*D[5][5];

        BtD[3][0] =  B[0][3]*D[0][0];
        BtD[3][1] =  B[0][3]*D[0][1];
        BtD[3][2] =  B[0][3]*D[0][2];
        BtD[3][3] =  B[3][3]*D[3][3];
        BtD[3][4] =  B[4][3]*D[4][4];
        //BtD[3][5] =  0.0;

        BtD[4][0] =  B[1][4]*D[1][0];
        BtD[4][1] =  B[1][4]*D[1][1];
        BtD[4][2] =  B[1][4]*D[1][2];
        BtD[4][3] =  B[3][4]*D[3][3];
        //BtD[4][4] =  0.0;
        BtD[4][5] =  B[5][4]*D[5][5];

        BtD[5][0] =  B[2][5]*D[2][0];
        BtD[5][1] =  B[2][5]*D[2][1];
        BtD[5][2] =  B[2][5]*D[2][2];
        //BtD[5][3] =  0.0;
        BtD[5][4] =  B[4][5]*D[4][4];
        BtD[5][5] =  B[5][5]*D[5][5];

        BtD[6][0] =  B[0][6]*D[0][0];
        BtD[6][1] =  B[0][6]*D[0][1];
        BtD[6][2] =  B[0][6]*D[0][2];
        BtD[6][3] =  B[3][6]*D[3][3];
        BtD[6][4] =  B[4][6]*D[4][4];
        //BtD[6][5] =  0.0;

        BtD[7][0] =  B[1][7]*D[1][0];
        BtD[7][1] =  B[1][7]*D[1][1];
        BtD[7][2] =  B[1][7]*D[1][2];
        BtD[7][3] =  B[3][7]*D[3][3];
        //BtD[7][4] =  0.0;
        BtD[7][5] =  B[5][7]*D[5][5];

        BtD[8][0] =  B[2][8]*D[2][0];
        BtD[8][1] =  B[2][8]*D[2][1];
        BtD[8][2] =  B[2][8]*D[2][2];
        //BtD[8][3] =  0.0;
        BtD[8][4] =  B[4][8]*D[4][4];
        BtD[8][5] =  B[5][8]*D[5][5];

        BtD[9][0] =  B[0][9]*D[0][0];
        BtD[9][1] =  B[0][9]*D[0][1];
        BtD[9][2] =  B[0][9]*D[0][2];
        BtD[9][3] =  B[3][9]*D[3][3];
        BtD[9][4] =  B[4][9]*D[4][4];
        //BtD[9][5] =  0.0;

        BtD[10][0] = B[1][10]*D[1][0];
        BtD[10][1] = B[1][10]*D[1][1];
        BtD[10][2] = B[1][10]*D[1][2];
        BtD[10][3] = B[3][10]*D[3][3];
        //BtD[10][4] = 0.0;
        BtD[10][5] = B[5][10]*D[5][5];

        BtD[11][0] = B[2][11]*D[2][0];
        BtD[11][1] = B[2][11]*D[2][1];
        BtD[11][2] = B[2][11]*D[2][2];
        //BtD[11][3] = 0.0;
        BtD[11][4] = B[4][11]*D[4][4];
        BtD[11][5] = B[5][11]*D[5][5];

        BtD[12][0] = B[0][12]*D[0][0];
        BtD[12][1] = B[0][12]*D[0][1];
        BtD[12][2] = B[0][12]*D[0][2];
        BtD[12][3] = B[3][12]*D[3][3];
        BtD[12][4] = B[4][12]*D[4][4];
        //BtD[12][5] = 0.0;

        BtD[13][0] = B[1][13]*D[1][0];
        BtD[13][1] = B[1][13]*D[1][1];
        BtD[13][2] = B[1][13]*D[1][2];
        BtD[13][3] = B[3][13]*D[3][3];
        //BtD[13][4] = 0.0;
        BtD[13][5] = B[5][13]*D[5][5];

        BtD[14][0] = B[2][14]*D[2][0];
        BtD[14][1] = B[2][14]*D[2][1];
        BtD[14][2] = B[2][14]*D[2][2];
        //BtD[14][3] = 0.0;
        BtD[14][4] = B[4][14]*D[4][4];
        BtD[14][5] = B[5][14]*D[5][5];

        BtD[15][0] = B[0][15]*D[0][0];
        BtD[15][1] = B[0][15]*D[0][1];
        BtD[15][2] = B[0][15]*D[0][2];
        BtD[15][3] = B[3][15]*D[3][3];
        BtD[15][4] = B[4][15]*D[4][4];
        //BtD[15][5] = 0.0;

        BtD[16][0] = B[1][16]*D[1][0];
        BtD[16][1] = B[1][16]*D[1][1];
        BtD[16][2] = B[1][16]*D[1][2];
        BtD[16][3] = B[3][16]*D[3][3];
        //BtD[16][4] = 0.0;
        BtD[16][5] = B[5][16]*D[5][5];

        BtD[17][0] = B[2][17]*D[2][0];
        BtD[17][1] = B[2][17]*D[2][1];
        BtD[17][2] = B[2][17]*D[2][2];
        //BtD[17][3] = 0.0;
        BtD[17][4] = B[4][17]*D[4][4];
        BtD[17][5] = B[5][17]*D[5][5];

        BtD[18][0] = B[0][18]*D[0][0];
        BtD[18][1] = B[0][18]*D[0][1];
        BtD[18][2] = B[0][18]*D[0][2];
        BtD[18][3] = B[3][18]*D[3][3];
        BtD[18][4] = B[4][18]*D[4][4];
        //BtD[18][5] = 0.0;

        BtD[19][0] = B[1][19]*D[1][0];
        BtD[19][1] = B[1][19]*D[1][1];
        BtD[19][2] = B[1][19]*D[1][2];
        BtD[19][3] = B[3][19]*D[3][3];
        //BtD[19][4] = 0.0;
        BtD[19][5] = B[5][19]*D[5][5];

        BtD[20][0] = B[2][20]*D[2][0];
        BtD[20][1] = B[2][20]*D[2][1];
        BtD[20][2] = B[2][20]*D[2][2];
        //BtD[20][3] = 0.0;
        BtD[20][4] = B[4][20]*D[4][4];
        BtD[20][5] = B[5][20]*D[5][5];

        BtD[21][0] = B[0][21]*D[0][0];
        BtD[21][1] = B[0][21]*D[0][1];
        BtD[21][2] = B[0][21]*D[0][2];
        BtD[21][3] = B[3][21]*D[3][3];
        BtD[21][4] = B[4][21]*D[4][4];
        //BtD[21][5] = 0.0;

        BtD[22][0] = B[1][22]*D[1][0];
        BtD[22][1] = B[1][22]*D[1][1];
        BtD[22][2] = B[1][22]*D[1][2];
        BtD[22][3] = B[3][22]*D[3][3];
        //BtD[22][4] = 0.0;
        BtD[22][5] = B[5][22]*D[5][5];

        BtD[23][0] = B[2][23]*D[2][0];
        BtD[23][1] = B[2][23]*D[2][1];
        BtD[23][2] = B[2][23]*D[2][2];
        //BtD[23][3] = 0.0;
        BtD[23][4] = B[4][23]*D[4][4];
        BtD[23][5] = B[5][23]*D[5][5];

        for (unsigned int i=0; i<24; i++)
            for (unsigned int j=0; j<24; j++)
                vBtDB[i][j] += BtD[i][0]*B[0][j] + BtD[i][1]*B[1][j] + BtD[i][2]*B[2][j] + BtD[i][3]*B[3][j] + BtD[i][4]*B[4][j] + BtD[i][5]*B[5][j];
    }

    return vBtDB;
}



bool FEA::MatrixAssembly()
{
    // To easily insert the nodes into the stiffness matrix, the node indexes are used to
    // set their position in the matrix. NULL positions will be set to 0 by default.
    // The size of the stiffness matrix equals the size of the MPs vector, the positions
    // corresponding to a offvertex are set to zero in the stiffness matrix
    int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    Ksize = 3*nTotalNodes;
    cout << "                - MatAssembly (hex, 24 DoF per el. Curr: " << quads_t.size() << " quads, " << nTotalNodes << " nodes, Ksize = " << Ksize << ")" << endl;
    // DoF = nEl*nodesperelement*DOFnode

    // Make sure we have points before continue.
    if (Ksize<=3)
        return false;

    K.clear();
    vector<float> ki = vector<float>(Ksize,0.0);
    for (unsigned int i=0; i<Ksize; i++)
        K.push_back(ki);

    float fTimer = 0.0;
    float freqKei = cv::getTickFrequency();
    float ft1 = 0.0;
    float ft2 = 0.0;

    for (unsigned int i=0; i<quads_t.size(); i++) {
        vector<int> nodes;
        nodes.push_back(quads_t[i][0]);
        nodes.push_back(quads_t[i][1]);
        nodes.push_back(quads_t[i][2]);
        nodes.push_back(quads_t[i][3]);
        nodes.push_back(quads_t[i][0] + vMPsXYZN_t.size());
        nodes.push_back(quads_t[i][1] + vMPsXYZN_t.size());
        nodes.push_back(quads_t[i][2] + vMPsXYZN_t.size());
        nodes.push_back(quads_t[i][3] + vMPsXYZN_t.size());

        vector<vector<float> > vfPts;
        for (unsigned int j=0; j<4; j++) {
            vector<float> vfPtsi;
            vfPtsi.push_back(vMPsXYZN_t[nodes[j]][0]);
            vfPtsi.push_back(vMPsXYZN_t[nodes[j]][1]);
            vfPtsi.push_back(vMPsXYZN_t[nodes[j]][2]);
            vfPts.push_back(vfPtsi);
        }

        for (unsigned int j=0; j<4; j++) {
            vector<float> vfPtsi;
            vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][0]);
            vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][1]);
            vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][2]);
            vfPts.push_back(vfPtsi);
        }

        vector<vector<float> > vBtDB = ComputeKei(vfPts);

        vector<int> matpos;
        for (unsigned j=0; j<nodes.size(); j++)
            matpos.push_back(nodes[j]*3);

        for (unsigned int node1=0; node1<nodes.size(); node1++)
            for (unsigned int node2=0; node2<nodes.size(); node2++)
                for (unsigned int m=0; m<3; m++)
                    for (unsigned int n=0; n<3; n++)
                        K[matpos[node1]+m][matpos[node2]+n] += vBtDB[node1*3+m][node2*3+n];
    }

    return true;
}



bool FEA::ComputeK1() {
    // Método Gauss Jordan para cálculo de matriz inversa (A|I) -> (I|A1)

    cout << "                - Compute K inverse " << endl;

    vector<vector<float> > Ktemp = K;

    vector<vector<int> > nBodies;
    int nTmp = 2;
    for (unsigned int i=0; i<Ktemp.size(); i++) {
        if (nTmp==2) {
            if (Ktemp[i][i]==0) {
                nTmp = 0;
                vector<int> nBody;
                nBody.push_back(nTmp);
                nBody.push_back(i);
                nBodies.push_back(nBody);
            }
            else {
                nTmp = 1;
                vector<int> nBody;
                nBody.push_back(nTmp);
                nBody.push_back(i);
                nBodies.push_back(nBody);
            }
        }
        else {
            if (Ktemp[i][i]==0 && nTmp!=0) {
                nTmp = 0;
                nBodies[nBodies.size()-1].push_back(i-1);
                vector<int> nBody;
                nBody.push_back(nTmp);
                nBody.push_back(i);
                nBodies.push_back(nBody);
            }
            if (Ktemp[i][i]!=0 && nTmp!=1) {
                nTmp = 1;
                nBodies[nBodies.size()-1].push_back(i-1);
                vector<int> nBody;
                nBody.push_back(nTmp);
                nBody.push_back(i);
                nBodies.push_back(nBody);
            }
            if (i==Ktemp.size()-1)
                nBodies[nBodies.size()-1].push_back(i-1);
        }
    }

    vector<int> nSizeBodies;
    int nIdBiggerBody = -1;
    int nBiggerBodySize = -1;
    for (unsigned int i=0; i<nBodies.size(); i++) {
        int nBodySize = nBodies[i][2]-nBodies[i][1];
        nSizeBodies.push_back(nBodySize);
        if (nBodySize>nBiggerBodySize) {
            nBiggerBodySize = nBodySize;
            nIdBiggerBody = i;
        }
    }

    vector<float> fBodySizeRatios;
    vector<bool> bBodiesToCompute = vector<bool>(nSizeBodies.size(), true);
    for (unsigned int i=0; i<nSizeBodies.size(); i++) {
        float fBodyRatio = nBiggerBodySize/nSizeBodies[i];
        fBodySizeRatios.push_back(fBodyRatio);
        if (fBodyRatio>5)
            bBodiesToCompute[i] = false;
    }

    K1.clear();
    K1 = vector<vector<float> >(Ksize,vector<float>(Ksize,(0.0)));
    for (unsigned int i=0; i<K1.size(); i++)
        K1[i][i] = 1.0;

    for (unsigned int i=0; i<Ksize; i++) {
        float fPivot = Ktemp[i][i];

        if (fPivot==0) {
            unsigned int nRowObj = i;
            if (i<(Ksize-1))
                for (unsigned int j=i+1; j<Ksize; j++) {
                    if (Ktemp[j][i]!=0) {
                        nRowObj = j;
                        break;
                    }
                }

            if (nRowObj==i) {
                cout << "                  Fail! No row to switch! Process stopped at row " << i << endl;
                return false;
            }

            vector<int> vSwitch;
            vSwitch.push_back(i);
            vSwitch.push_back(nRowObj);
            vnSwitchedRows.push_back(vSwitch);

            vector<float> vfAux = Ktemp[i];

            for (unsigned int aa=0; aa<Ksize; aa++) {
                Ktemp[i][aa] = Ktemp[nRowObj][aa];
                Ktemp[nRowObj][aa] = vfAux[aa];
            }
        }

        fPivot = Ktemp[i][i];

        if (fPivot!=1)
            for (unsigned int j=0; j<Ksize; j++) {
                if (Ktemp[i][j]!=0)
                    Ktemp[i][j] /= fPivot;
                if (K1[i][j]!=0)
                    K1[i][j] /= fPivot;
            }


        for (unsigned int j=0; j<Ksize; j++) {
            if (j==i)
                continue;
            float fPivotRow = Ktemp[j][i];
            for (unsigned int s=0; s<Ksize; s++) {
                if (Ktemp[i][s]!=0)
                    Ktemp[j][s] -= fPivotRow*Ktemp[i][s];
                if (K1[i][s]!=0)
                    K1[j][s] -= fPivotRow*K1[i][s];
            }
        }
    }

    cout << "                  K1 computed and saved " << endl;

    return true;
}



/*
void FEA::Set_uf(vector<vector<float> > vPoints)
{
    vector<vector<float> > top_layer;
    //top_layer.resize(u0.size());
    top_layer.resize(vMPsXYZN_t);

    for (unsigned int i=0; i<vPoints.size(); i++)
    {
        vector<float> vfMP_top;

        float fPos1 = vPoints[i][0];
        float fPos2 = vPoints[i][1];
        float fPos3 = vPoints[i][2];

        vfMP_top.push_back(fPos1);
        vfMP_top.push_back(fPos2);
        vfMP_top.push_back(fPos3);

        top_layer[i] = vfMP_top;
    }


    for (unsigned int i=0; i<triangles_t.size(); i++)
    {
        // Get vertices and indices for reference triangle
        int v0idx = triangles_t[i][0];
        vector<float> v0 = top_layer[v0idx];

        int v1idx = triangles_t[i][1];
        vector<float> v1 = top_layer[v1idx];

        int v2idx = triangles_t[i][2];
        vector<float> v2 = top_layer[v2idx];

        // Set variables to store the solution
        vector<float> m01;
        int m01idx = vIdxNewVertices[i][0];
        m01.push_back((v0[0]+v1[0])/2);
        m01.push_back((v0[1]+v1[1])/2);
        m01.push_back((v0[2]+v1[2])/2);
        top_layer[m01idx] = m01;

        vector<float> m02;
        int m02idx = vIdxNewVertices[i][1];
        m02.push_back((v0[0]+v2[0])/2);
        m02.push_back((v0[1]+v2[1])/2);
        m02.push_back((v0[2]+v2[2])/2);
        top_layer[m02idx] = m02;

        vector<float> m12;
        int m12idx = vIdxNewVertices[i][2];
        m12.push_back((v1[0]+v2[0])/2);
        m12.push_back((v1[1]+v2[1])/2);
        m12.push_back((v1[2]+v2[2])/2);
        top_layer[m12idx] = m12;

        vector<float> g;
        int gidx = vIdxNewVertices[i][3];
        g.push_back((v0[0]+v1[0]+v2[0])/3);
		g.push_back((v0[1]+v1[1]+v2[1])/3);
		g.push_back((v0[2]+v1[2]+v2[2])/3);
		top_layer[gidx] = g;
    }

    cout << "size top_layer = " << top_layer.size() << endl;

    vMPsXYZN_t.clear();
    vMPsXYZN_t = top_layer;

    vMPsXYZN_t2.clear();
    vMPsXYZN_t2.resize(top_layer.size());

    // Set second layer
    vector<vector<float> > bottom_layer;
    bottom_layer.resize(vMPsXYZN_t2.size());

    for (unsigned int i=0; i<top_layer.size(); i++)
	{
	    if (top_layer[i].empty())
            continue;

		vector<float> point;
		float pos1 = top_layer[i][0] - h;
		float pos2 = top_layer[i][1] - h;
		float pos3 = top_layer[i][2] - h;
		point.push_back(pos1);
		point.push_back(pos2);
		point.push_back(pos3);

		bottom_layer[i] = point;
	}

	vMPsXYZN_t2 = bottom_layer;





    uf.clear();
    uf.resize(u0);

    for (unsigned int i=0; i<top_layer.size(); i++)
    {
        if (top_layer[i].empty())
        {
            uf.push_back(0.0);
            uf.push_back(0.0);
            uf.push_back(0.0);
            continue;
        }

        for (unsigned int j=0; j<3; j++)
            uf.push_back(top_layer[i][j]);
    }

    for (unsigned int i=0; i<bottom_layer.size(); i++)
    {
        if (bottom_layer[i].empty())
        {
            uf.push_back(0.0);
            uf.push_back(0.0);
            uf.push_back(0.0);
            continue;
        }

        for (unsigned int j=0; j<3; j++)
            uf.push_back(bottom_layer[i][j]);
    }
}
*/



void FEA::Set_u0()
{
    u0.clear();

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
        for (unsigned int j=0; j<3; j++)
            u0.push_back(vMPsXYZN_t[i][j]);

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
        for (unsigned int j=0; j<3; j++)
            u0.push_back(vMPsXYZN_t2[i][j]);
}


/*
void FEA::Set_uf(vector<vector<float> > vPoints)
{
    //vector<vector<float> > top_layer;
    //top_layer.resize(u0.size());
    //top_layer.resize(vMPsXYZN_t);

    //int size_vMPsXYZN_t = vMPsXYZN_t.size();
    vMPsXYZN_t.clear();
    vMPsXYZN_t.resize(vPoints.size());

    for (unsigned int i=0; i<vPoints.size(); i++)
    {
        vector<float> vfMP_top;

        float fPos1 = vPoints[i][0];
        float fPos2 = vPoints[i][1];
        float fPos3 = vPoints[i][2];

        vfMP_top.push_back(fPos1);
        vfMP_top.push_back(fPos2);
        vfMP_top.push_back(fPos3);

        vMPsXYZN_t[i] = vfMP_top;
    }

    //cout << "vNewPointsBase.size() = " << vNewPointsBase.size() << endl;


    for (unsigned int i=0; i<vNewPointsBase.size(); i++)
    {
        vector<float> mi;

        if (vNewPointsBase[i].size()==2)
        {
            int v0idx = vNewPointsBase[i][0];
            int v1idx = vNewPointsBase[i][1];

            vector<float> v0 = vMPsXYZN_t[v0idx];
            vector<float> v1 = vMPsXYZN_t[v1idx];

            mi.push_back((v0[0]+v1[0])/2);
            mi.push_back((v0[1]+v1[1])/2);
            mi.push_back((v0[2]+v1[2])/2);
        }
        if(vNewPointsBase[i].size()==3)
        {
            int v0idx = vNewPointsBase[i][0];
            int v1idx = vNewPointsBase[i][1];
            int v2idx = vNewPointsBase[i][2];

            vector<float> v0 = vMPsXYZN_t[v0idx];
            vector<float> v1 = vMPsXYZN_t[v1idx];
            vector<float> v2 = vMPsXYZN_t[v2idx];

            mi.push_back((v0[0]+v1[0]+v2[0])/3);
            mi.push_back((v0[1]+v1[1]+v2[1])/3);
            mi.push_back((v0[2]+v1[2]+v2[2])/3);
        }

        vMPsXYZN_t.push_back(mi);
    }

    //cout << "size top_layer(vMPsXYZN_t) = " << vMPsXYZN_t.size() << endl;

    vMPsXYZN_t2.clear();
    vMPsXYZN_t2.resize(vMPsXYZN_t.size());

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
	{
	    if (vMPsXYZN_t[i].empty())
            continue;

		vector<float> point;
		float pos1 = vMPsXYZN_t[i][0] - h;
		float pos2 = vMPsXYZN_t[i][1] - h;
		float pos3 = vMPsXYZN_t[i][2] - h;
		point.push_back(pos1);
		point.push_back(pos2);
		point.push_back(pos3);

		vMPsXYZN_t2[i] = point;
	}

	//cout << "size bottom_layer(vMPsXYZN_t2) = " << vMPsXYZN_t2.size() << endl;

    uf.clear();
    //uf.resize(u0);

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
    {
        if (vMPsXYZN_t[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t[i][j]);
    }

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
    {
        if (vMPsXYZN_t2[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t2[i][j]);
    }

    //cout << "u0 size = " << u0.size() << endl;
    //cout << "uf size = " << uf.size() << endl;
}
*/

void FEA::Set_uf(vector<vector<float> > vPoints)
{
    vMPsXYZN_t.clear();
    vMPsXYZN_t.resize(vPoints.size());

    for (unsigned int i=0; i<vPoints.size(); i++)
    {
        vector<float> vfMP;

        vfMP.push_back(vPoints[i][0]);
        vfMP.push_back(vPoints[i][1]);
        vfMP.push_back(vPoints[i][2]);

        vMPsXYZN_t[i] = vfMP;
    }

    //cout << "vNewPointsBase.size() = " << vNewPointsBase.size() << endl;

    for (unsigned int i=0; i<vNewPointsBase.size(); i++)
    {
        vector<float> mi;

        if (vNewPointsBase[i].size()==2)
        {
            int v0idx = vNewPointsBase[i][0];
            int v1idx = vNewPointsBase[i][1];

            vector<float> v0 = vMPsXYZN_t[v0idx];
            vector<float> v1 = vMPsXYZN_t[v1idx];

            mi.push_back((v0[0]+v1[0])/2);
            mi.push_back((v0[1]+v1[1])/2);
            mi.push_back((v0[2]+v1[2])/2);
        }
        if(vNewPointsBase[i].size()==3)
        {
            int v0idx = vNewPointsBase[i][0];
            int v1idx = vNewPointsBase[i][1];
            int v2idx = vNewPointsBase[i][2];

            vector<float> v0 = vMPsXYZN_t[v0idx];
            vector<float> v1 = vMPsXYZN_t[v1idx];
            vector<float> v2 = vMPsXYZN_t[v2idx];

            mi.push_back((v0[0]+v1[0]+v2[0])/3);
            mi.push_back((v0[1]+v1[1]+v2[1])/3);
            mi.push_back((v0[2]+v1[2]+v2[2])/3);
        }

        vMPsXYZN_t.push_back(mi);
    }

    //cout << "size top_layer(vMPsXYZN_t) = " << vMPsXYZN_t.size() << endl;

    vMPsXYZN_t2.clear();
    vMPsXYZN_t2.resize(vMPsXYZN_t.size());

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
	{
	    if (vMPsXYZN_t[i].empty())
            continue;

		vector<float> point;
		point.push_back(vMPsXYZN_t[i][0] - h);
		point.push_back(vMPsXYZN_t[i][1] - h);
		point.push_back(vMPsXYZN_t[i][2] - h);

		vMPsXYZN_t2[i] = point;
	}

	//cout << "size bottom_layer(vMPsXYZN_t2) = " << vMPsXYZN_t2.size() << endl;

    uf.clear();

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
    {
        if (vMPsXYZN_t[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t[i][j]);
    }

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
    {
        if (vMPsXYZN_t2[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t2[i][j]);
    }

    //cout << "u0 size = " << u0.size() << endl;
    //cout << "uf size = " << uf.size() << endl;
}


void FEA::ComputeDisplacement()
{
    a.clear();

    for (unsigned int i=0; i<u0.size(); i++)
        a.push_back(uf[i]-u0[i]);
}


void FEA::ComputeForces()   // K · a = f
{
    f.clear();
    f.resize(Ksize);
    for (unsigned int i=0; i<Ksize; i++)
    {
        float fi = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            if (K[i][j]!=0 && a[j]!=0)
                fi += K[i][j] * a[j];
        f[i] = fi;
    }
}


float FEA::ComputeStrainEnergy()
{
    // sE = a' · K · a
    vector<float> aK = vector<float>(Ksize,0.0);
    for (unsigned int i=0; i<Ksize; i++)
    {
        float prod = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            if (a[j]!=0)
                if (K[j][i]!=0)
                    prod += a[j] * K[j][i];
        aK[i] = prod;
    }

    sE = 0.0;
    for (unsigned int i=0; i<Ksize; i++)
        sE += aK[i] * a[i];

    CurrentSE = sE;

    return sE;
}

float FEA::NormalizeStrainEnergy()
{
    // Normalize: divide total sE by number of elements (Ksize/3).
    // nEls (tracked and extra) won't matter anymore.
    float invNormFactor = 1/fNormFactor;
    nsE = sE*invNormFactor;

    int nEl = Ksize/3;
    nsE = sE / nEl;

    return nsE;
}


float FEA::GetStrainEnergy() { return sE; }

float FEA::GetNormalizedStrainEnergy() { return nsE; }

void FEA::setbfea2(bool bSet) { bInFEA2 = bSet; }

void FEA::setbfea(bool bSet) { bInFEA = bSet; }

//void FEA::setCurrEdge(int input) { nCurrEdge = input; }
