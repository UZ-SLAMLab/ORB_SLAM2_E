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

// Invert matrix
//#include <cstdlib>
//#include <iomanip>


//FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, bool bSetDebug)
FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg(input_fg1), bDebugMode(bSetDebug)
{
    Ksize = 0.0;
    sE = 0.0;

    lambda = ( nu*E ) / ( (1+nu) * (1-2*nu) );
    G = E / ( 2 * (1+nu) );

    D.clear();
    D.push_back(vector<float>({ (lambda+2*G) ,    lambda    ,    lambda    ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({    lambda    , (lambda+2*G) ,    lambda    ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({    lambda    ,    lambda    , (lambda+2*G) ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,      G      ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,     0.0     ,      G      ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,     0.0     ,     0.0     ,      G      }));
 
    gs.clear();
    gs.push_back(vector<float>({-fg,-fg,-fg}));
    gs.push_back(vector<float>({+fg,-fg,-fg}));
    gs.push_back(vector<float>({+fg,+fg,-fg}));
    gs.push_back(vector<float>({-fg,+fg,-fg}));
    gs.push_back(vector<float>({-fg,-fg,+fg}));
    gs.push_back(vector<float>({+fg,-fg,+fg}));
    gs.push_back(vector<float>({+fg,+fg,+fg}));
    gs.push_back(vector<float>({-fg,+fg,+fg}));
}


FEA::~FEA()
{
}


vector<vector<float> > FEA::ComputeKei(vector<vector<float> > vfPts) {

    vector<vector<float> > vBtDB = vector<vector<float> >(24,vector<float>(24,0.0));

    for (unsigned int ops=0; ops<gs.size(); ops++)
    {
        float xi = gs[ops][0];
        float eta = gs[ops][1];
        float zeta = gs[ops][2];

        float dN1dxi = -0.125 * ( (1-eta) * (1-zeta) );     float dN1deta = -0.125 * ( (1-xi) * (1-zeta) );     float dN1dzeta = -0.125 * ( (1-xi) * (1-eta) ); 
        float dN2dxi = +0.125 * ( (1-eta) * (1-zeta) );     float dN2deta = -0.125 * ( (1+xi) * (1-zeta) );     float dN2dzeta = -0.125 * ( (1+xi) * (1-eta) ); 
        float dN3dxi = +0.125 * ( (1+eta) * (1-zeta) );     float dN3deta = +0.125 * ( (1+xi) * (1-zeta) );     float dN3dzeta = -0.125 * ( (1+xi) * (1+eta) ); 
        float dN4dxi = -0.125 * ( (1+eta) * (1-zeta) );     float dN4deta = +0.125 * ( (1-xi) * (1-zeta) );     float dN4dzeta = -0.125 * ( (1-xi) * (1+eta) ); 
        float dN5dxi = -0.125 * ( (1-eta) * (1+zeta) );     float dN5deta = -0.125 * ( (1-xi) * (1+zeta) );     float dN5dzeta = +0.125 * ( (1-xi) * (1-eta) ); 
        float dN6dxi = +0.125 * ( (1-eta) * (1+zeta) );     float dN6deta = -0.125 * ( (1+xi) * (1+zeta) );     float dN6dzeta = +0.125 * ( (1+xi) * (1-eta) ); 
        float dN7dxi = +0.125 * ( (1+eta) * (1+zeta) );     float dN7deta = +0.125 * ( (1+xi) * (1+zeta) );     float dN7dzeta = +0.125 * ( (1+xi) * (1+eta) ); 
        float dN8dxi = -0.125 * ( (1+eta) * (1+zeta) );     float dN8deta = +0.125 * ( (1-xi) * (1+zeta) );     float dN8dzeta = +0.125 * ( (1-xi) * (1+eta) ); 
                
        /*dxdxi*/   float J_00 = dN1dxi*vfPts[0][0] + dN2dxi*vfPts[1][0] + dN3dxi*vfPts[2][0] + dN4dxi*vfPts[3][0] + dN5dxi*vfPts[4][0] + dN6dxi*vfPts[5][0] + dN7dxi*vfPts[6][0] + dN8dxi*vfPts[7][0];
        /*dydxi*/   float J_01 = dN1dxi*vfPts[0][1] + dN2dxi*vfPts[1][1] + dN3dxi*vfPts[2][1] + dN4dxi*vfPts[3][1] + dN5dxi*vfPts[4][1] + dN6dxi*vfPts[5][1] + dN7dxi*vfPts[6][1] + dN8dxi*vfPts[7][1];
        /*dzdxi*/   float J_02 = dN1dxi*vfPts[0][2] + dN2dxi*vfPts[1][2] + dN3dxi*vfPts[2][2] + dN4dxi*vfPts[3][2] + dN5dxi*vfPts[4][2] + dN6dxi*vfPts[5][2] + dN7dxi*vfPts[6][2] + dN8dxi*vfPts[7][2];
        /*dxdeta*/  float J_10 = dN1deta*vfPts[0][0] + dN2deta*vfPts[1][0] + dN3deta*vfPts[2][0] + dN4deta*vfPts[3][0] + dN5deta*vfPts[4][0] + dN6deta*vfPts[5][0] + dN7deta*vfPts[6][0] + dN8deta*vfPts[7][0];
        /*dydeta*/  float J_11 = dN1deta*vfPts[0][1] + dN2deta*vfPts[1][1] + dN3deta*vfPts[2][1] + dN4deta*vfPts[3][1] + dN5deta*vfPts[4][1] + dN6deta*vfPts[5][1] + dN7deta*vfPts[6][1] + dN8deta*vfPts[7][1];
        /*dzdeta*/  float J_12 = dN1deta*vfPts[0][2] + dN2deta*vfPts[1][2] + dN3deta*vfPts[2][2] + dN4deta*vfPts[3][2] + dN5deta*vfPts[4][2] + dN6deta*vfPts[5][2] + dN7deta*vfPts[6][2] + dN8deta*vfPts[7][2];
        /*dxdzeta*/ float J_20 = dN1dzeta*vfPts[0][0] + dN2dzeta*vfPts[1][0] + dN3dzeta*vfPts[2][0] + dN4dzeta*vfPts[3][0] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0] + dN7dzeta*vfPts[6][0] + dN8dzeta*vfPts[7][0];
        /*dydzeta*/ float J_21 = dN1dzeta*vfPts[0][1] + dN2dzeta*vfPts[1][1] + dN3dzeta*vfPts[2][1] + dN4dzeta*vfPts[3][1] + dN5dzeta*vfPts[4][1] + dN6dzeta*vfPts[5][1] + dN7dzeta*vfPts[6][1] + dN8dzeta*vfPts[7][1];
        /*dzdzeta*/ float J_22 = dN1dzeta*vfPts[0][2] + dN2dzeta*vfPts[1][2] + dN3dzeta*vfPts[2][2] + dN4dzeta*vfPts[3][2] + dN5dzeta*vfPts[4][2] + dN6dzeta*vfPts[5][2] + dN7dzeta*vfPts[6][2] + dN8dzeta*vfPts[7][2];

        float Jac = J_00*J_11*J_22 + J_01*J_12*J_20 + J_10*J_21*J_02 - J_20*J_11*J_02 - J_10*J_01*J_22 - J_21*J_12*J_00;

        float J1_00 = (+1) * ( (J_11*J_22) - (J_21*J_12) ) / Jac;     float J1_01 = (-1) * ( (J_01*J_22) - (J_21*J_02) ) / Jac;     float J1_02 = (-1) * ( (J_01*J_12) - (J_11*J_02) ) / Jac;
        float J1_10 = (-1) * ( (J_10*J_22) - (J_20*J_12) ) / Jac;     float J1_11 = (-1) * ( (J_00*J_22) - (J_20*J_02) ) / Jac;     float J1_12 = (-1) * ( (J_00*J_12) - (J_10*J_02) ) / Jac;
        float J1_20 = (+1) * ( (J_10*J_21) - (J_20*J_11) ) / Jac;     float J1_21 = (-1) * ( (J_00*J_21) - (J_20*J_01) ) / Jac;     float J1_22 = (-1) * ( (J_00*J_11) - (J_10*J_01) ) / Jac;

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

        float BtD[24][6] = {0.0};

        for (unsigned int i=0; i<24; i++)
            for (unsigned int j=0; j<6; j++)
                BtD[i][j] = B[0][i]*D[0][j] + B[1][i]*D[1][j] + B[2][i]*D[2][j] + B[3][i]*D[3][j] + B[4][i]*D[4][j] + B[5][i]*D[5][j];

        for (unsigned int i=0; i<24; i++)
            for (unsigned int j=0; j<24; j++){
                float vBtDBaux = BtD[i][0]*B[0][j] + BtD[i][1]*B[1][j] + BtD[i][2]*B[2][j] + BtD[i][3]*B[3][j] + BtD[i][4]*B[4][j] + BtD[i][5]*B[5][j];
                vBtDB[i][j] += vBtDBaux * Jac;
            }
    }

    return vBtDB;
}


bool FEA::MatrixAssembly() {
    int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    Ksize = 3*nTotalNodes;
    cout << "                - MatAssembly (hex, 24 DoF per el. Curr: " << quads_t.size() << " quads, " << nTotalNodes << " nodes, Ksize = " << Ksize << ")" << endl;

    if (Ksize<=3)
        return false;

    K = vector<vector<float> >(Ksize,vector<float>(Ksize,0.0));

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

        vector<vector<float> > Kei = ComputeKei(vfPts);

        vector<int> mn;
        for (unsigned j=0; j<nodes.size(); j++) {
            mn.push_back(nodes[j]*3);
        }

        for (unsigned int ni=0; ni<nodes.size(); ni++){
            for (unsigned int nj=0; nj<nodes.size(); nj++){
                for (unsigned int m=0; m<3; m++){
                    for (unsigned int n=0; n<3; n++){
                        K[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                    }
                }
            }
        }
    }
    return true;
}



vector<vector<float> > FEA::InvertMatrixEigen(vector<vector<float> > m1) { 
	Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A;
	A.resize(m1.size(),m1.size());

	for (unsigned int i=0; i<m1.size(); i++){
		for (unsigned int j=0; j<m1[i].size(); j++){
			A(i,j) = m1[i][j];
		}
	}

    float detA = A.determinant();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A1 = A.inverse();

    vector<vector<float> > Ainv;
    if (detA!=0.0){
        for (unsigned int i=0; i<m1.size(); i++){
            vector<float> Ainvi;
            for (unsigned int j=0; j<m1[i].size(); j++){
                Ainvi.push_back(A1(i,j));
            }
            Ainv.push_back(Ainvi);
        }
    }
    else{
        vector<float> Ainvi;
        Ainvi.push_back(0.0);
        Ainv.push_back(Ainvi);
    }

    return Ainv;
} 


vector<vector<float> > FEA::MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2) {
    vector<vector<float> > pr = vector<vector<float> >(m1.size(),vector<float>(m2[0].size(),0.0));

    if (m1[0].size() != m2.size()) {
        pr = vector<vector<float> >(1,vector<float>(1,0.0));
        return pr;
    }

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A;
	A.resize(m1.size(),m1[0].size());
	for (unsigned int i=0; i<m1.size(); i++){
		for (unsigned int j=0; j<m1[i].size(); j++){
			A(i,j) = m1[i][j];
		}
	}

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> B;
	B.resize(m2.size(),m2[0].size());
	for (unsigned int i=0; i<m2.size(); i++){
		for (unsigned int j=0; j<m2[i].size(); j++){
			B(i,j) = m2[i][j];
		}
	}

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> C;
	C.resize(pr.size(),pr[0].size());
    C = A*B;

    for (unsigned int i=0; i<pr.size(); i++){
        for (unsigned int j=0; j<pr[i].size(); j++){
            pr[i][j] = C(i,j);
        }
    }

    return pr;
}


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


void FEA::Set_uf(vector<vector<float> > vPoints){
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

    vMPsXYZN_t2.clear();
    vMPsXYZN_t2.resize(vMPsXYZN_t.size());

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++){
	    if (vMPsXYZN_t[i].empty())
            continue;

		vector<float> point;
		point.push_back(vMPsXYZN_t[i][0] - h);
		point.push_back(vMPsXYZN_t[i][1] - h);
		point.push_back(vMPsXYZN_t[i][2] - h);

		vMPsXYZN_t2[i] = point;
	}

    uf.clear();

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++){
        if (vMPsXYZN_t[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t[i][j]);
    }

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++){
        if (vMPsXYZN_t2[i].empty())
            for (unsigned int j=0; j<3; j++)
                uf.push_back(0.0);
        else
            for (unsigned int j=0; j<3; j++)
                uf.push_back(vMPsXYZN_t2[i][j]);
    }
}


void FEA::ComputeDisplacement(){
    a.clear();
    for (unsigned int i=0; i<u0.size(); i++)
        a.push_back(uf[i]-u0[i]);
}


void FEA::ComputeForces(){
    f.clear();
    f.resize(Ksize);
    for (unsigned int i=0; i<Ksize; i++){
        float fi = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            if (K[i][j]!=0 && a[j]!=0)
                fi += K[i][j] * a[j];
        f[i] = fi;
    }
}


float FEA::ComputeStrainEnergy(){
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

float FEA::NormalizeStrainEnergy(){
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
