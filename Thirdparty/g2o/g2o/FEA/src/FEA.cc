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


//FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, bool bSetDebug)
FEA::FEA(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg(input_fg1), bDebugMode(bSetDebug){
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


FEA::~FEA(){
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


void FEA::Set_u0(){
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

    for (unsigned int i=0; i<vPoints.size(); i++){
        vector<float> vfMP;

        vfMP.push_back(vPoints[i][0]);
        vfMP.push_back(vPoints[i][1]);
        vfMP.push_back(vPoints[i][2]);

        vMPsXYZN_t[i] = vfMP;
    }

    for (unsigned int i=0; i<vNewPointsBase.size(); i++){
        vector<float> mi;

        if (vNewPointsBase[i].size()==2){
            int v0idx = vNewPointsBase[i][0];
            int v1idx = vNewPointsBase[i][1];

            vector<float> v0 = vMPsXYZN_t[v0idx];
            vector<float> v1 = vMPsXYZN_t[v1idx];

            mi.push_back((v0[0]+v1[0])/2);
            mi.push_back((v0[1]+v1[1])/2);
            mi.push_back((v0[2]+v1[2])/2);
        }
        if(vNewPointsBase[i].size()==3){
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
    vva.clear();
    for (unsigned int i=0; i<u0.size(); i++){
        vector<float> va;
        va.push_back(uf[i]-u0[i]);
        vva.push_back(va);
    }
}


void FEA::ComputeForces(){
    vvf.clear();
    vvf = MultiplyMatricesEigen(K,vva);
}


float FEA::ComputeStrainEnergy(){
    // sE = a' · K · a = a' · F

    vector<vector<float> > vvat;
    vector<float> vvati;
    for (unsigned int i=0; i<vva.size(); i++){
        vvati.push_back(vva[i][0]);
    }
    vvat.push_back(vvati);

    vector<vector<float> > vvsE = MultiplyMatricesEigen(vvat,vvf);
    sE = vvsE[0][0];

    CurrentSE = sE;
    return sE;
}

float FEA::NormalizeStrainEnergy(){
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
