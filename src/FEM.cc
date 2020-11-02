/**
* This file is an add-on to ORB-SLAM2.
*
* Its goal is to create a mesh from the 3D points within the computed map.
* Also, it will assembly a rigidity matrix for such mesh, and compute the
* strain deformation energy.
*
* Copyright (C) 2017 Íñigo Cirauqui Viloria <icirauqui at gmail dot com> (University of Zaragoza)
*
* This add-on, as ORB-SLAM2, is free software: you can redistribute it and/or
* modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your
* option) any later version.
*
* This add-on is distributed in the hope that it will be useful,
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

#include "FEM.h"

#include <iostream>
#include <fstream>

#include <opencv2/highgui/highgui.hpp>

#include <thread>

// Invert matrix
#include <cstdlib>
#include <iomanip>


namespace ORB_SLAM2
{


/// 0 - CONSTRUCTOR
FEM::FEM(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg1(input_fg1), fg2(-input_fg1), fg(input_fg1), bDebugMode(bSetDebug)
{
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

    gs.push_back(vector<float>({+fg,+fg,+fg}));
    gs.push_back(vector<float>({+fg,+fg,-fg}));
    gs.push_back(vector<float>({+fg,-fg,+fg}));
    gs.push_back(vector<float>({+fg,-fg,-fg}));
    gs.push_back(vector<float>({-fg,+fg,+fg}));
    gs.push_back(vector<float>({-fg,+fg,-fg}));
    gs.push_back(vector<float>({-fg,-fg,+fg}));
    gs.push_back(vector<float>({-fg,-fg,-fg}));

}


FEM::~FEM()
{
}


/// 1 - MESHING

bool FEM::Compute(int nMode) {

    GetMapPointCoordinates(nMode);

    if (LoadMPsIntoCloud(nMode)) {
        cout << "           - MLS" << endl;
        if (MLS(nMode)) {
            cout << "           - Tri mesh" << endl;
            if(ComputeMesh(nMode)) {
                tri2quad(nMode);
                SetSecondLayer(nMode);

                if (nMode == 1)
                    Set_u0(vpMPs_t,nMode);
                else if (nMode == 2)
                    Set_u0(vpMPs_ut,nMode);

                MatrixAssembly(nMode);

                if (nMode == 2)
                    ComputeK1();

                return true;
            }
        }
    }
    return false;


    /*

    GetMapPointCoordinates(1);
    GetMapPointCoordinates(2);

    if (LoadMPsIntoCloud(1) && LoadMPsIntoCloud(2))
    {
        cout << "           - MLS" << endl;
        if (MLS(1) && MLS(2))
        {
            cout << "           - Tri mesh" << endl;
            if(ComputeMesh(1) && ComputeMesh(2))
            {
                tri2quad_t();
                tri2quad_u(); //to be reviewed
                //tri2quad_ut();
                SetSecondLayer(1);
                SetSecondLayer(2);
                Set_u0(vpMPs_t,1);
                Set_u0(vpMPs_ut,2);

                MatrixAssembly(1);
                MatrixAssembly(2);
                ComputeK1();

                return true;
            }
        }
    }
    return false;

    */
}


void FEM::ReadyFEA(FEA* pFEA, int nMode) {
    if (nMode == 1) {
        pFEA->triangles_t = triangles_t;
        pFEA->quads_t = quads_t;
        pFEA->vIdxNewVertices = vIdxNewVertices;
        pFEA->vNewPointsBase = vNewPointsBase;
    }


    /*
    pFEA->Ksize = Ksize;
    pFEA->u0 = u0;
    pFEA->vMPsXYZN_t  = vMPsXYZN_t;
    pFEA->vMPsXYZN_t2 = vMPsXYZN_t2;

    bool bMatAssembly = pFEA->MatrixAssembly();
    if (bMatAssembly && nMode==1)
        bool bMatInverse = pFEA->ComputeK1();
    */
}


/// New functions

void FEM::GetMapPointCoordinates(int nMode) {
    if (nMode == 1) {
        for (unsigned int i=0; i<vpMPs_t.size(); i++) {
            MapPoint* pMP = vpMPs_t[i];

            vector<float> vfMP;

            cv::Mat pos = pMP->GetWorldPos();
            cv::Vec3f vPos = pos;
            vfMP.push_back(vPos.val[0]);
            vfMP.push_back(vPos.val[1]);
            vfMP.push_back(vPos.val[2]);

            vMPsXYZN_t.push_back(vfMP);
            vbtMPsActive.push_back(true);
            vMPsGlobalID_t.push_back(pMP->mnId);
        }
    }
    else {
        for (unsigned int i=0; i<vpMPs_ut.size(); i++)
        {
            MapPoint* pMP = vpMPs_ut[i];

            vector<float> vfMP;

            cv::Mat pos = pMP->GetWorldPos();
            cv::Vec3f vPos = pos;
            vfMP.push_back(vPos.val[0]);
            vfMP.push_back(vPos.val[1]);
            vfMP.push_back(vPos.val[2]);

            vMPsXYZN_ut.push_back(vfMP);
            vbuMPsActive.push_back(true);
            vMPsGlobalID_u.push_back(pMP->mnId);
        }
    }
}


bool FEM::LoadMPsIntoCloud(int nMode) {
    if (nMode == 1) {
        // Load top layer into cloud
        pc_t_0.width = vMPsXYZN_t.size(); //cloudTracked.width    = 5;
        pc_t_0.height   = 1;
        pc_t_0.is_dense = false;
        pc_t_0.points.resize (pc_t_0.width * pc_t_0.height);

        for (size_t i = 0; i < pc_t_0.points.size(); i++) {
            pc_t_0.points[i].x = vMPsXYZN_t[i][0];
            pc_t_0.points[i].y = vMPsXYZN_t[i][1];
            pc_t_0.points[i].z = vMPsXYZN_t[i][2];
        }

        if (pc_t_0.width>0) {
            pcPath_t = pcPath_t + std::to_string(nFrameId) + ".pcd";
            pcl::io::savePCDFileASCII (pcPath_t, pc_t_0);
            return true;
        }
        else
            return false;
    }

    else if (nMode == 2) {
        pc_u_0.width = vMPsXYZN_ut.size(); //cloudTracked.width    = 5;
        pc_u_0.height   = 1;
        pc_u_0.is_dense = false;
        pc_u_0.points.resize (pc_u_0.width * pc_u_0.height);

        for (size_t i = 0; i < vMPsXYZN_ut.size(); i++) {
            pc_u_0.points[i].x = vMPsXYZN_ut[i][0];
            pc_u_0.points[i].y = vMPsXYZN_ut[i][1];
            pc_u_0.points[i].z = vMPsXYZN_ut[i][2];
        }

        if (pc_u_0.width>0) {
            pcPath_ut = pcPath_ut + std::to_string(nFrameId) + ".pcd";
            pcl::io::savePCDFileASCII (pcPath_ut, pc_u_0);
            return true;
        }
        else
            return false;
    }
}


bool FEM::MLS(int nMode) {
    if (nMode == 1 ) {
        PointCloud<PointXYZ>::Ptr pc_t_mls_0 (new pcl::PointCloud<PointXYZ> (pc_t_0));

        // Create KD-Tree & MLS & set parameters
        pcl::search::KdTree<pcl::PointXYZ>::Ptr tree1 (new pcl::search::KdTree<pcl::PointXYZ>);
        pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls1;
        mls1.setComputeNormals (true);
        mls1.setInputCloud (pc_t_mls_0);
        mls1.setPolynomialFit (true);
        mls1.setSearchMethod (tree1);
        mls1.setSearchRadius (900000);	//100 //0.25    // Original 0.03 // Default 0.00
        mls1.setPolynomialOrder(3); //50 // Default 2

        mls1.process (pc_t_1);

        // Get corresponding indexes: for each output point, returns the index of the input one.
        PointIndicesPtr pIdx1 = mls1.getCorrespondingIndices();
        for (unsigned int i=0; i<pIdx1->indices.size(); i++)
            vtindices.push_back(pIdx1->indices[i]);

        int idxit = 0;
        for (unsigned int i=0; i<vMPsXYZN_t.size(); i++) {
            int currentpos = i;
            if (currentpos==vtindices[idxit]) {
                vbtMPsActive[i] = true;
                idxit++;
            }
            else
                vbtMPsActive[i] = false;
        }

        /// Save smoothed cloud
        if (pc_t_1.width>0) {
            pcsPath_t = pcsPath_t + std::to_string(nFrameId) + ".pcd";
            pcl::io::savePCDFileASCII (pcsPath_t, pc_t_1);
            return true;
        }
        else
            return false;
    }

    else if (nMode == 2) {
        PointCloud<PointXYZ>::Ptr pc_u_mls_0 (new pcl::PointCloud<PointXYZ> (pc_u_0));

        // Create KD-Tree & MLS & set parameters
        pcl::search::KdTree<pcl::PointXYZ>::Ptr tree1 (new pcl::search::KdTree<pcl::PointXYZ>);
        pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls1;
        mls1.setComputeNormals (true);
        mls1.setInputCloud (pc_u_mls_0);
        mls1.setPolynomialFit (true);
        mls1.setSearchMethod (tree1);
        mls1.setSearchRadius (900000);	//100 //0.25    // Original 0.03 // Default 0.00
        mls1.setPolynomialOrder(3); //50 // Default 2

        mls1.process (pc_u_1);

        // Get corresponding indexes: for each output point, returns the index of the input one.
        PointIndicesPtr pIdx1 = mls1.getCorrespondingIndices();
        for (unsigned int i=0; i<pIdx1->indices.size(); i++)
            vuindices.push_back(pIdx1->indices[i]);

        int idxit = 0;
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++) {
            int currentpos = i;
            if (currentpos==vuindices[idxit]) {
                vbuMPsActive[i] = true;
                idxit++;
            }
            else
                vbuMPsActive[i] = false;
        }

        /// Save smoothed cloud
        if (pc_u_1.width>0) {
            pcsPath_ut = pcsPath_ut + std::to_string(nFrameId) + ".pcd";
            pcl::io::savePCDFileASCII (pcsPath_ut, pc_u_1);
            return true;
        }
        else
            return false;
    }

}


bool FEM::ComputeMesh(int nMode) {
    // TRIANGULATION, GREEDY PROJECTION
    if (nMode == 1) {
        PointCloud<PointNormal>::Ptr ptr_pc_t_mesh_1 (new pcl::PointCloud<PointNormal> (pc_t_1));

        // Search tree
        pcl::search::KdTree<pcl::PointNormal>::Ptr tree (new pcl::search::KdTree<pcl::PointNormal>);
        tree->setInputCloud (ptr_pc_t_mesh_1);

        // Greedy Projection parameters
        int nInMu = 5;
        int nSearchRad = 30;
        int nMaxNeig = 70;
        int nSurAng = 150;
        int nMinAng = 45;
        int nMaxAng = 90;
        CalculateGP3Parameters(ptr_pc_t_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);

        // Greedy Projection
        static pcl::GreedyProjectionTriangulation<pcl::PointNormal> GreedyProj3;
        GreedyProj3.setMu (nInMu);    // original 2.8
        GreedyProj3.setSearchRadius(nSearchRad);      // original 0.50
        GreedyProj3.setMaximumNearestNeighbors (nMaxNeig);        // original 150
        GreedyProj3.setMaximumSurfaceAngle(nSurAng); // 45 degrees
        GreedyProj3.setMinimumAngle(nMinAng); // 10 degrees
        GreedyProj3.setMaximumAngle(nMaxAng); // 120 degrees 2/3
        GreedyProj3.setNormalConsistency(true);

        /*cout << "GP3 Params: " << nInMu << " " << nSearchRad << " " << nMaxNeig << endl;*/

        bool bGo = true;
        if (bGo == true)
        {
            GreedyProj3.setMu (5.5);    // original 2.8
            GreedyProj3.setSearchRadius(0.75);      // original 0.50
            GreedyProj3.setMaximumNearestNeighbors (50);        // original 150
            GreedyProj3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
            GreedyProj3.setMinimumAngle(M_PI/18); // 10 degrees
            GreedyProj3.setMaximumAngle(9*M_PI/10); // 120 degrees 2/3
            GreedyProj3.setNormalConsistency(true);
        }

        GreedyProj3.setInputCloud (ptr_pc_t_mesh_1);
        GreedyProj3.setSearchMethod (tree);
        GreedyProj3.reconstruct (mesh_t);

        if (mesh_t.polygons.size()>0) {
            // Get the vertex indexes of each triangle
            for (unsigned int i=0; i<mesh_t.polygons.size(); i++) {
                unsigned int nver0 = mesh_t.polygons[i].vertices[0];
                unsigned int nver1 = mesh_t.polygons[i].vertices[1];
                unsigned int nver2 = mesh_t.polygons[i].vertices[2];

                vector<int> triangle;
                triangle.push_back(nver0);
                triangle.push_back(nver1);
                triangle.push_back(nver2);
                triangles_t.push_back(triangle);

                vector<MapPoint*> vpMP2Draw;
                vpMP2Draw.push_back(vpMPs_t[vtindices[nver0]]);
                vpMP2Draw.push_back(vpMPs_t[vtindices[nver1]]);
                vpMP2Draw.push_back(vpMPs_t[vtindices[nver2]]);
                vpMPs2Draw.push_back(vpMP2Draw);

                vector<cv::KeyPoint*> vpKP2Draw;
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver0]]);
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver1]]);
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver2]]);
                vpKPs2Draw.push_back(vpKP2Draw);
            }

            // Save reconstruction
            pcrPath_t = pcrPath_t + std::to_string(nFrameId) + ".vtk";
            pcl::io::saveVTKFile (pcrPath_t, mesh_t);
            return true;
        }
        else {
            cout << "No data" << endl;
            return false;
        }
    }

    else if (nMode == 2) {
        PointCloud<PointNormal>::Ptr ptr_pc_u_mesh_1 (new pcl::PointCloud<PointNormal> (pc_u_1));

        // Search tree
        pcl::search::KdTree<pcl::PointNormal>::Ptr tree (new pcl::search::KdTree<pcl::PointNormal>);
        tree->setInputCloud (ptr_pc_u_mesh_1);

        // Greedy Projection parameters
        int nInMu = 5;
        int nSearchRad = 30;
        int nMaxNeig = 70;
        int nSurAng = 150;
        int nMinAng = 45;
        int nMaxAng = 90;
        CalculateGP3Parameters(ptr_pc_u_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);

        // Greedy Projection
        static pcl::GreedyProjectionTriangulation<pcl::PointNormal> GreedyProj3;
        GreedyProj3.setMu (nInMu);    // original 2.8
        GreedyProj3.setSearchRadius(nSearchRad);      // original 0.50
        GreedyProj3.setMaximumNearestNeighbors (nMaxNeig);        // original 150
        GreedyProj3.setMaximumSurfaceAngle(nSurAng); // 45 degrees
        GreedyProj3.setMinimumAngle(nMinAng); // 10 degrees
        GreedyProj3.setMaximumAngle(nMaxAng); // 120 degrees 2/3
        GreedyProj3.setNormalConsistency(true);

        bool bGo = true;
        if (bGo == true)
        {
            GreedyProj3.setMu (5.5);    // original 2.8
            GreedyProj3.setSearchRadius(0.75);      // original 0.50
            GreedyProj3.setMaximumNearestNeighbors (50);        // original 150
            GreedyProj3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
            GreedyProj3.setMinimumAngle(M_PI/18); // 10 degrees
            GreedyProj3.setMaximumAngle(9*M_PI/10); // 120 degrees 2/3
            GreedyProj3.setNormalConsistency(true);
        }

        GreedyProj3.setInputCloud (ptr_pc_u_mesh_1);
        GreedyProj3.setSearchMethod (tree);
        GreedyProj3.reconstruct (mesh_u);

        if (mesh_u.polygons.size()>0) {
            for (unsigned int i=0; i<mesh_u.polygons.size(); i++) {
                unsigned int nver0 = mesh_u.polygons[i].vertices[0];
                unsigned int nver1 = mesh_u.polygons[i].vertices[1];
                unsigned int nver2 = mesh_u.polygons[i].vertices[2];

                vector<int> triangle;
                triangle.push_back(nver0);
                triangle.push_back(nver1);
                triangle.push_back(nver2);
                triangles_u.push_back(triangle);

                vector<MapPoint*> vpMP2Draw;
                vpMP2Draw.push_back(vpMPs_ut[vuindices[nver0]]);
                vpMP2Draw.push_back(vpMPs_ut[vuindices[nver1]]);
                vpMP2Draw.push_back(vpMPs_ut[vuindices[nver2]]);
                vpMPs2Drawu.push_back(vpMP2Draw);

                /*
                vector<cv::KeyPoint*> vpKP2Draw;
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver0]]);
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver1]]);
                vpKP2Draw.push_back(vpKPs_t[vtindices[nver2]]);
                vpKPs2Drawu.push_back(vpKP2Draw);
                */
            }

            // Save reconstruction
            pcrPath_u = pcrPath_u + std::to_string(nFrameId) + ".vtk";
            pcl::io::saveVTKFile (pcrPath_u, mesh_u);
            return true;
        }
        else {
            cout << "No data" << endl;
            return false;
        }
    }

}


void FEM::CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc, int *pnInMu, int *pnSearchRad, int *pnMaxNeig, int *pnSurAng, int *pnMinAng, int *pnMaxAng) {

    vector<vector<float> > vit1;
    for (unsigned int i=0; i<ppc->points.size(); i++) {
        vector<float> viti;
        viti.push_back(ppc->points[i].x);
        viti.push_back(ppc->points[i].y);
        viti.push_back(ppc->points[i].z);
        vit1.push_back(viti);
    }

    vector<float> vtot, vdmin, vdmax;
    float dd, dx, dy, dz;
    for (unsigned int i=0; i<vit1.size(); i++) {
        float dmin = 1000.0;
        float dmax = 0.0;
        for (unsigned int j=0; j<vit1.size(); j++) {
            if (j==i)
                continue;

            dx = vit1[j][0] - vit1[i][0];
            dy = vit1[j][1] - vit1[i][1];
            dz = vit1[j][2] - vit1[i][2];
            dx *= dx;
            dy *= dy;
            dz *= dz;
            dd = dx + dy + dz;
            dd = sqrt(dd);

            if (dd < dmin)
                dmin = dd;
            if (dd > dmax)
                dmax = dd;

            vtot.push_back(dd);
        }
        vdmin.push_back(dmin);
        vdmax.push_back(dmax);
    }


    float fmedtot = 0.0;
    float fsigmatot = 0.0;
    int ntot = 0;
    for (unsigned int i=0; i<vtot.size(); i++) {
        fmedtot += vtot[i]/vtot.size();
    }
    for (unsigned int i=0; i<vtot.size(); i++) {
        fsigmatot += (vtot[i]-fmedtot)*(vtot[i]-fmedtot);
    }
    fsigmatot /= vtot.size();
    fsigmatot = sqrt(fsigmatot);


    float ftminmin = 1000.0;
    float ftmaxmin = 0.0;
    float ftminmax = 1000.0;
    float ftmaxmax = 0.0;
    for (unsigned int i=0; i<vdmin.size(); i++) {
        if (vdmin[i]<ftminmin)
            ftminmin = vdmin[i];
        if (vdmin[i]>ftmaxmin)
            ftmaxmin = vdmin[i];
        if (vdmax[i]<ftminmax)
            ftminmax = vdmax[i];
        if (vdmax[i]>ftmaxmax)
            ftmaxmax = vdmax[i];
    }

    float fmedmin = (ftmaxmin-ftminmin)/2;
    float fmedmax = (ftmaxmax-ftminmax)/2;
    float fsigmamin = 0.0;
    float fsigmamax = 0.0;
    for (unsigned int i=0; i<vdmin.size(); i++) {
        fsigmamin += (vdmin[i]-fmedmin)*(vdmin[i]-fmedmin);
        fsigmamax += (vdmax[i]-fmedmax)*(vdmax[i]-fmedmax);
        fsigmamin /= vdmin.size();
        fsigmamax /= vdmax.size();
        fsigmamin = sqrt(fsigmamin);
        fsigmamax = sqrt(fsigmamax);
    }

    /*
    cout << " -- Min = " << fmedmin << " - " << fsigmamin << endl;
    cout << " -- Max = " << fmedmax << " - " << fsigmamax << endl;
    cout << " -- Tot = " << fmedtot << " - " << fsigmatot << endl;
    cout << " Search Radious " << ceil(fmedtot) << endl;
    */

    *pnInMu = ceil(fmedtot/fmedmin);
    *pnSearchRad = ceil(fmedtot);
    *pnMaxNeig = vtot.size() / 4;//fmedtot + fsigmatot;
    *pnSurAng = 150;
    *pnMinAng = 45;
    *pnMaxAng = 180-2*45;
}


void FEM::tri2quad(int nMode) {
        if (nMode == 1)
            tri2quad_t();
        else if (nMode == 2)
            tri2quad_u();
}


void FEM::tri2quad_t() {
  	// Get all vertices in mesh
  	vector<int> vertices;
  	for (unsigned int i=0; i<triangles_t.size(); i++)
  		for (unsigned int j=0; j<triangles_t[i].size(); j++)
  			if ( !( std::find(vertices.begin(), vertices.end(), triangles_t[i][j]) != vertices.end() ) )
            	vertices.push_back(triangles_t[i][j]);

	int maxvertex = vertices.back();

  	vector<int> offvertices;
  	for (int i=0; i<maxvertex; i++)
  		if ( !( std::find(vertices.begin(), vertices.end(), i) != vertices.end() ) )
            offvertices.push_back(i);

    for (unsigned int i=0; i<vertices.size(); i++)
    	vbtMPsActive[vtindices[vertices[i]]] = true;

    for (unsigned int i=0; i<offvertices.size(); i++)
    	vbtMPsActive[vtindices[offvertices[i]]] = false;

  	// Get true indexes
  	vector<vector<int> > tempidx;
  	for (unsigned int i=0; i<triangles_t.size(); i++) {
  		vector<int> triangle;
  		for (unsigned int j=0; j<triangles_t[i].size(); j++)
  			triangle.push_back(vtindices[triangles_t[i][j]]);
  		tempidx.push_back(triangle);
  	}

  	triangles_t.clear();
  	for (unsigned int i=0; i<tempidx.size(); i++)
  		triangles_t.push_back(tempidx[i]);


  	// Record the 3 edges_t in each triangle
    for (unsigned int i=0; i<triangles_t.size(); i++) {
        if (triangles_t[i].size() < 3)
            continue;

        vector<int> edge1;
        vector<int> edge2;
        vector<int> edge3;
        edge1.push_back(triangles_t[i][0]);
        edge1.push_back(triangles_t[i][1]);
        edge2.push_back(triangles_t[i][0]);
        edge2.push_back(triangles_t[i][2]);
        edge3.push_back(triangles_t[i][1]);
        edge3.push_back(triangles_t[i][2]);

        vector<vector<int> > tempedges;
        tempedges.push_back(edge1);
        tempedges.push_back(edge2);
        tempedges.push_back(edge3);

        edgespertriangle_t.push_back(tempedges);
    }



	vector<vector<vector<int> > > vCommonEdges;
	vCommonEdges.resize(triangles_t.size());
	for (unsigned int i=0; i<vCommonEdges.size(); i++)
		vCommonEdges[i].resize(3);
		/*triangulo1
			commonedge1
				edge1 triangulo2 edge2
			commonedge2
				edge1 triangulo2 edge2
			commonedge3
				edge1 triangulo2 edge2
		triangulo2
		triangulo3
		...*/

	vIdxNewVertices.clear();
	vIdxNewVertices.resize(triangles_t.size());
	//for(unsigned int i=0; i<vIdxNewVertices.size(); i++)
		//for (unsigned int j=0; j<4; j++)
			//vIdxNewVertices[i].push_back(-1);
		/*triangulo1
			idxvertexedge0 idxvertexedge1 idxvertexedge2 idxbaricentro
		triangulo2
		triangulo3
		...*/

	vector<bool> vbComputedTriangle = vector<bool>(triangles_t.size(), false);
		/*triangulo1
		triangulo2
		triangulo3
		...*/

	quads_t.clear();

	// Fill vCommonEdges with the function originally designed to join common edges.
    for (unsigned int i=0; i<edgespertriangle_t.size(); i++) {
        // Won't happen, but just in case
        if (edgespertriangle_t[i].size()<3)
            continue;

        // For each triangle, we look for common edges with all the other polygons
        for (unsigned int j=0; j<edgespertriangle_t.size(); j++) {
            // Don't look within itself
            if (j==i)
                continue;

            // Again, won't happen, but just in case
            if (edgespertriangle_t[j].size()<3)
                continue;

            // For all the edges within the reference triangle
            for (unsigned int k=0; k<edgespertriangle_t[i].size(); k++) {
                // Compare the first edge in the reference triangle with all the three(k)
                // edges in the candidate
                if  (
                        ( (edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 0 & 1
						vCommonEdges[i][0].push_back(0);
						vCommonEdges[i][0].push_back(j);
						vCommonEdges[i][0].push_back(k);
                    }

                // Compare the second edge in the reference triangle with all the three(k)
                // edges in the candidate
                if  (
                        ( (edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 0 & 2
						vCommonEdges[i][1].push_back(1);
						vCommonEdges[i][1].push_back(j);
						vCommonEdges[i][1].push_back(k);
                    }

                // Compare the third edge in the reference triangle with all the three(k)
                // edges in the candidate
                if  (
                        ( (edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 1 & 2
						vCommonEdges[i][2].push_back(2);
						vCommonEdges[i][2].push_back(j);
						vCommonEdges[i][2].push_back(k);
                    }
                    // k changes here, so in each loop, we compare all the edges in
                    // the reference triangle with only one edge of the candidate
            }
        }   // End of the loop through all the candidates for a reference triangle
    }   // End of the loop through all the triangles set as reference

    //cout << "edgespertrianglecomputed" << endl;

    vNewPointsBase.clear();

	for (unsigned int i=0; i<triangles_t.size(); i++) {
		// Get vertices and indices for refferece triangle
		int v0idx = triangles_t[i][0];
		vector<float> v0 = vMPsXYZN_t[v0idx];

		int v1idx = triangles_t[i][1];
		vector<float> v1 = vMPsXYZN_t[v1idx];

		int v2idx = triangles_t[i][2];
		vector<float> v2 = vMPsXYZN_t[v2idx];

		// Set variables to store the ressult
		vector<float> m01;
		int m01idx;
		vector<float> m02;
		int m02idx;
		vector<float> m12;
		int m12idx;
		vector<float> g;
		int gidx;
		vector<int> quad1;
		vector<int> quad2;
		vector<int> quad3;

		// Compute edge 0
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 0) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m01idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m01 = vMPsXYZN_t[m01idx];
					}
					else {
						m01idx = vMPsXYZN_t.size();

						m01.push_back((v0[0]+v1[0])/2);
		 				m01.push_back((v0[1]+v1[1])/2);
		  				m01.push_back((v0[2]+v1[2])/2);

						vMPsXYZN_t.push_back(m01);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v1idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m01idx = vMPsXYZN_t.size();

				m01.push_back((v0[0]+v1[0])/2);
		 		m01.push_back((v0[1]+v1[1])/2);
		  		m01.push_back((v0[2]+v1[2])/2);

				vMPsXYZN_t.push_back(m01);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v1idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute edge 1
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 1) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m02idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m02 = vMPsXYZN_t[m02idx];
					}
					else {
						m02idx = vMPsXYZN_t.size();

					  	m02.push_back((v0[0]+v2[0])/2);
					  	m02.push_back((v0[1]+v2[1])/2);
					  	m02.push_back((v0[2]+v2[2])/2);

						vMPsXYZN_t.push_back(m02);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m02idx = vMPsXYZN_t.size();

			 	m02.push_back((v0[0]+v2[0])/2);
				m02.push_back((v0[1]+v2[1])/2);
				m02.push_back((v0[2]+v2[2])/2);

				vMPsXYZN_t.push_back(m02);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute edge 2
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 2) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m12idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m12 = vMPsXYZN_t[m12idx];
					}
					else {
						m12idx = vMPsXYZN_t.size();

					  	m12.push_back((v1[0]+v2[0])/2);
					  	m12.push_back((v1[1]+v2[1])/2);
					  	m12.push_back((v1[2]+v2[2])/2);

						vMPsXYZN_t.push_back(m12);

						NewPointsBasei.push_back(v1idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m12idx = vMPsXYZN_t.size();

				m12.push_back((v1[0]+v2[0])/2);
				m12.push_back((v1[1]+v2[1])/2);
				m12.push_back((v1[2]+v2[2])/2);

				vMPsXYZN_t.push_back(m12);

				NewPointsBasei.push_back(v1idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute barycentre in each triangle
		g.push_back((v0[0]+v1[0]+v2[0])/3);
		g.push_back((v0[1]+v1[1]+v2[1])/3);
		g.push_back((v0[2]+v1[2]+v2[2])/3);
		gidx = vMPsXYZN_t.size();
		vMPsXYZN_t.push_back(g);

		vector<int> NewPointsBasei;
        NewPointsBasei.push_back(v0idx);
        NewPointsBasei.push_back(v1idx);
        NewPointsBasei.push_back(v2idx);
        vNewPointsBase.push_back(NewPointsBasei);

		// Set the three quads
		quad1.push_back(v0idx);
		quad1.push_back(m01idx);
		quad1.push_back(gidx);
		quad1.push_back(m02idx);

		quad2.push_back(v1idx);
		quad2.push_back(m12idx);
		quad2.push_back(gidx);
		quad2.push_back(m01idx);

		quad3.push_back(v2idx);
		quad3.push_back(m02idx);
		quad3.push_back(gidx);
		quad3.push_back(m12idx);

		// Add the quads to the global vector
		quads_t.push_back(quad1);
		quads_t.push_back(quad2);
		quads_t.push_back(quad3);

		// Set vertices idx
		vIdxNewVertices[i].push_back(m01idx);
		vIdxNewVertices[i].push_back(m02idx);
		vIdxNewVertices[i].push_back(m12idx);
		vIdxNewVertices[i].push_back(gidx);

		// Set triangle as already computed
		vbComputedTriangle[i] = true;
	}
}


void FEM::tri2quad_u() {

  	vector<int> vertices, offvertices;

  	for (unsigned int i=0; i<triangles_u.size(); i++)
  		for (unsigned int j=0; j<triangles_u[i].size(); j++)
  			if ( !( std::find(vertices.begin(), vertices.end(), triangles_u[i][j]) != vertices.end() ) )
            	vertices.push_back(triangles_u[i][j]);

  	for (int i=0; i<vertices.back(); i++)
  		if ( !( std::find(vertices.begin(), vertices.end(), i) != vertices.end() ) )
            offvertices.push_back(i);

    for (unsigned int i=0; i<vertices.size(); i++)
    	vbuMPsActive[vuindices[vertices[i]]] = true;

    for (unsigned int i=0; i<offvertices.size(); i++)
    	vbuMPsActive[vuindices[offvertices[i]]] = false;

  	vector<vector<int> > tempidx;
  	for (unsigned int i=0; i<triangles_u.size(); i++) {
  		vector<int> triangle;
  		for (unsigned int j=0; j<triangles_u[i].size(); j++)
  			triangle.push_back(vuindices[triangles_u[i][j]]);
  		tempidx.push_back(triangle);
  	}

  	triangles_u.clear();
  	for (unsigned int i=0; i<tempidx.size(); i++)
  		triangles_u.push_back(tempidx[i]);

    for (unsigned int i=0; i<triangles_u.size(); i++) {
        if (triangles_u[i].size() < 3)
            continue;

        vector<int> edge1;
        vector<int> edge2;
        vector<int> edge3;
        edge1.push_back(triangles_u[i][0]);
        edge1.push_back(triangles_u[i][1]);
        edge2.push_back(triangles_u[i][0]);
        edge2.push_back(triangles_u[i][2]);
        edge3.push_back(triangles_u[i][1]);
        edge3.push_back(triangles_u[i][2]);

        vector<vector<int> > tempedges;
        tempedges.push_back(edge1);
        tempedges.push_back(edge2);
        tempedges.push_back(edge3);

        edgespertriangle_u.push_back(tempedges);
    }

	vector<vector<vector<int> > > vCommonEdges;
	vCommonEdges.resize(triangles_u.size());
	for (unsigned int i=0; i<vCommonEdges.size(); i++)
		vCommonEdges[i].resize(3);

	vIdxNewVertices_u.clear();
	vIdxNewVertices_u.resize(triangles_u.size());
	vector<bool> vbComputedTriangle = vector<bool>(triangles_u.size(), false);
	quads_u.clear();

    for (unsigned int i=0; i<edgespertriangle_u.size(); i++) {
        if (edgespertriangle_u[i].size()<3)
            continue;

        for (unsigned int j=0; j<edgespertriangle_u.size(); j++) {
            if (j==i)
                continue;

            if (edgespertriangle_u[j].size()<3)
                continue;

            for (unsigned int k=0; k<edgespertriangle_u[i].size(); k++) {
                if  (
                        ( (edgespertriangle_u[i][0][0] == edgespertriangle_u[j][k][0])
                       && (edgespertriangle_u[i][0][1] == edgespertriangle_u[j][k][1]) )
                        ||
                        ( (edgespertriangle_u[i][0][0] == edgespertriangle_u[j][k][1])
                       && (edgespertriangle_u[i][0][1] == edgespertriangle_u[j][k][0]) )
                    )
                    {
						vCommonEdges[i][0].push_back(0);
						vCommonEdges[i][0].push_back(j);
						vCommonEdges[i][0].push_back(k);
                    }

                if  (
                        ( (edgespertriangle_u[i][1][0] == edgespertriangle_u[j][k][0])
                       && (edgespertriangle_u[i][1][1] == edgespertriangle_u[j][k][1]) )
                        ||
                        ( (edgespertriangle_u[i][1][0] == edgespertriangle_u[j][k][1])
                       && (edgespertriangle_u[i][1][1] == edgespertriangle_u[j][k][0]) )
                    )
                    {
						vCommonEdges[i][1].push_back(1);
						vCommonEdges[i][1].push_back(j);
						vCommonEdges[i][1].push_back(k);
                    }

                if  (
                        ( (edgespertriangle_u[i][2][0] == edgespertriangle_u[j][k][0])
                       && (edgespertriangle_u[i][2][1] == edgespertriangle_u[j][k][1]) )
                        ||
                        ( (edgespertriangle_u[i][2][0] == edgespertriangle_u[j][k][1])
                       && (edgespertriangle_u[i][2][1] == edgespertriangle_u[j][k][0]) )
                    )
                    {
						vCommonEdges[i][2].push_back(2);
						vCommonEdges[i][2].push_back(j);
						vCommonEdges[i][2].push_back(k);
                    }
            }
        }
    }


    vNewPointsBase_u.clear();

	for (unsigned int i=0; i<triangles_u.size(); i++) {
		int v0idx = triangles_u[i][0];
		vector<float> v0 = vMPsXYZN_ut[v0idx];
		int v1idx = triangles_u[i][1];
		vector<float> v1 = vMPsXYZN_ut[v1idx];
		int v2idx = triangles_u[i][2];
		vector<float> v2 = vMPsXYZN_ut[v2idx];

		vector<float> m01;
		int m01idx;
		vector<float> m02;
		int m02idx;
		vector<float> m12;
		int m12idx;
		vector<float> g;
		int gidx;
		vector<int> quad1;
		vector<int> quad2;
		vector<int> quad3;

		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 0) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m01idx = vIdxNewVertices_u[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m01 = vMPsXYZN_ut[m01idx];
					}
					else {
						m01idx = vMPsXYZN_ut.size();

						m01.push_back((v0[0]+v1[0])/2);
		 				m01.push_back((v0[1]+v1[1])/2);
		  				m01.push_back((v0[2]+v1[2])/2);

						vMPsXYZN_ut.push_back(m01);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v1idx);
						vNewPointsBase_u.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m01idx = vMPsXYZN_ut.size();

				m01.push_back((v0[0]+v1[0])/2);
		 		m01.push_back((v0[1]+v1[1])/2);
		  		m01.push_back((v0[2]+v1[2])/2);

				vMPsXYZN_ut.push_back(m01);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v1idx);
                vNewPointsBase_u.push_back(NewPointsBasei);
			}
		}

		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 1) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m02idx = vIdxNewVertices_u[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m02 = vMPsXYZN_ut[m02idx];
					}
					else {
						m02idx = vMPsXYZN_ut.size();

					  	m02.push_back((v0[0]+v2[0])/2);
					  	m02.push_back((v0[1]+v2[1])/2);
					  	m02.push_back((v0[2]+v2[2])/2);

						vMPsXYZN_ut.push_back(m02);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase_u.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m02idx = vMPsXYZN_ut.size();

			 	m02.push_back((v0[0]+v2[0])/2);
				m02.push_back((v0[1]+v2[1])/2);
				m02.push_back((v0[2]+v2[2])/2);

				vMPsXYZN_ut.push_back(m02);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase_u.push_back(NewPointsBasei);
			}
		}

		for (unsigned int j=0; j<vCommonEdges[i].size(); j++) {
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty()) {
				if (vCommonEdges[i][j][0] == 2) {
					if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
						m12idx = vIdxNewVertices_u[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m12 = vMPsXYZN_ut[m12idx];
					}
					else {
						m12idx = vMPsXYZN_ut.size();

					  	m12.push_back((v1[0]+v2[0])/2);
					  	m12.push_back((v1[1]+v2[1])/2);
					  	m12.push_back((v1[2]+v2[2])/2);

						vMPsXYZN_ut.push_back(m12);

						NewPointsBasei.push_back(v1idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase_u.push_back(NewPointsBasei);
					}
				}
			}
			else {
				m12idx = vMPsXYZN_ut.size();

				m12.push_back((v1[0]+v2[0])/2);
				m12.push_back((v1[1]+v2[1])/2);
				m12.push_back((v1[2]+v2[2])/2);

				vMPsXYZN_ut.push_back(m12);

				NewPointsBasei.push_back(v1idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase_u.push_back(NewPointsBasei);
			}
		}

		g.push_back((v0[0]+v1[0]+v2[0])/3);
		g.push_back((v0[1]+v1[1]+v2[1])/3);
		g.push_back((v0[2]+v1[2]+v2[2])/3);
		gidx = vMPsXYZN_ut.size();
		vMPsXYZN_ut.push_back(g);

		vector<int> NewPointsBasei;
        NewPointsBasei.push_back(v0idx);
        NewPointsBasei.push_back(v1idx);
        NewPointsBasei.push_back(v2idx);
        vNewPointsBase_u.push_back(NewPointsBasei);

		quad1.push_back(v0idx);
		quad1.push_back(m01idx);
		quad1.push_back(gidx);
		quad1.push_back(m02idx);

		quad2.push_back(v1idx);
		quad2.push_back(m12idx);
		quad2.push_back(gidx);
		quad2.push_back(m01idx);

		quad3.push_back(v2idx);
		quad3.push_back(m02idx);
		quad3.push_back(gidx);
		quad3.push_back(m12idx);

		quads_u.push_back(quad1);
		quads_u.push_back(quad2);
		quads_u.push_back(quad3);

		vIdxNewVertices_u[i].push_back(m01idx);
		vIdxNewVertices_u[i].push_back(m02idx);
		vIdxNewVertices_u[i].push_back(m12idx);
		vIdxNewVertices_u[i].push_back(gidx);

		vbComputedTriangle[i] = true;
	}
}


void FEM::SetSecondLayer(int nMode) {
    if (nMode == 1) {
        for (unsigned int i=0; i<vMPsXYZN_t.size(); i++) {
            vector<float> point;
            float pos1 = vMPsXYZN_t[i][0] - h;
            float pos2 = vMPsXYZN_t[i][1] - h;
            float pos3 = vMPsXYZN_t[i][2] - h;
            point.push_back(pos1);
            point.push_back(pos2);
            point.push_back(pos3);
            vMPsXYZN_t2.push_back(point);
        }
    }
    else if (nMode == 2) {
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++) {
            vector<float> point;
            float pos1 = vMPsXYZN_ut[i][0] - h;
            float pos2 = vMPsXYZN_ut[i][1] - h;
            float pos3 = vMPsXYZN_ut[i][2] - h;
            point.push_back(pos1);
            point.push_back(pos2);
            point.push_back(pos3);
            vMPsXYZN_ut2.push_back(point);
        }
    }
}


void FEM::Set_u0(vector<MapPoint*> vpMPs, int nMode) {
    if (nMode == 1 ) {
        u0.clear();
        for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
            for (unsigned int j=0; j<3; j++)
                u0.push_back(vMPsXYZN_t[i][j]);
        for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
            for (unsigned int j=0; j<3; j++)
                u0.push_back(vMPsXYZN_t2[i][j]);
    }
    else if (nMode == 2) {
        u0u.clear();
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++)
            for (unsigned int j=0; j<3; j++)
                u0u.push_back(vMPsXYZN_ut[i][j]);
        for (unsigned int i=0; i<vMPsXYZN_ut2.size(); i++)
            for (unsigned int j=0; j<3; j++)
                u0u.push_back(vMPsXYZN_ut2[i][j]);
    }
}


vector<vector<float> > FEM::ComputeKei(vector<vector<float> > vfPts) {
    vector<vector<float> > vBtDB;
    vector<float> vBtDBi = vector<float>(24,0.0);
        for (unsigned int i=0; i<24; i++)
            vBtDB.push_back(vBtDBi);

    for (unsigned int ops=0; ops<gs.size(); ops++) {
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


bool FEM::MatrixAssembly(int nMode) {
    // To easily insert the nodes into the stiffness matrix, the node indexes are used to
    // set their position in the matrix. NULL positions will be set to 0 by default.
    // The size of the stiffness matrix equals the size of the MPs vector, the positions
    // corresponding to a offvertex are set to zero in the stiffness matrix

    if (nMode == 1) {
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
    else if (nMode == 2) {
        int nTotalNodes = vMPsXYZN_ut.size() + vMPsXYZN_ut2.size();
        Kusize = 3*nTotalNodes;
        cout << "                - MatAssembly (hex, 24 DoF per el. Curr: " << quads_u.size() << " quads, " << nTotalNodes << " nodes, Kusize = " << Kusize << ")" << endl;

        if (Kusize<=3)
            return false;

        Ku.clear();
        vector<float> ki = vector<float>(Kusize,0.0);
        for (unsigned int i=0; i<Kusize; i++)
            Ku.push_back(ki);

        for (unsigned int i=0; i<quads_u.size(); i++) {
            vector<int> nodes;
            nodes.push_back(quads_u[i][0]);
            nodes.push_back(quads_u[i][1]);
            nodes.push_back(quads_u[i][2]);
            nodes.push_back(quads_u[i][3]);
            nodes.push_back(quads_u[i][0] + vMPsXYZN_ut.size());
            nodes.push_back(quads_u[i][1] + vMPsXYZN_ut.size());
            nodes.push_back(quads_u[i][2] + vMPsXYZN_ut.size());
            nodes.push_back(quads_u[i][3] + vMPsXYZN_ut.size());

            vector<vector<float> > vfPts;
            for (unsigned int j=0; j<4; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][2]);
                vfPts.push_back(vfPtsi);
            }

            for (unsigned int j=0; j<4; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][2]);
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
                            Ku[matpos[node1]+m][matpos[node2]+n] += vBtDB[node1*3+m][node2*3+n];
        }
        return true;
    }
}


/*
float FEM::ComputeDet(vector<vector<float> > Kin) {
    float fs = Kin.size();

    for (unsigned int i=0; i<fs; i++);





}
*/


bool FEM::ComputeK1() {
    // Método Gauss Jordan para cálculo de matriz inversa (A|I) -> (I|A1)
    cout << "                - Compute K inverse " << endl;
    vector<vector<float> > Ktemp = Ku;
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
    K1 = vector<vector<float> >(Kusize,vector<float>(Kusize,(0.0)));
    for (unsigned int i=0; i<K1.size(); i++)
        K1[i][i] = 1.0;

    for (unsigned int i=0; i<Kusize; i++) {
        float fPivot = Ktemp[i][i];

        if (fPivot==0) {
            unsigned int nRowObj = i;
            if (i<(Kusize-1))
                for (unsigned int j=i+1; j<Kusize; j++) {
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

            for (unsigned int aa=0; aa<Kusize; aa++) {
                Ktemp[i][aa] = Ktemp[nRowObj][aa];
                Ktemp[nRowObj][aa] = vfAux[aa];
            }
        }

        fPivot = Ktemp[i][i];

        if (fPivot!=1)
            for (unsigned int j=0; j<Kusize; j++) {
                if (Ktemp[i][j]!=0)
                    Ktemp[i][j] /= fPivot;
                if (K1[i][j]!=0)
                    K1[i][j] /= fPivot;
            }


        for (unsigned int j=0; j<Kusize; j++) {
            if (j==i)
                continue;
            float fPivotRow = Ktemp[j][i];
            for (unsigned int s=0; s<Kusize; s++) {
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





















/// Old functions

bool FEM::LoadMPsIntoCloud() {
    bool b1 = LoadTrackedMPsIntoCloud();
    bool b2 = LoadUnTrackedMPsIntoCloud();

    if (b1&&b2)
        return true;
    else
        return false;
}

bool FEM::LoadTrackedMPsIntoCloud() {
    // Load top layer into cloud
    pc_t_0.width = vMPsXYZN_t.size(); //cloudTracked.width    = 5;
    pc_t_0.height   = 1;
    pc_t_0.is_dense = false;
    pc_t_0.points.resize (pc_t_0.width * pc_t_0.height);

    //unsigned int sizetoplayer = pc_t_0.points.size();

    //for (size_t i = 0; i < sizetoplayer; i++)
    for (size_t i = 0; i < pc_t_0.points.size(); i++)
    {
        pc_t_0.points[i].x = vMPsXYZN_t[i][0];
        pc_t_0.points[i].y = vMPsXYZN_t[i][1];
        pc_t_0.points[i].z = vMPsXYZN_t[i][2];
    }

    if (pc_t_0.width>0)
    {
        pcPath_t = pcPath_t + std::to_string(nFrameId) + ".pcd";
        pcl::io::savePCDFileASCII (pcPath_t, pc_t_0);
        return true;
    }
    else
        return false;
}

bool FEM::LoadUnTrackedMPsIntoCloud() {
    cloudUnTracked.width = vMPsXYZN_t.size() + vMPsXYZN_ut.size(); //cloudTracked.width    = 5;
    cloudUnTracked.height   = 1;
    cloudUnTracked.is_dense = false;
    cloudUnTracked.points.resize (cloudUnTracked.width * cloudUnTracked.height);

    for (size_t i = 0; i < vMPsXYZN_t.size(); ++i)
    {
        cloudUnTracked.points[i].x = vMPsXYZN_t[i][0];
        cloudUnTracked.points[i].y = vMPsXYZN_t[i][1];
        cloudUnTracked.points[i].z = vMPsXYZN_t[i][2];
    }
    for (size_t i = vMPsXYZN_t.size(), j = 0 ; i < vMPsXYZN_t.size() + vMPsXYZN_ut.size(); ++i, ++j)
    {
        cloudUnTracked.points[i].x = vMPsXYZN_ut[j][0];
        cloudUnTracked.points[i].y = vMPsXYZN_ut[j][1];
        cloudUnTracked.points[i].z = vMPsXYZN_ut[j][2];
    }

    cout << "pcut" << endl;
    if (cloudUnTracked.width>0)
    {
        pcPath_ut = pcPath_ut + std::to_string(nFrameId) + ".pcd";
        pcl::io::savePCDFileASCII (pcPath_ut, cloudUnTracked);
        return true;
    }
    else
        return false;
}

bool FEM::MLS() {
    bool b1 = MLS_t();
    //bool b2 = MLS_ut();

    if (b1)
        return true;
    else
        return false;
}

bool FEM::MLS_t() {

    PointCloud<PointXYZ>::Ptr pc_t_mls_0 (new pcl::PointCloud<PointXYZ> (pc_t_0));

    /// MLS 1st iteration

    // Create a KD-Tree
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree1 (new pcl::search::KdTree<pcl::PointXYZ>);

    // Set parameters
    pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls1;
    mls1.setComputeNormals (true);
    mls1.setInputCloud (pc_t_mls_0);
    mls1.setPolynomialFit (true);
    mls1.setSearchMethod (tree1);
    mls1.setSearchRadius (900000);	//100 //0.25    // Original 0.03 // Default 0.00
    mls1.setPolynomialOrder(3); //50 // Default 2

    mls1.process (pc_t_1);
    pcl::PointCloud<pcl::PointNormal>::Ptr ptr_pc_t_mls_1 (new pcl::PointCloud<pcl::PointNormal> (pc_t_1));

    // Get corresponding indexes: for each output point, returns the index of the input one.
    PointIndicesPtr pIdx1 = mls1.getCorrespondingIndices();
    for (unsigned int i=0; i<pIdx1->indices.size(); i++)
        vtindices.push_back(pIdx1->indices[i]);

    int idxit = 0;
    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
    {
        int currentpos = i;
        if (currentpos==vtindices[idxit])
        {
            vbtMPsActive[i] = true;
            idxit++;
        }
        else
            vbtMPsActive[i] = false;
    }

    /// Save smoothed cloud
    if (pc_t_1.width>0)
    {
        pcsPath_t = pcsPath_t + std::to_string(nFrameId) + ".pcd";
        pcl::io::savePCDFileASCII (pcsPath_t, pc_t_1);
        return true;
    }
    else
        return false;
}

bool FEM::MLS_ut() {
    pcl::PointCloud<pcl::PointXYZ>::Ptr slCloud (new pcl::PointCloud<pcl::PointXYZ> ());
    pcl::io::loadPCDFile (pcPath_ut, *slCloud);

    if (slCloud->width == 0)
        return false;

    // Create a KD-Tree
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;

    // Set parameters
    mls.setComputeNormals (true);
    mls.setInputCloud (slCloud);
    mls.setPolynomialFit (true);
    mls.setSearchMethod (tree);
    mls.setSearchRadius (0.03);

    // Apply
    mls.process (sCloud_ut);

    // Save smoothed cloud
    pcsPath_ut += std::to_string(nFrameId);
    pcsPath_ut += ".pcd";

    if (sCloud_ut.width==0)
    {
        cout << "No data to save" << endl;
        return false;
    }
    else
    {
        pcl::io::savePCDFileASCII (pcsPath_ut, sCloud_ut);
        return true;
    }
}

bool FEM::ComputeMesh() {
    bool b3 = ComputeMesh_t();

    if (b3)
        return true;
    else
        return false;
}

bool FEM::ComputeMesh_t() {
    // TRIANGULATION, GREEDY PROJECTION

    // Cloud
    PointCloud<PointNormal>::Ptr ptr_pc_t_mesh_1 (new pcl::PointCloud<PointNormal> (pc_t_1));

  	// Search tree
  	pcl::search::KdTree<pcl::PointNormal>::Ptr tree (new pcl::search::KdTree<pcl::PointNormal>);
  	tree->setInputCloud (ptr_pc_t_mesh_1);

    // Greedy Projection
    int nInMu = 5;
    int nSearchRad = 30;
    int nMaxNeig = 70;
    int nSurAng = 150;
    int nMinAng = 45;
    int nMaxAng = 90;
    CalculateGP3Parameters(ptr_pc_t_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);

    static pcl::GreedyProjectionTriangulation<pcl::PointNormal> GreedyProj3;
    GreedyProj3.setMu (nInMu);    // original 2.8
    GreedyProj3.setSearchRadius(nSearchRad);      // original 0.50
    GreedyProj3.setMaximumNearestNeighbors (nMaxNeig);        // original 150
    GreedyProj3.setMaximumSurfaceAngle(nSurAng); // 45 degrees
    GreedyProj3.setMinimumAngle(nMinAng); // 10 degrees
    GreedyProj3.setMaximumAngle(nMaxAng); // 120 degrees 2/3
    GreedyProj3.setNormalConsistency(true);

    cout << "GP3 Params: " << nInMu << " " << nSearchRad << " " << nMaxNeig << endl;

    bool bGo = true;
    if (bGo == true)
    {
        GreedyProj3.setMu (5.5);    // original 2.8
        GreedyProj3.setSearchRadius(0.75);      // original 0.50
        GreedyProj3.setMaximumNearestNeighbors (50);        // original 150
        GreedyProj3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
        GreedyProj3.setMinimumAngle(M_PI/18); // 10 degrees
        GreedyProj3.setMaximumAngle(9*M_PI/10); // 120 degrees 2/3
        GreedyProj3.setNormalConsistency(true);
    }



    GreedyProj3.setInputCloud (ptr_pc_t_mesh_1);
    GreedyProj3.setSearchMethod (tree);
    GreedyProj3.reconstruct (mesh_t);

    if (mesh_t.polygons.size()>0)
    {
        // Get the vertex indexes of each triangle
        for (unsigned int i=0; i<mesh_t.polygons.size(); i++)
        {
            unsigned int nver0 = mesh_t.polygons[i].vertices[0];
            unsigned int nver1 = mesh_t.polygons[i].vertices[1];
            unsigned int nver2 = mesh_t.polygons[i].vertices[2];

            vector<int> triangle;
            triangle.push_back(nver0);
            triangle.push_back(nver1);
            triangle.push_back(nver2);
            triangles_t.push_back(triangle);

            //vector<vector<MapPoint*> > vpMPs2Draw;
            vector<MapPoint*> vpMP2Draw;
            vpMP2Draw.push_back(vpMPs_t[vtindices[nver0]]);
            vpMP2Draw.push_back(vpMPs_t[vtindices[nver1]]);
            vpMP2Draw.push_back(vpMPs_t[vtindices[nver2]]);
            vpMPs2Draw.push_back(vpMP2Draw);

            vector<cv::KeyPoint*> vpKP2Draw;
            vpKP2Draw.push_back(vpKPs_t[vtindices[nver0]]);
            vpKP2Draw.push_back(vpKPs_t[vtindices[nver1]]);
            vpKP2Draw.push_back(vpKPs_t[vtindices[nver2]]);
            vpKPs2Draw.push_back(vpKP2Draw);
        }

        // Save reconstruction
        pcrPath_t = pcrPath_t + std::to_string(nFrameId) + ".vtk";
        pcl::io::saveVTKFile (pcrPath_t, mesh_t);
        return true;
    }
    else
    {
        cout << "No data" << endl;
        return false;
    }
}


/*
void FEM::tri2quad() {
    tri2quad_t();
    tri2quad_ut();
}
*/

/*
void FEM::tri2quad_t()
{
  	// Get all vertices in mesh
  	//cout << triangles_t.size() << endl;
  	vector<int> vertices;
  	for (unsigned int i=0; i<triangles_t.size(); i++)
  		for (unsigned int j=0; j<triangles_t[i].size(); j++)
  			if ( !( std::find(vertices.begin(), vertices.end(), triangles_t[i][j]) != vertices.end() ) )
            	vertices.push_back(triangles_t[i][j]);

	int maxvertex = vertices.back();

  	vector<int> offvertices;
  	for (int i=0; i<maxvertex; i++)
  		if ( !( std::find(vertices.begin(), vertices.end(), i) != vertices.end() ) )
            offvertices.push_back(i);

    for (unsigned int i=0; i<vertices.size(); i++)
    {
    	vbtMPsActive[vtindices[vertices[i]]] = true;
    }

    for (unsigned int i=0; i<offvertices.size(); i++)
    {
    	vbtMPsActive[vtindices[offvertices[i]]] = false;
    }

  	// Get true indexes
  	vector<vector<int> > tempidx;
  	for (unsigned int i=0; i<triangles_t.size(); i++)
  	{
  		vector<int> triangle;
  		for (unsigned int j=0; j<triangles_t[i].size(); j++)
  			triangle.push_back(vtindices[triangles_t[i][j]]);
  		tempidx.push_back(triangle);
  	}

  	triangles_t.clear();
  	for (unsigned int i=0; i<tempidx.size(); i++)
  		triangles_t.push_back(tempidx[i]);


  	// Record the 3 edges_t in each triangle
    for (unsigned int i=0; i<triangles_t.size(); i++)
    {
        if (triangles_t[i].size() < 3)
            continue;

        vector<int> edge1;
        vector<int> edge2;
        vector<int> edge3;
        edge1.push_back(triangles_t[i][0]);
        edge1.push_back(triangles_t[i][1]);
        edge2.push_back(triangles_t[i][0]);
        edge2.push_back(triangles_t[i][2]);
        edge3.push_back(triangles_t[i][1]);
        edge3.push_back(triangles_t[i][2]);

        vector<vector<int> > tempedges;
        tempedges.push_back(edge1);
        tempedges.push_back(edge2);
        tempedges.push_back(edge3);

        edgespertriangle_t.push_back(tempedges);
    }



	vector<vector<vector<int> > > vCommonEdges;
	vCommonEdges.resize(triangles_t.size());
	for (unsigned int i=0; i<vCommonEdges.size(); i++)
		vCommonEdges[i].resize(3);
		//triangulo1
		//	commonedge1
		//		edge1 triangulo2 edge2
		//	commonedge2
		//		edge1 triangulo2 edge2
		//	commonedge3
		//		edge1 triangulo2 edge2
		//triangulo2
		//triangulo3
		//...

	//vector<vector<int> > vIdxNewVertices;
	vIdxNewVertices.clear();
	vIdxNewVertices.resize(triangles_t.size());
	//for(unsigned int i=0; i<vIdxNewVertices.size(); i++)
		//for (unsigned int j=0; j<4; j++)
			//vIdxNewVertices[i].push_back(-1);
		//triangulo1
		//	idxvertexedge0 idxvertexedge1 idxvertexedge2 idxbaricentro
		//triangulo2
		//triangulo3
		//...

	vector<bool> vbComputedTriangle = vector<bool>(triangles_t.size(), false);
		//triangulo1
		//triangulo2
		//triangulo3
		//...

	quads_t.clear();

	// Fill vCommonEdges with the function originally designed to join common edges.
    for (unsigned int i=0; i<edgespertriangle_t.size(); i++)
    {
        // Won't happen, but just in case
        if (edgespertriangle_t[i].size()<3)
            continue;

        // For each triangle, we look for common edges with all the other polygons
        for (unsigned int j=0; j<edgespertriangle_t.size(); j++)
        {
            // Don't look within itself
            if (j==i)
                continue;

            // Again, won't happen, but just in case
            if (edgespertriangle_t[j].size()<3)
                continue;

            // For all the edges within the reference triangle
            for (unsigned int k=0; k<edgespertriangle_t[i].size(); k++)
            {
                // Compare the first edge in the reference triangle with all the three(k)
                // edges in the candidate
                if
                    (
                        ( (edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 0 & 1
						vCommonEdges[i][0].push_back(0);
						vCommonEdges[i][0].push_back(j);
						vCommonEdges[i][0].push_back(k);
                    }

                // Compare the second edge in the reference triangle with all the three(k)
                // edges in the candidate
                if
                    (
                        ( (edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 0 & 2
						vCommonEdges[i][1].push_back(1);
						vCommonEdges[i][1].push_back(j);
						vCommonEdges[i][1].push_back(k);
                    }

                // Compare the third edge in the reference triangle with all the three(k)
                // edges in the candidate
                if
                    (
                        ( (edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][0])
                       && (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][1]) )
                        ||
                        ( (edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][1])
                       && (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][0]) )
                    )
                    {
                        // edge1 = 1 & 2
						vCommonEdges[i][2].push_back(2);
						vCommonEdges[i][2].push_back(j);
						vCommonEdges[i][2].push_back(k);
                    }
                    // k changes here, so in each loop, we compare all the edges in
                    // the reference triangle with only one edge of the candidate
            }
        }   // End of the loop through all the candidates for a reference triangle
    }   // End of the loop through all the triangles set as reference

    //cout << "edgespertrianglecomputed" << endl;

    vNewPointsBase.clear();

	for (unsigned int i=0; i<triangles_t.size(); i++)
	{
		// Get vertices and indices for refferece triangle
		int v0idx = triangles_t[i][0];
		vector<float> v0 = vMPsXYZN_t[v0idx];

		int v1idx = triangles_t[i][1];
		vector<float> v1 = vMPsXYZN_t[v1idx];

		int v2idx = triangles_t[i][2];
		vector<float> v2 = vMPsXYZN_t[v2idx];

		// Set variables to store the ressult
		vector<float> m01;
		int m01idx;
		vector<float> m02;
		int m02idx;
		vector<float> m12;
		int m12idx;
		vector<float> g;
		int gidx;
		vector<int> quad1;
		vector<int> quad2;
		vector<int> quad3;

		// Compute edge 0
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++)
		{
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty())
			{
				if (vCommonEdges[i][j][0] == 0)
				{
					if (vbComputedTriangle[vCommonEdges[i][j][1]])
					{
						m01idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m01 = vMPsXYZN_t[m01idx];
					}
					else
					{
						m01idx = vMPsXYZN_t.size();

						m01.push_back((v0[0]+v1[0])/2);
		 				m01.push_back((v0[1]+v1[1])/2);
		  				m01.push_back((v0[2]+v1[2])/2);

						vMPsXYZN_t.push_back(m01);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v1idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else
			{
				m01idx = vMPsXYZN_t.size();

				m01.push_back((v0[0]+v1[0])/2);
		 		m01.push_back((v0[1]+v1[1])/2);
		  		m01.push_back((v0[2]+v1[2])/2);

				vMPsXYZN_t.push_back(m01);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v1idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute edge 1
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++)
		{
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty())
			{
				if (vCommonEdges[i][j][0] == 1)
				{
					if (vbComputedTriangle[vCommonEdges[i][j][1]])
					{
						m02idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m02 = vMPsXYZN_t[m02idx];
					}
					else
					{
						m02idx = vMPsXYZN_t.size();

					  	m02.push_back((v0[0]+v2[0])/2);
					  	m02.push_back((v0[1]+v2[1])/2);
					  	m02.push_back((v0[2]+v2[2])/2);

						vMPsXYZN_t.push_back(m02);

						NewPointsBasei.push_back(v0idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else
			{
				m02idx = vMPsXYZN_t.size();

			 	m02.push_back((v0[0]+v2[0])/2);
				m02.push_back((v0[1]+v2[1])/2);
				m02.push_back((v0[2]+v2[2])/2);

				vMPsXYZN_t.push_back(m02);

				NewPointsBasei.push_back(v0idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute edge 2
		for (unsigned int j=0; j<vCommonEdges[i].size(); j++)
		{
		    vector<int> NewPointsBasei;
			if (!vCommonEdges[i][j].empty())
			{
				if (vCommonEdges[i][j][0] == 2)
				{
					if (vbComputedTriangle[vCommonEdges[i][j][1]])
					{
						m12idx = vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
						m12 = vMPsXYZN_t[m12idx];
					}
					else
					{
						m12idx = vMPsXYZN_t.size();

					  	m12.push_back((v1[0]+v2[0])/2);
					  	m12.push_back((v1[1]+v2[1])/2);
					  	m12.push_back((v1[2]+v2[2])/2);

						vMPsXYZN_t.push_back(m12);

						NewPointsBasei.push_back(v1idx);
						NewPointsBasei.push_back(v2idx);
						vNewPointsBase.push_back(NewPointsBasei);
					}
				}
			}
			else
			{
				m12idx = vMPsXYZN_t.size();

				m12.push_back((v1[0]+v2[0])/2);
				m12.push_back((v1[1]+v2[1])/2);
				m12.push_back((v1[2]+v2[2])/2);

				vMPsXYZN_t.push_back(m12);

				NewPointsBasei.push_back(v1idx);
                NewPointsBasei.push_back(v2idx);
                vNewPointsBase.push_back(NewPointsBasei);
			}
		}

		// Compute barycentre in each triangle
		g.push_back((v0[0]+v1[0]+v2[0])/3);
		g.push_back((v0[1]+v1[1]+v2[1])/3);
		g.push_back((v0[2]+v1[2]+v2[2])/3);
		gidx = vMPsXYZN_t.size();
		vMPsXYZN_t.push_back(g);

		vector<int> NewPointsBasei;
        NewPointsBasei.push_back(v0idx);
        NewPointsBasei.push_back(v1idx);
        NewPointsBasei.push_back(v2idx);
        vNewPointsBase.push_back(NewPointsBasei);

		// Set the three quads
		quad1.push_back(v0idx);
		quad1.push_back(m01idx);
		quad1.push_back(gidx);
		quad1.push_back(m02idx);

		quad2.push_back(v1idx);
		quad2.push_back(m12idx);
		quad2.push_back(gidx);
		quad2.push_back(m01idx);

		quad3.push_back(v2idx);
		quad3.push_back(m02idx);
		quad3.push_back(gidx);
		quad3.push_back(m12idx);

		// Add the quads to the global vector
		quads_t.push_back(quad1);
		quads_t.push_back(quad2);
		quads_t.push_back(quad3);

		// Set vertices idx
		vIdxNewVertices[i].push_back(m01idx);
		vIdxNewVertices[i].push_back(m02idx);
		vIdxNewVertices[i].push_back(m12idx);
		vIdxNewVertices[i].push_back(gidx);

		// Set triangle as already computed
		vbComputedTriangle[i] = true;

	}
}

*/


void FEM::tri2quad_ut() {
    for (unsigned int i=0; i<edgespertriangle_ut.size(); i++)
    {
        if (edgespertriangle_ut[i].size()<3)
            continue;

        for (unsigned int j=0; j<edgespertriangle_ut.size(); j++)
        {
            if (j==i)
                continue;

            if (edgespertriangle_ut[j].size()<3)
                continue;

            for (unsigned int k=0; k<edgespertriangle_ut[i].size(); k++)
            {
                if
                    (
                        ( (edgespertriangle_ut[i][0][0] == edgespertriangle_ut[j][k][0])
                       && (edgespertriangle_ut[i][0][1] == edgespertriangle_ut[j][k][1]) )
                        ||
                        ( (edgespertriangle_ut[i][0][0] == edgespertriangle_ut[j][k][1])
                       && (edgespertriangle_ut[i][0][1] == edgespertriangle_ut[j][k][0]) )
                    )
                    {
                        if ( (edgespertriangle_ut[i][0][0] == 0) && (edgespertriangle_ut[i][0][1] == 0) )
                            continue;
                        if ( (edgespertriangle_ut[j][0][0] == 0) && (edgespertriangle_ut[j][0][1] == 0) )
                            continue;

                        // edge1 = 0 & 1
                        int finalvertex = 0;
                        if (k==0) finalvertex = edgespertriangle_ut[j][1][1];
                        else if (k==1) finalvertex = edgespertriangle_ut[j][0][1];
                        else if (k==2) finalvertex = edgespertriangle_ut[j][0][0];
                        vector<int> quad;
                        quad.push_back(edgespertriangle_ut[i][1][1]);
                        quad.push_back(edgespertriangle_ut[i][0][0]);
                        quad.push_back(finalvertex);
                        quad.push_back(edgespertriangle_ut[i][0][1]);
                        quads_ut.push_back(quad);

                        for (unsigned int l=0; l<edgespertriangle_ut[i].size(); l++)
                        {
                            edgespertriangle_ut[i][l][0] = 0;
                            edgespertriangle_ut[i][l][1] = 0;
                            edgespertriangle_ut[j][l][0] = 0;
                            edgespertriangle_ut[j][l][1] = 0;
                        }

                        break;
                    }
                if
                    (
                        ( (edgespertriangle_ut[i][1][0] == edgespertriangle_ut[j][k][0])
                       && (edgespertriangle_ut[i][1][1] == edgespertriangle_ut[j][k][1]) )
                        ||
                        ( (edgespertriangle_ut[i][1][0] == edgespertriangle_ut[j][k][1])
                       && (edgespertriangle_ut[i][1][1] == edgespertriangle_ut[j][k][0]) )
                    )
                    {
                        if ( (edgespertriangle_ut[i][1][0] == 0) && (edgespertriangle_ut[i][1][1] == 0) )
                            continue;
                        if ( (edgespertriangle_ut[j][1][0] == 0) && (edgespertriangle_ut[j][1][1] == 0) )
                            continue;

                        // edge1 = 0 & 2
                        int finalvertex = 0;
                        if (k==0) finalvertex = edgespertriangle_ut[j][1][1];
                        else if (k==1) finalvertex = edgespertriangle_ut[j][0][1];
                        else if (k==2) finalvertex = edgespertriangle_ut[j][0][0];
                        vector<int> quad;
                        quad.push_back(edgespertriangle_ut[i][0][1]);
                        quad.push_back(edgespertriangle_ut[i][1][0]);
                        quad.push_back(finalvertex);
                        quad.push_back(edgespertriangle_ut[i][1][1]);
                        quads_ut.push_back(quad);

                        for (unsigned int l=0; l<edgespertriangle_ut[i].size(); l++)
                        {
                            edgespertriangle_ut[i][l][0] = 0;
                            edgespertriangle_ut[i][l][1] = 0;
                            edgespertriangle_ut[j][l][0] = 0;
                            edgespertriangle_ut[j][l][1] = 0;
                        }

                        break;
                    }
                if
                    (
                        ( (edgespertriangle_ut[i][2][0] == edgespertriangle_ut[j][k][0])
                       && (edgespertriangle_ut[i][2][1] == edgespertriangle_ut[j][k][1]) )
                        ||
                        ( (edgespertriangle_ut[i][2][0] == edgespertriangle_ut[j][k][1])
                       && (edgespertriangle_ut[i][2][1] == edgespertriangle_ut[j][k][0]) )
                    )
                    {
                        if ( (edgespertriangle_ut[i][2][0] == 0) && (edgespertriangle_ut[i][2][1] == 0) )
                            continue;
                        if ( (edgespertriangle_ut[j][2][0] == 0) && (edgespertriangle_ut[j][2][1] == 0) )
                            continue;

                        // edge1 = 1 & 2
                        int finalvertex = 0;
                        if (k==0) finalvertex = edgespertriangle_ut[j][1][1];
                        else if (k==1) finalvertex = edgespertriangle_ut[j][0][1];
                        else if (k==2) finalvertex = edgespertriangle_ut[j][0][0];
                        vector<int> quad;
                        quad.push_back(edgespertriangle_ut[i][0][0]);
                        quad.push_back(edgespertriangle_ut[i][2][0]);
                        quad.push_back(finalvertex);
                        quad.push_back(edgespertriangle_ut[i][2][1]);
                        quads_ut.push_back(quad);

                        for (unsigned int l=0; l<edgespertriangle_ut[i].size(); l++)
                        {
                            edgespertriangle_ut[i][l][0] = 0;
                            edgespertriangle_ut[i][l][1] = 0;
                            edgespertriangle_ut[j][l][0] = 0;
                            edgespertriangle_ut[j][l][1] = 0;
                        }

                        break;
                    }
            }
        }
    }

    vector<int> positionstoerase;
    for (unsigned int i=0; i<quads_ut.size(); i++)
        if (quads_ut[i].size() < 4)
            positionstoerase.push_back(i);

    for (unsigned int i=0; i<positionstoerase.size(); i++)
        quads_ut.erase(quads_ut.begin() + positionstoerase[i]);
}

void FEM::SetSecondLayer_t() {
	for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
	{
		vector<float> point;
		float pos1 = vMPsXYZN_t[i][0] - h;
		float pos2 = vMPsXYZN_t[i][1] - h;
		float pos3 = vMPsXYZN_t[i][2] - h;
		point.push_back(pos1);
		point.push_back(pos2);
		point.push_back(pos3);
		vMPsXYZN_t2.push_back(point);
	}

	// Display bottom layer
	vector<vector<float> > allpoints;

	for (unsigned int i=0; i<vMPsXYZN_t.size(); i++) {
        allpoints.push_back(vMPsXYZN_t[i]);
    }

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++) {
        allpoints.push_back(vMPsXYZN_t2[i]);
    }

	sCloud_quads12.width = allpoints.size();
	sCloud_quads12.height = 1;
	sCloud_quads12.is_dense = false;
	sCloud_quads12.points.resize(sCloud_quads12.width*sCloud_quads12.height);

	for (unsigned int i=0; i<allpoints.size(); i++)
	{
        sCloud_quads12.points[i].x = allpoints[i][0];
        sCloud_quads12.points[i].y = allpoints[i][1];
        sCloud_quads12.points[i].z = allpoints[i][2];
	}

    // Save both layers
    pcPath_t12 += std::to_string(nFrameId);
    pcPath_t12 += ".pcd";

    if (sCloud_quads12.width==0)
        cout << "No data to save" << endl;
    else
        pcl::io::savePCDFileASCII (pcPath_t12, sCloud_quads12);

}





/// 2 - FEM CALCULATIONS

/*
void FEM::PrintBehaviorMatrix()
{
    cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
    cout << "lambda = " << lambda << endl;
    cout << "G = " << G << endl;
    cout << "D = " << endl;
    for (unsigned int i=0; i<6; i++)
    {
        for (unsigned int j=0; j<6; j++)
        {
            cout << D[i][j] << "\t";
        }
        cout << endl;
    }
    cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
}



void FEM::SetGaussLegendre (float input_fg1, float input_fg2)
{
    fg1 = input_fg1;
    fg2 = input_fg2;
    fg = input_fg1;

    vector<float> g0;     g0.push_back(+fg);  g0.push_back(+fg);  g0.push_back(+fg);     gs.push_back(g0);
    vector<float> g1;     g1.push_back(+fg);  g1.push_back(+fg);  g1.push_back(-fg);     gs.push_back(g1);
    vector<float> g2;     g2.push_back(+fg);  g2.push_back(-fg);  g2.push_back(+fg);     gs.push_back(g2);
    vector<float> g3;     g3.push_back(+fg);  g3.push_back(-fg);  g3.push_back(-fg);     gs.push_back(g3);
    vector<float> g4;     g4.push_back(-fg);  g4.push_back(+fg);  g4.push_back(+fg);     gs.push_back(g4);
    vector<float> g5;     g5.push_back(-fg);  g5.push_back(+fg);  g5.push_back(-fg);     gs.push_back(g5);
    vector<float> g6;     g6.push_back(-fg);  g6.push_back(-fg);  g6.push_back(+fg);     gs.push_back(g6);
    vector<float> g7;     g7.push_back(-fg);  g7.push_back(-fg);  g7.push_back(-fg);     gs.push_back(g7);
}



void FEM::SetElementDepth (float input_h)
{
    h = input_h;
}
*/




bool FEM::MatrixAssembly() {
    cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl
         << "- - - - - - - - - - - - - MatAssembly - - - - - - - - - - - - -" << endl
         << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl
         << "- DOFe = 24 degrees of freedom per element" << endl
         << "- nEl = " << quads_t.size() << " elements (quads)" << endl
         << "- Ksize' = " << quads_t.size()*8*3 << " total number of DoF in the domain" << endl;  // DoF = nEl*nodesperelement*DOFnode

    // To easily insert the nodes into the stiffness matrix, they are not simply put in
    // one after the other. Instead, we use the node indexes to set their position in
    // the matrix. This way we will get a slightly bigger stiffness matrix, but this fact
    // won't affect the calculations, since NULL positions will be set to 0 by default.

    // First we set the size of the stiffness matrix, which will be equal to the size of
    // the MPs vector, the positions corresponding to a offvertex will be set to zero in
    // the stiffness matrix
    int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    Ksize = 3*nTotalNodes;
    //cout << "  nTotalNodes/Ksize = " << nTotalNodes << "/" << Ksize << endl;

    // Make sure we have points before continue.
    if (Ksize<=3)
        return false;

    int nSize1stLayer = vMPsXYZN_t.size();

    K.clear();
    vector<float> ki = vector<float>(Ksize,0.0);
    for (unsigned int i=0; i<Ksize; i++)
        K.push_back(ki);

    //cout << " - Assembly - K[" << K.size() << "][" << K[0].size() << "] initialized to 0.0" << endl;

    for (unsigned int i=0; i<quads_t.size(); i++)
    {
        vector<int> nodes;
        nodes.push_back(quads_t[i][0]);
        nodes.push_back(quads_t[i][1]);
        nodes.push_back(quads_t[i][2]);
        nodes.push_back(quads_t[i][3]);
        nodes.push_back(quads_t[i][0] + nSize1stLayer);
        nodes.push_back(quads_t[i][1] + nSize1stLayer);
        nodes.push_back(quads_t[i][2] + nSize1stLayer);
        nodes.push_back(quads_t[i][3] + nSize1stLayer);

        vector<vector<float> > vfPts;
        for (unsigned int j=0; j<nodes.size(); j++)
            vfPts.push_back(vMPsXYZN_t[j]); //modded

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

void FEM::MatrixAssembly_ut() {
    // Ke
    // Assembly applying connections

    //cout << " - Assembly - Begin " << endl;

    //unsigned int n = 24;  // number of degrees of freedom per element
    //unsigned int nEl = quads_ut.size();	// number of elements
    //unsigned int Ksize = nEl * 8 * 3;	// total number of degrees of freedom in the domain. numberofelements * nodesperelement * DoFpernode

    //int nQuads = quads_ut.size() * 4;

    int maxId = 0;
    for (unsigned int i=0; i<quads_ut.size(); i++)
        for (unsigned int j=0; j<quads_ut[i].size(); j++)
            if (quads_ut[i][j]>maxId)
                maxId = quads_ut[i][j];

    int nMaxId1stLayer = maxId;
    int nMaxId2ndLayer = 2*nMaxId1stLayer;

    Ksize = 3*nMaxId2ndLayer + 3;

    K.clear();
    vector<float> ki = vector<float>(Ksize,0.0);
    for (unsigned int i=0; i<Ksize; i++)
        K.push_back(ki);

    for (unsigned int i=0; i<quads_ut.size(); i++)
    {
        vector<int> nodes;
        nodes.push_back(quads_ut[i][0]);
        nodes.push_back(quads_ut[i][1]);
        nodes.push_back(quads_ut[i][2]);
        nodes.push_back(quads_ut[i][3]);
        nodes.push_back(quads_ut[i][0] + nMaxId1stLayer);
        nodes.push_back(quads_ut[i][1] + nMaxId1stLayer);
        nodes.push_back(quads_ut[i][2] + nMaxId1stLayer);
        nodes.push_back(quads_ut[i][3] + nMaxId1stLayer);

        cout << " Quad " << i << "/" << quads_ut.size()-1 << ": ";
        for (unsigned int j=0; j<nodes.size(); j++)
            cout << nodes[j] << ".";
        cout << endl;

        vector<int> matpos;
        int matpos_max = 0;
        for (unsigned j=0; j<nodes.size(); j++)
        {
            matpos.push_back(nodes[j]*3);
            if (nodes[j]*3 > matpos_max)
                matpos_max = nodes[j]*3;
        }

        for (unsigned int node1=0; node1<nodes.size(); node1++)
            for (unsigned int node2=0; node2<nodes.size(); node2++)
                for (unsigned int m=0; m<3; m++)
                    for (unsigned int n=0; n<3; n++)
                        K[matpos[node1]+m][matpos[node2]+n] += Ke[node1*3+m][node2*3+n];
    }
}


/*
bool FEM::ComputeK1() {
    // Método Gauss Jordan para cálculo de matriz inversa
    // (A|I) -> (I|A1)

    vector<vector<float> > Ktemp = K;

    K1.clear();
    K1 = vector<vector<float> >(Ksize,vector<float>(Ksize,(0.0)));
    for (unsigned int i=0; i<K1.size(); i++)
        K1[i][i] = 1.0;

    for (unsigned int i=0; i<Ksize; i++)
    {
        float fPivot = Ktemp[i][i];

        if (fPivot==0)
        {
            unsigned int nRowObj = i;
            if (i<(Ksize-1))
                for (unsigned int j=i+1; j<Ksize; j++)
                {
                    if (Ktemp[j][i]!=0)
                    {
                        nRowObj = j;
                        break;
                    }
                }

            //cout << "Ksize = " << Ksize << endl;
            //cout << "nRowObj = " << nRowObj << endl;

            if (nRowObj==i)
                return false;

            vector<int> vSwitch;
            vSwitch.push_back(i);
            vSwitch.push_back(nRowObj);
            vnSwitchedRows.push_back(vSwitch);

            vector<float> vfAux = Ktemp[i];

            for (unsigned int aa=0; aa<Ksize; aa++)
            {
                Ktemp[i][aa] = Ktemp[nRowObj][aa];
                Ktemp[nRowObj][aa] = vfAux[aa];
            }
        }

        fPivot = Ktemp[i][i];

        if (fPivot!=1)
            for (unsigned int j=0; j<Ksize; j++)
            {
                Ktemp[i][j] /= fPivot;
                K1[i][j] /= fPivot;
            }


        for (unsigned int j=0; j<Ksize; j++)
        {
            if (j==i)
                continue;
            float fPivotRow = Ktemp[j][i];
            for (unsigned int s=0; s<Ksize; s++)
            {
                Ktemp[j][s] -= fPivotRow*Ktemp[i][s];
                K1[j][s] -= fPivotRow*K1[i][s];
            }
        }
    }

    return true;
}

*/

void FEM::Set_u0(vector<MapPoint*> vpMPs) {
    u0.clear();

    //int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
        for (unsigned int j=0; j<3; j++)
            u0.push_back(vMPsXYZN_t[i][j]);

    for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
        for (unsigned int j=0; j<3; j++)
            u0.push_back(vMPsXYZN_t2[i][j]);

    //cout << "nTotalNodes/Ksize/u0size = " << nTotalNodes << "/" << Ksize << "/" << u0.size() << endl;
}

void FEM::Set_uf(vector<MapPoint*> vpMPs) {
    uf.clear();

    vector<vector<float> > top_layer;
    top_layer.resize(vMPsXYZN_t.size());

    for (unsigned int i=0; i<vpMPs.size(); i++)
    {
        MapPoint* pMP = vpMPs[i];

        vector<float> vfMP_top;

        cv::Mat pos = pMP->GetWorldPos();

        cv::Vec3f vPos = pos;

        float fPos1 = vPos.val[0];
        float fPos2 = vPos.val[1];
        float fPos3 = vPos.val[2];

        vfMP_top.push_back(fPos1);
        vfMP_top.push_back(fPos2);
        vfMP_top.push_back(fPos3);

        top_layer[i] = vfMP_top;
    }

    for (unsigned int i=0; i<triangles_t.size(); i++)
    {
        // Get vertices and indices for refferece triangle
        int v0idx = triangles_t[i][0];
        vector<float> v0 = top_layer[v0idx];

        int v1idx = triangles_t[i][1];
        vector<float> v1 = top_layer[v1idx];

        int v2idx = triangles_t[i][2];
        vector<float> v2 = top_layer[v2idx];

        // Set variables to store the ressult
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

    int nTotalNodes1 = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    cout << "             nTotalNodes/Ksize/u0size = " << nTotalNodes1 << "/" << Ksize << "/" << u0.size() << endl;

    int nTotalNodes2 = top_layer.size() + bottom_layer.size();
    int ksize2 = 3*nTotalNodes2;
    cout << "             nTotalNodes/Ksize/ufsize = " << nTotalNodes2 << "/" << ksize2 << "/" << uf.size() << endl;
}

void FEM::ComputeDisplacement() {
    a.clear();
    for (unsigned int i=0; i<u0.size(); i++)
        a.push_back(uf[i]-u0[i]);
}

void FEM::ComputeForces() {
    f.clear();
    f.resize(Ksize);
    for (unsigned int i=0; i<Ksize; i++)
    {
        float fi = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            fi += K[i][j] * a[j];
        f[i] = fi;
    }
}

void FEM::SetForces() {
    int nForcesKnown = f.size();

    for (unsigned int i=nForcesKnown; i<Ksize; i++)
        f.push_back(0.0);
}

void FEM::ComputeUnknownDisplacements() {
    // a = K1 · f
    a.clear();
    a.resize(Ksize);

    for (unsigned int i=0; i<Ksize; i++)
    {
        float ai = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            ai += K1[i][j] * f[j];
        a[i] = ai;
    }
}

void FEM::ComputeStrainEnergy() {
    // sE = a' · K · a
    vector<float> aK = vector<float>(Ksize,0.0);
    for (unsigned int i=0; i<Ksize; i++)
    {
        float prod = 0.0;
        for (unsigned int j=0; j<Ksize; j++)
            prod += a[j] * K[j][i];
        aK[i] = prod;
    }

    sE = 0.0;
    for (unsigned int i=0; i<Ksize; i++)
        sE += aK[i] * a[i];
}



float FEM::GetStrainEnergy() { return sE; }

float FEM::GetH() { return h; }

/*vector<vector<float> > FEM::GetKe() { return K; }*/

vector<vector<float> > FEM::GetK() { return K; }

vector<vector<int> > FEM::GetTrianglesT() { return triangles_t; }

vector<vector<int> > FEM::GetQuadsT() { return quads_t; }

vector<float> FEM::GetU0() { return u0; }

vector<vector<float> > FEM::GetPointsL1() { return vMPsXYZN_t; }

vector<vector<float> > FEM::GetPointsL2() { return vMPsXYZN_t2; }

vector<vector<int> > FEM::GetIdxNewVtx() { return vIdxNewVertices; }

vector<vector<int> > FEM::GetNewPointsBase() { return vNewPointsBase; }

vector<long unsigned int> FEM::GetMPsGlobalID() { return vMPsGlobalID_t; }

int FEM::GetKsize() { return Ksize; }









/*
void FEM::SaveQuadMesh_t()
{
    /// Record mesh

    pcl::PointCloud<pcl::PointXYZ> pc_2save;

    // Load top layer into cloud
    pc_2save.width = vMPsXYZN_t.size()/2; //cloudTracked.width    = 5;
    pc_2save.height   = 1;
    pc_2save.is_dense = false;
    pc_2save.points.resize (pc_2save.width * pc_2save.height);

    unsigned int sizetoplayer = pc_2save.points.size();

    for (size_t i = 0; i < sizetoplayer; i++)
    {
        pc_2save.points[i].x = vMPsXYZN_t[i][0];
        pc_2save.points[i].y = vMPsXYZN_t[i][1];
        pc_2save.points[i].z = vMPsXYZN_t[i][2];
    }

    string pcPath_2save = "../../../output/PointClouds/pc_2save";
    pcPath_2save += std::to_string(nFrameId);
    pcPath_2save += ".pcd";

    if (pc_2save.width>0)
        pcl::io::savePCDFileASCII (pcPath_2save, pc_2save);
    else
        cout << "No data to save" << endl;


    /// Record connections

    string pcPath_2savetxt = "../../../output/PointClouds/pc_2save";
    pcPath_2savetxt += std::to_string(nFrameId);
    pcPath_2savetxt += ".txt";


    //int maxsize = vMPsXYZN_t.size();
    vector<vector<int> > cn2save;

    for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
    {
        vector<int> cni;
        cni.push_back(i);
        cn2save.push_back(cni);
    }

    for (unsigned int i=0; i<quads_t.size(); i++)
    {
        int idx0 = quads_t[i][0];
        int idx1 = quads_t[i][1];
        int idx2 = quads_t[i][2];
        int idx3 = quads_t[i][3];

        if ( !( std::find(cn2save[idx0].begin(), cn2save[idx0].end(), idx1) != cn2save[idx0].end() ) )
            cn2save[idx0].push_back(idx1);
        if ( !( std::find(cn2save[idx0].begin(), cn2save[idx0].end(), idx3) != cn2save[idx0].end() ) )
            cn2save[idx0].push_back(idx3);

        if ( !( std::find(cn2save[idx1].begin(), cn2save[idx1].end(), idx0) != cn2save[idx1].end() ) )
            cn2save[idx1].push_back(idx0);
        if ( !( std::find(cn2save[idx1].begin(), cn2save[idx1].end(), idx2) != cn2save[idx1].end() ) )
            cn2save[idx1].push_back(idx3);

        if ( !( std::find(cn2save[idx2].begin(), cn2save[idx2].end(), idx1) != cn2save[idx2].end() ) )
            cn2save[idx2].push_back(idx1);
        if ( !( std::find(cn2save[idx2].begin(), cn2save[idx2].end(), idx3) != cn2save[idx2].end() ) )
            cn2save[idx2].push_back(idx3);

        if ( !( std::find(cn2save[idx3].begin(), cn2save[idx3].end(), idx0) != cn2save[idx3].end() ) )
            cn2save[idx3].push_back(idx0);
        if ( !( std::find(cn2save[idx3].begin(), cn2save[idx3].end(), idx2) != cn2save[idx3].end() ) )
            cn2save[idx3].push_back(idx3);
    }


    ofstream cntxt;
    cntxt.open(pcPath_2savetxt);

    for (unsigned int i=0; i<cn2save.size(); i++)
    {
        unsigned int idxs = 0;
        for (unsigned int j=0; j<cn2save[i].size(); j++)
        {
            cntxt << cn2save[i][j];
            idxs++;
            if (idxs<(cn2save[i].size()-1))
                cntxt << " ";
        }
        cntxt << endl;
    }

    cntxt.close();
}
*/


} //namespace ORB_SLAM2
