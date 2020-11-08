/**
* This file is an addon to ORB-SLAM2.
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

#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <thread>

#include "MapPoint.h"
#include "Thirdparty/g2o/g2o/FEA/include/FEA.h"

// OPENCV LIBRARIES
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

// PCL LIBRARIES
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/point_types.h>

#include <pcl/common/common.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/mls.h>
#include <pcl/surface/impl/mls.hpp>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/cloud_viewer.h>

#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/surface/gp3.h>

// OTHER LIBRARIES
#include <Eigen/StdVector>
#include <algorithm>
#include <pcl/visualization/vtk.h>

using namespace std;
using namespace pcl;

namespace ORB_SLAM2
{
class FEM
{
public:	// FUNCTIONS

	// Constructor & Destructor
	FEM(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug);
    ~FEM();

    // Build mesh and compute K
    bool Compute(int nMode);

    // Get coordinates and normals from MapPoints
    void GetMapPointCoordinates(int nMode);

    // Load MapPoints into cloud
    bool LoadMPsIntoCloud(int nMode);

    // Moving Least Squares
    bool MLS(int nMode);

    // Triangulation, greedy projection (meshing)
    bool ComputeMesh(int nMode);

    void CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc, int *pnInMu, int *pnSearchRad, int *pnMaxNeig, int *pnSurAng, int *pnMinAng, int *pnMaxAng);

    // Turns a triangle mesh into a quadrilateral mesh, 1 tri = 3 quads
    // The function adds a node in the middle of each side of the original triangle, and an extra node in the ortocenter.
    void tri2quad(int nMode);
    void tri2quad_t();
    void tri2quad_u();

    // Duplicates the previously generated mesh at a distance, builds a 3D structure.
    void SetSecondLayer(int nMode);

    void ReadyFEA(FEA* pFEA, int nMode);

    vector<vector<float> > ComputeKei(vector<vector<float> > vfPts);

    bool MatrixAssembly(int nMode);

    void ImposeDirichletEncastre(vector<vector<int> > vD, float Klarge);

    vector<vector<float> > InvertMatrixEigen(vector<vector<float> > m1);
    vector<vector<float> > MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2);

    void Set_u0(vector<MapPoint*> vpMPs, int nMode);

    void Set_uf(vector<MapPoint*> vpMPs);

    void ComputeDisplacement();

    void ComputeForces(); // K·a = f

    // Access functions
    float GetStrainEnergy();
    vector<vector<float> > GetK();
    vector<vector<int> > GetTrianglesT();
    vector<vector<int> > GetQuadsT();
    vector<float> GetU0();
    vector<vector<float> > GetPointsL1();
    vector<vector<float> > GetPointsL2();
    vector<vector<int> > GetIdxNewVtx();
    vector<vector<int> > GetNewPointsBase();
    int GetKsize();
    vector<long unsigned int> GetMPsGlobalID();

public:	// VARIABLES

    vector<MapPoint*> vpMPs_t;
    vector<MapPoint*> vpMPs_ut;
    vector<cv::KeyPoint*> vpKPs_t;

    // IndexTrackedMpInFEM - IndexTrackedMpInFrame
    vector<int> idxMpF;

    vector<vector<float> > vMPsXYZN_t;
    vector<vector<float> > vMPsXYZN_t2;

    vector<vector<float> > vMPsXYZN_ut;
    vector<vector<float> > vMPsXYZN_ut2;

    vector<vector<MapPoint*> > vpMPs2Draw;
    vector<vector<MapPoint*> > vpMPs2Drawu;
    vector<vector<cv::KeyPoint*> > vpKPs2Draw;
    vector<vector<cv::KeyPoint*> > vpKPs2Drawu;

    bool bEInverse = false;


private:

    // Set MapPoints as active or not
    vector<bool> vbtMPsActive;
    vector<bool> vbuMPsActive;

    // Staring PointCloud
    pcl::PointCloud<pcl::PointXYZ> pc_t_0;
    pcl::PointCloud<pcl::PointXYZ> pc_u_0;

    // Smoothed PointCloud (MLS)
    pcl::PointCloud<pcl::PointNormal> pc_t_1;
    pcl::PointCloud<pcl::PointNormal> pc_u_1;
    vector<int> vtindices;
    vector<int> vuindices;

    // Reconstructed mesh
    PolygonMesh mesh_t;
    PolygonMesh mesh_u;

    vector<vector<int> > triangles_t;
    vector<vector<int> > triangles_u;
    vector<vector<vector<int> > > edgespertriangle_t;
    vector<vector<vector<int> > > edgespertriangle_u;

    vector<vector<int> > vIdxNewVertices;
    vector<vector<int> > vIdxNewVertices_u;
    vector<vector<int> > vNewPointsBase;
    vector<vector<int> > vNewPointsBase_u;

    vector<vector<int> > quads_t;
    vector<vector<int> > quads_u;

    // Frame
    unsigned int nFrameId;

    // Young Modulus [Pa] & Poisson Coefficient
    unsigned int E;
    float nu;

    // Lamé parameters & Behaviour matrix
    float lambda = 0.0;
    float G = 0.0;
    vector<vector<float> > D;

    // Gauss Points
    float fg;
    vector<vector<float> > gs;

    // Element depth
    float h;

    // Elemental matrix
    vector<vector<float> > Ke;

    // Stiffness matrix
    unsigned int Ksize = 0;
    vector<vector<float> > K;
    float DetK = 0.0;

    unsigned int Kusize = 0;
    vector<vector<float> > Ku;
    float DetKu = 0.0;

    vector<vector<float> > K1;
    float DetK1 = 0.0;
    vector<vector<int> > vnSwitchedRows;

    // Displacement
    vector<float> u0;  // Starting position of mappoints
    vector<float> u0u;  // Starting position of mappoints
    vector<float> uf;  // New position after being modified by the optimization
    vector<float> u;   // u = uf-u0
    vector<float> a = vector<float>(Ksize,0.0);     // VectorAssembly(u)

    // Forces
    vector<float> f = vector<float>(Ksize,0.0);

    // Strain energy
    float sE = 0.0;

    // Debug mode
    bool bDebugMode = false;

}; // class FEM
} // namespace ORB_SLAM
#endif // FEM_H