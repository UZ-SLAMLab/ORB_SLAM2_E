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

#include "MapPoint.h"
#include "Thirdparty/g2o/g2o/FEA/include/FEA.h"

// OPENCV LIBRARIES
#include <opencv2/opencv.hpp>

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

	//Constructor
	FEM(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug);
    //FEM(unsigned int input_nFrameId, unsigned int input_E, float input_nu, bool bSetDebug);

    //Destructor
    ~FEM();

    bool Compute(int nMode);

    // Get coordinates and normals from MapPoints
    void GetMapPointCoordinates(int nMode);

    // Load MapPoints into cloud
    bool LoadMPsIntoCloud();
    bool LoadMPsIntoCloud(int nMode);
    bool LoadTrackedMPsIntoCloud();
    bool LoadUnTrackedMPsIntoCloud();


    // Moving Least Squares
    bool MLS();
    bool MLS(int nMode);
    bool MLS_t();
    bool MLS_ut();

    // Triangulation, greedy projection (meshing)
    bool ComputeMesh();
    bool ComputeMesh(int nMode);
    bool ComputeMesh_t();

    void CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc, int *pnInMu, int *pnSearchRad, int *pnMaxNeig, int *pnSurAng, int *pnMinAng, int *pnMaxAng);

    // Get mesh connections
    void GetMeshConnections();
    void GetMeshConnections_t();
    void GetMeshConnections_ut();

    /// tri2quad
    // Turns a triangle mesh into a quadrilateral mesh
    // The function adds a node in the middle of each side of the original triangle, and an extra node in the ortocenter.
    // Therefore each triangle is converted into 3 quadrilaterals.
    /*void tri2quad();*/
    void tri2quad(vector<vector<int> > triangles, vector<bool> vbMPsActive, vector<int> vinindices, vector<vector<int> > vIdxNewVertices, vector<vector<int> > quads, vector<vector<int> > vNewPointsBase, vector<vector<float> > vMPsXYZN);
    void tri2quad(int nMode);
    void tri2quad_t();
    void tri2quad_u();
    void tri2quad_ut();


    /// SaveQuadMesh
    // Saves the quad mesh into a vtk file for offline use.
    void SaveQuadMesh_t();


    /// SetSecondLayer
    // Duplicates the previously generated mesh at a certain distance, generating a second layer similar to the original
    // one but placed at a distance, allowing for the construction of a 3D structure.
    //void SetSecondLayer();
    // Generates a second layer from both matched and unmatched map points.
    void SetSecondLayer_t();
    // Generates a second layer from the matched points.
    //void SetSecondLayer_ut();
    // Generates a second layer from the unmatched points.
    void SetSecondLayer(int nMode);




    void ReadyFEA(FEA* pFEA, int nMode);


    void PrintBehaviorMatrix();

    void SetGaussLegendre (float input_fg1, float input_fg2);
    void SetElementDepth (float input_h);

    vector<vector<float> > GetKe();

    //float AproxFunction(int i, int j, float y1, float y2, float y3, float z2, float z3, float h, float Xi, float Eta, float Zeta);
    vector<vector<float> > ComputeKei(vector<vector<float> > vfPts);


    bool MatrixAssembly();
    bool MatrixAssembly(int nMode);
    void MatrixAssembly_ut();

    float ComputeDet(vector<vector<float> > Kin);

    bool ComputeK1();

    void Set_u0(vector<MapPoint*> vpMPs);
    void Set_u0(vector<MapPoint*> vpMPs, int nMode);

    void Set_uf(vector<MapPoint*> vpMPs);

    void ComputeDisplacement();

    void ComputeForces(); // K·a = f

    void SetForces();

    void ComputeUnknownDisplacements();

    void ComputeStrainEnergy();

    // Access functions
    float GetStrainEnergy();
    float GetH();
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

	// Tracked MapPoints
    vector<MapPoint*> vpMPs_t;
    vector<MapPoint*> vpMPs_ut;
    vector<cv::KeyPoint*> vpKPs_t;

    // IndexTrackedMpInFEM - IndexTrackedMpInFrame
    vector<int> idxMpF;

    //MapPoint coordinates and normals
    vector<vector<float> > vMPsXYZN_t;
    vector<vector<float> > vMPsXYZN_t2;

    vector<long unsigned int> vMPsGlobalID_t;
    vector<long unsigned int> vMPsGlobalID_u;

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
    pcl::PointCloud<pcl::PointNormal> cloudUnTracked;
    string pcPath_t = "output/PointClouds/pc_t_f";
    string pcPath_ut = "output/PointClouds/pc_ut_f";

    // Smoothed PointCloud (MLS)
    pcl::PointCloud<pcl::PointNormal> pc_t_1;
    pcl::PointCloud<pcl::PointNormal> pc_u_1;
    vector<int> vtindices;
    vector<int> vuindices;
    pcl::PointCloud<pcl::PointNormal> sCloud_ut;
    string pcsPath_t = "output/PointClouds/pcs_t_f";
    string pcsPath_ut = "output/PointClouds/pcs_ut_f";

    // Reconstructed mesh
    PolygonMesh mesh_t;
    PolygonMesh mesh_u;
    PolygonMesh mesh_ut;
    pcl::PointCloud<pcl::PointXYZ>::Ptr reconstructed_cloud_t;// (new pcl::PointCloud<pcl::PointXYZ>);
    string pcrPath_t = "output/PointClouds/pcr_t_f";
    string pcrPath_u = "output/PointClouds/pcr_t_f";
    string pcrPath_ut = "output/PointClouds/pcr_ut_f";
    int nVerticesInMesh_t = 0;
    int nVerticesInMesh_ut = 0;

    // Final Point cloud with both layers
    pcl::PointCloud<pcl::PointXYZ> sCloud_quads12;
    string pcPath_t12 = "output/PointClouds/pc2_t_f";

    vector<vector<int> > triangles_t;
    vector<vector<int> > triangles_u;
    vector<vector<int> > triangles_ut;
    vector<vector<int> > edges_t;
    vector<vector<int> > edges_ut;
    vector<vector<vector<int> > > edgespertriangle;
    vector<vector<vector<int> > > edgespertriangle_t;
    vector<vector<vector<int> > > edgespertriangle_u;
    vector<vector<vector<int> > > edgespertriangle_ut;
    vector<int> nNotT2Q_t;
    vector<int> nNotT2Q_ut;

    vector<vector<int> > vIdxNewVertices;
    vector<vector<int> > vIdxNewVertices_t;
    vector<vector<int> > vIdxNewVertices_u;
    vector<vector<int> > vNewPointsBase;
    vector<vector<int> > vNewPointsBase_t;
    vector<vector<int> > vNewPointsBase_u;

    vector<vector<int> > quads_t;
    vector<vector<int> > quads_u;
    vector<vector<int> > quads_ut;

    vector<vector<int> > edgesq_t;
    vector<vector<int> > edgesq_ut;

    // Frame
    unsigned int nFrameId;

    // Young Modulus [Pa]
    unsigned int E;

    // Poisson Coefficient
    float nu;

    // Behavior matrix
    float D[6][6]= {};

    // Lamé parameters
    float lambda = 0.0;
    float G = 0.0;

    // Gauss Points
    float fg1 = 0.0;
    float fg2 = 0.0;
    float fg = 0.0;
    vector<vector<float> > gs;

    // Element depth
    float h = 0.0;

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
