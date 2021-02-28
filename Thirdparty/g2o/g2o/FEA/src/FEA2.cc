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

#include "../include/FEA2.h"


FEA2::FEA2(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, int input_nElType, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg(input_fg1), nElType(input_nElType), bDebugMode(bSetDebug){
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


FEA2::~FEA2(){
}


bool FEA2::Compute(int nMode) {

    GetMapPointCoordinates(nMode);

    if (LoadMPsIntoCloud(nMode)) {
        if (bDebugMode) cout << "           - MLS" << endl;
        if (MLS(nMode)) {
            if (bDebugMode) cout << "           - Tri mesh" << endl;
            if(ComputeMesh(nMode)) {
                if (nElType==1)         //C3D6
                    tri2quad(nMode);

                SetSecondLayer(nMode);

                if (nMode == 1)
                    Set_u0(vpMPs_t,nMode);
                else if (nMode == 2)
                    Set_u0(vpMPs_ut,nMode);

                if (nElType==1)
                    MatrixAssemblyC3D8(nMode);
                if (nElType==2)
                    MatrixAssemblyC3D6(nMode);

                if (nMode == 1)
                    ImposeDirichletEncastre_K(nMode,vvDir_t,100000000.0);
                else if (nMode == 2)
                    ImposeDirichletEncastre_K(nMode,vvDir_u,100000000.0);


                if (nMode == 2){
                    K1 = InvertMatrixEigen(K);
                    UpdateForces();
                }
                    

                return true;
            }
        }
    }
    return false;
}


void FEA2::GetMapPointCoordinates(int nMode) {
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
        }
    }
}


bool FEA2::LoadMPsIntoCloud(int nMode) {
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

        if (pc_t_0.width>0)
            return true;
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

        if (pc_u_0.width>0)
            return true;
        else
            return false;
    }

    else{
        return false;
    }
}


bool FEA2::MLS(int nMode) {
    if (nMode == 1 ) {
        PointCloud<PointXYZ>::Ptr pc_t_mls_0 (new pcl::PointCloud<PointXYZ> (pc_t_0));

        // Create KD-Tree & MLS & set parameters
        pcl::search::KdTree<pcl::PointXYZ>::Ptr tree1 (new pcl::search::KdTree<pcl::PointXYZ>);
        pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls1;
        mls1.setComputeNormals (true);
        mls1.setInputCloud (pc_t_mls_0);
        //mls1.setPolynomialFit (true);
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
        if (pc_t_1.width>0)
            return true;
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
        //mls1.setPolynomialFit (true);
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
        if (pc_u_1.width>0)
            return true;
        else
            return false;
    }
    else {
        return false;
    }
}


bool FEA2::ComputeMesh(int nMode) {
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
        //CalculateGP3Parameters(ptr_pc_t_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);
        //cout << "InMu = " << nInMu << endl;
        //cout << "nSearchRad = " << nSearchRad << endl;
        //cout << "nMaxNeig = " << nMaxNeig << endl;
        //cout << "nSurAng = " << nSurAng << endl;
        //cout << "nMinAng = " << nMinAng << endl;
        //cout << "nMaxAng = " << nMaxAng << endl;

        //nSearchRad = 30;
        //nSurAng = 120;
        //nMaxNeig = 1000000;

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
        //CalculateGP3Parameters(ptr_pc_u_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);

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
            }

            return true;
        }
        else {
            cout << "No data" << endl;
            return false;
        }
    }
    else{
        return false;
    }

}


void FEA2::CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc, int *pnInMu, int *pnSearchRad, int *pnMaxNeig, int *pnSurAng, int *pnMinAng, int *pnMaxAng) {
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


void FEA2::tri2quad(int nMode) {
        if (nMode == 1)
            tri2quad_t();
        else if (nMode == 2)
            tri2quad_u();
}


void FEA2::tri2quad_t() {
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


void FEA2::tri2quad_u() {
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


void FEA2::SetSecondLayer(int nMode) {
    if (nMode == 1) {
        vvDir_t.clear();
        for (unsigned int i=0; i<vMPsXYZN_t.size(); i++) {
            vector<float> point;
            float pos1 = vMPsXYZN_t[i][0] - h;
            float pos2 = vMPsXYZN_t[i][1] - h;
            float pos3 = vMPsXYZN_t[i][2] - h;
            point.push_back(pos1);
            point.push_back(pos2);
            point.push_back(pos3);
            vMPsXYZN_t2.push_back(point);

            vector<int> vDiri;
            vDiri.push_back(vMPsXYZN_t.size() + i);
            vvDir_t.push_back(vDiri);
        }
    }
    else if (nMode == 2) {
        vvDir_u.clear();
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++) {
            vector<float> point;
            float pos1 = vMPsXYZN_ut[i][0] - h;
            float pos2 = vMPsXYZN_ut[i][1] - h;
            float pos3 = vMPsXYZN_ut[i][2] - h;
            point.push_back(pos1);
            point.push_back(pos2);
            point.push_back(pos3);
            vMPsXYZN_ut2.push_back(point);

            vector<int> vDiri;
            vDiri.push_back(vMPsXYZN_ut.size() + i);
            vvDir_u.push_back(vDiri);
        }
    }
}


void FEA2::Set_u0(vector<MapPoint*> vpMPs, int nMode) {
    /*
    cout << "Set_u0" << endl;
    cout << "nMode=" << nMode << endl;
    if (nMode == 1 ) {
        cout << "1" << endl;
        u0.clear();
        vector<vector<float> > vt1 = vector_resize_cols(vMPsXYZN_t,1);
        vector<vector<float> > vt2 = vector_resize_cols(vMPsXYZN_t2,1);
        u0.resize(vt1.size()+vt2.size());
        cout << "u0=" << u0.size() << "  vt1=" << vt1.size() << "  vt2=" << vt2.size() << endl;
        for (unsigned int i=0; i<vt1.size(); i++)
            u0[i] = vt1[i][0];
        for (unsigned int i=0; i<vt2.size(); i++)
            u0[i] = vt2[i][0];
    }
    else if (nMode == 2) {
        cout << "2" << endl;
        u0u.clear();
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_ut[i]; j++)
                u0u.push_back(vMPsXYZN_ut[i][j]);
        for (unsigned int i=0; i<vMPsXYZN_ut2.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_ut2[i]; j++)
                u0u.push_back(vMPsXYZN_ut2[i][j]);
    }
    */



    if (nMode == 1) {
        u0.clear();
        for (unsigned int i=0; i<vMPsXYZN_t.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_t[i].size(); j++)
                u0.push_back(vMPsXYZN_t[i][j]);
        for (unsigned int i=0; i<vMPsXYZN_t2.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_t2[i].size(); j++)
                u0.push_back(vMPsXYZN_t2[i][j]);
    }
    else if (nMode == 2) {
        u0u.clear();
        for (unsigned int i=0; i<vMPsXYZN_ut.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_ut[i].size(); j++)
                u0u.push_back(vMPsXYZN_ut[i][j]);
        for (unsigned int i=0; i<vMPsXYZN_ut2.size(); i++)
            for (unsigned int j=0; j<vMPsXYZN_ut2[i].size(); j++)
                u0u.push_back(vMPsXYZN_ut2[i][j]);
    }
}


vector<vector<float> > FEA2::ComputeKeiC3D8(vector<vector<float> > vfPts) {

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


vector<vector<float> > FEA2::ComputeKeiC3D6(vector<vector<float> > vfPts) {

    vector<vector<float> > vBtDB = vector<vector<float> >(18,vector<float>(18,0.0));

    for (unsigned int ops=0; ops<gs.size(); ops++)
    {
        float xi = gs[ops][0];
        float eta = gs[ops][1];
        float zeta = gs[ops][2];

        float dN1dxi = -(1 + zeta)/2;    float dN1deta = -(1 + zeta)/2;    float dN1dzeta =  (1-xi-eta)/2;
        float dN2dxi =  (1 + zeta)/2;    float dN2deta =   0.0;            float dN2dzeta =  xi/2;
        float dN3dxi =  0.0;             float dN3deta =  (1 + zeta)/2;    float dN3dzeta =  eta/2;
        float dN4dxi = -(1 - zeta)/2;    float dN4deta = -(1 - zeta)/2;    float dN4dzeta = -(1-xi-eta)/2;
        float dN5dxi =  (1 - zeta)/2;    float dN5deta =   0.0;            float dN5dzeta = -xi/2;
        float dN6dxi =  0.0;             float dN6deta =  (1 - zeta)/2;    float dN6dzeta = -eta/2;

        /*dxdxi*/   float J_00 = dN1dxi*vfPts[0][0]   + dN2dxi*vfPts[1][0]   + dN3dxi*vfPts[2][0]   + dN4dxi*vfPts[3][0]   + dN5dxi*vfPts[4][0]   + dN6dxi*vfPts[5][0];
        /*dydxi*/   float J_01 = dN1dxi*vfPts[0][1]   + dN2dxi*vfPts[1][1]   + dN3dxi*vfPts[2][1]   + dN4dxi*vfPts[3][1]   + dN5dxi*vfPts[4][1]   + dN6dxi*vfPts[5][1];
        /*dzdxi*/   float J_02 = dN1dxi*vfPts[0][2]   + dN2dxi*vfPts[1][2]   + dN3dxi*vfPts[2][2]   + dN4dxi*vfPts[3][2]   + dN5dxi*vfPts[4][2]   + dN6dxi*vfPts[5][2];
        /*dxdeta*/  float J_10 = dN1deta*vfPts[0][0]  + dN2deta*vfPts[1][0]  + dN3deta*vfPts[2][0]  + dN4deta*vfPts[3][0]  + dN5deta*vfPts[4][0]  + dN6deta*vfPts[5][0];
        /*dydeta*/  float J_11 = dN1deta*vfPts[0][1]  + dN2deta*vfPts[1][1]  + dN3deta*vfPts[2][1]  + dN4deta*vfPts[3][1]  + dN5deta*vfPts[4][1]  + dN6deta*vfPts[5][1];
        /*dzdeta*/  float J_12 = dN1deta*vfPts[0][2]  + dN2deta*vfPts[1][2]  + dN3deta*vfPts[2][2]  + dN4deta*vfPts[3][2]  + dN5deta*vfPts[4][2]  + dN6deta*vfPts[5][2];
        /*dxdzeta*/ float J_20 = dN1dzeta*vfPts[0][0] + dN2dzeta*vfPts[1][0] + dN3dzeta*vfPts[2][0] + dN4dzeta*vfPts[3][0] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0];
        /*dydzeta*/ float J_21 = dN1dzeta*vfPts[0][1] + dN2dzeta*vfPts[1][1] + dN3dzeta*vfPts[2][1] + dN4dzeta*vfPts[3][1] + dN5dzeta*vfPts[4][1] + dN6dzeta*vfPts[5][1];
        /*dzdzeta*/ float J_22 = dN1dzeta*vfPts[0][2] + dN2dzeta*vfPts[1][2] + dN3dzeta*vfPts[2][2] + dN4dzeta*vfPts[3][2] + dN5dzeta*vfPts[4][2] + dN6dzeta*vfPts[5][2];

        float Jac = J_00*J_11*J_22 + J_01*J_12*J_20 + J_10*J_21*J_02 - J_20*J_11*J_02 - J_10*J_01*J_22 - J_21*J_12*J_00;

        //cout << "Jac = " << Jac << endl;

        float J1_00 = (+1) * ( (J_11*J_22) - (J_21*J_12) ) / Jac;     float J1_01 = (-1) * ( (J_01*J_22) - (J_21*J_02) ) / Jac;     float J1_02 = (-1) * ( (J_01*J_12) - (J_11*J_02) ) / Jac;
        float J1_10 = (-1) * ( (J_10*J_22) - (J_20*J_12) ) / Jac;     float J1_11 = (-1) * ( (J_00*J_22) - (J_20*J_02) ) / Jac;     float J1_12 = (-1) * ( (J_00*J_12) - (J_10*J_02) ) / Jac;
        float J1_20 = (+1) * ( (J_10*J_21) - (J_20*J_11) ) / Jac;     float J1_21 = (-1) * ( (J_00*J_21) - (J_20*J_01) ) / Jac;     float J1_22 = (-1) * ( (J_00*J_11) - (J_10*J_01) ) / Jac;

        float dN1dx = J1_00*dN1dxi + J1_01*dN1deta + J1_02*dN1dzeta;     float dN1dy = J1_10*dN1dxi + J1_11*dN1deta + J1_12*dN1dzeta;     float dN1dz = J1_20*dN1dxi + J1_21*dN1deta + J1_22*dN1dzeta;
        float dN2dx = J1_00*dN2dxi + J1_01*dN2deta + J1_02*dN2dzeta;     float dN2dy = J1_10*dN2dxi + J1_11*dN2deta + J1_12*dN2dzeta;     float dN2dz = J1_20*dN2dxi + J1_21*dN2deta + J1_22*dN2dzeta;
        float dN3dx = J1_00*dN3dxi + J1_01*dN3deta + J1_02*dN3dzeta;     float dN3dy = J1_10*dN3dxi + J1_11*dN3deta + J1_12*dN3dzeta;     float dN3dz = J1_20*dN3dxi + J1_21*dN3deta + J1_22*dN3dzeta;
        float dN4dx = J1_00*dN4dxi + J1_01*dN4deta + J1_02*dN4dzeta;     float dN4dy = J1_10*dN4dxi + J1_11*dN4deta + J1_12*dN4dzeta;     float dN4dz = J1_20*dN4dxi + J1_21*dN4deta + J1_22*dN4dzeta;
        float dN5dx = J1_00*dN5dxi + J1_01*dN5deta + J1_02*dN5dzeta;     float dN5dy = J1_10*dN5dxi + J1_11*dN5deta + J1_12*dN5dzeta;     float dN5dz = J1_20*dN5dxi + J1_21*dN5deta + J1_22*dN5dzeta;
        float dN6dx = J1_00*dN6dxi + J1_01*dN6deta + J1_02*dN6dzeta;     float dN6dy = J1_10*dN6dxi + J1_11*dN6deta + J1_12*dN6dzeta;     float dN6dz = J1_20*dN6dxi + J1_21*dN6deta + J1_22*dN6dzeta;

        float B[6][18] = {  dN1dx ,  0.0  ,  0.0  , dN2dx ,  0.0  ,  0.0  , dN3dx ,  0.0  ,  0.0  , dN4dx ,  0.0  ,  0.0  , dN5dx ,  0.0  ,  0.0  , dN6dx ,  0.0  ,  0.0  ,
                            0.0  , dN1dy ,  0.0  ,  0.0  , dN2dy ,  0.0  ,  0.0  , dN3dy ,  0.0  ,  0.0  , dN4dy ,  0.0  ,  0.0  , dN5dy ,  0.0  ,  0.0  , dN6dy ,  0.0  ,
                            0.0  ,  0.0  , dN1dz ,  0.0  ,  0.0  , dN2dz ,  0.0  ,  0.0  , dN3dz ,  0.0  ,  0.0  , dN4dz ,  0.0  ,  0.0  , dN5dz ,  0.0  ,  0.0  , dN6dz ,
                            dN1dy , dN1dx ,  0.0  , dN2dy , dN2dx ,  0.0  , dN3dy , dN3dx ,  0.0  , dN4dy , dN4dx ,  0.0  , dN5dy , dN5dx ,  0.0  , dN6dy , dN6dx ,  0.0  ,
                            dN1dz ,  0.0  , dN1dx , dN2dz ,  0.0  , dN2dx , dN3dz ,  0.0  , dN3dx , dN4dz ,  0.0  , dN4dx , dN5dz ,  0.0  , dN5dx , dN6dz ,  0.0  , dN6dx ,
                            0.0  , dN1dz , dN1dy ,  0.0  , dN2dz , dN2dy ,  0.0  , dN3dz , dN3dy ,  0.0  , dN4dz , dN4dy ,  0.0  , dN5dz , dN5dy ,  0.0  , dN6dz , dN6dy  };

        float BtD[18][6] = {0.0};

        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<6; j++)
                BtD[i][j] = B[0][i]*D[0][j] + B[1][i]*D[1][j] + B[2][i]*D[2][j] + B[3][i]*D[3][j] + B[4][i]*D[4][j] + B[5][i]*D[5][j];

        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++){
                float vBtDBaux = BtD[i][0]*B[0][j] + BtD[i][1]*B[1][j] + BtD[i][2]*B[2][j] + BtD[i][3]*B[3][j] + BtD[i][4]*B[4][j] + BtD[i][5]*B[5][j];
                vBtDB[i][j] += vBtDBaux * Jac;
            }
    
    }

    return vBtDB;
}


bool FEA2::MatrixAssemblyC3D8(int nMode) {
    if (nMode == 1) {
        int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
        Ksize = 3*nTotalNodes;
        if (bDebugMode) cout << "                - MatAssembly (hex, 24 DoF per el. Curr: " << quads_t.size() << " quads, " << nTotalNodes << " nodes, Ksize = " << Ksize << ")" << endl;
        
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

            vector<vector<float> > Kei = ComputeKeiC3D8(vfPts);

            vector<int> mn;
            for (unsigned j=0; j<nodes.size(); j++) {
                mn.push_back(nodes[j]*3);
            }

            for (unsigned int ni=0; ni<nodes.size(); ni++){
                for (unsigned int nj=0; nj<nodes.size(); nj++){
                    for (unsigned int m=0; m<3; m++){
                        for (unsigned int n=0; n<3; n++){
                            if ( (mn[ni]+m)>=Ksize || (mn[nj]+n)>=Ksize )   continue;
                            K[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                        }
                    }
                }
            }
        }
        return true;
    }
    else if (nMode == 2) {
        int nTotalNodes = vMPsXYZN_ut.size() + vMPsXYZN_ut2.size();
        Kusize = 3*nTotalNodes;
        if (bDebugMode) cout << "                - MatAssembly (hex, 24 DoF per el. Curr: " << quads_u.size() << " quads, " << nTotalNodes << " nodes, Kusize = " << Kusize << ")" << endl;
        
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

            vector<vector<float> > Kei = ComputeKeiC3D8(vfPts);

            vector<int> mn;
            for (unsigned j=0; j<nodes.size(); j++) {
                mn.push_back(nodes[j]*3);
            }

            for (unsigned int ni=0; ni<nodes.size(); ni++){
                for (unsigned int nj=0; nj<nodes.size(); nj++){
                    for (unsigned int m=0; m<3; m++){
                        for (unsigned int n=0; n<3; n++){
                            if ( (mn[ni]+m)>=Ksize || (mn[nj]+n)>=Ksize )   continue;
                            Ku[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                        }
                    }
                }
            }
        }
        return true;
    }
    else{
        return false;
    }
}


bool FEA2::MatrixAssemblyC3D6(int nMode) {
    if (nMode == 1) {
        int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
        Ksize = 3*nTotalNodes;
        if (bDebugMode) cout << "                - MatAssembly (tri, 18 DoF per el. Curr: " << triangles_t.size() << " triangles, " << nTotalNodes << " nodes, Ksize = " << Ksize << ")" << endl;

        if (Ksize<=3)
            return false;

        K = vector<vector<float> >(Ksize,vector<float>(Ksize,0.0));

        //cout << " 0 " << endl;

        for (unsigned int i=0; i<triangles_t.size(); i++) {
            vector<int> nodes;
            nodes.push_back(triangles_t[i][0]);
            nodes.push_back(triangles_t[i][1]);
            nodes.push_back(triangles_t[i][2]);
            nodes.push_back(triangles_t[i][0] + vMPsXYZN_t.size());
            nodes.push_back(triangles_t[i][1] + vMPsXYZN_t.size());
            nodes.push_back(triangles_t[i][2] + vMPsXYZN_t.size());

            //cout << " 1 " << endl;

            vector<vector<float> > vfPts;
            for (unsigned int j=0; j<3; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_t[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_t[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_t[nodes[j]][2]);
                vfPts.push_back(vfPtsi);
            }

            //cout << " 2 " << endl;

            for (unsigned int j=0; j<3; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][2]);
                vfPts.push_back(vfPtsi);
            }

            //cout << " 3 " << endl;

            vector<vector<float> > Kei = ComputeKeiC3D6(vfPts);

            //cout << " 4 " << endl;

            vector<int> mn;
            for (unsigned j=0; j<nodes.size(); j++) {
                mn.push_back(nodes[j]*3);
            }

            //cout << "pre-k" << endl;
            for (unsigned int ni=0; ni<nodes.size(); ni++){
                for (unsigned int nj=0; nj<nodes.size(); nj++){
                    for (unsigned int m=0; m<3; m++){
                        for (unsigned int n=0; n<3; n++){
                            if ( (mn[ni]+m)>=Ksize || (mn[nj]+n)>=Ksize )   continue;
                            K[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                        }
                    }
                }
            }
            //cout << "6" << endl;
        }
        return true;
    }
    else if (nMode == 2) {
        int nTotalNodes = vMPsXYZN_ut.size() + vMPsXYZN_ut2.size();
        Kusize = 3*nTotalNodes;
        if (bDebugMode) cout << "                - MatAssembly (tri, 18 DoF per el. Curr: " << triangles_u.size() << " triangles, " << nTotalNodes << " nodes, Kusize = " << Kusize << ")" << endl;

        if (Kusize<=3)
            return false;

        Ku.clear();
        vector<float> ki = vector<float>(Kusize,0.0);
        for (unsigned int i=0; i<Kusize; i++)
            Ku.push_back(ki);

        for (unsigned int i=0; i<triangles_u.size(); i++) {
            vector<int> nodes;
            nodes.push_back(triangles_u[i][0]);
            nodes.push_back(triangles_u[i][1]);
            nodes.push_back(triangles_u[i][2]);
            nodes.push_back(triangles_u[i][0] + vMPsXYZN_ut.size());
            nodes.push_back(triangles_u[i][1] + vMPsXYZN_ut.size());
            nodes.push_back(triangles_u[i][2] + vMPsXYZN_ut.size());

            vector<vector<float> > vfPts;
            for (unsigned int j=0; j<3; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][2]);
                vfPts.push_back(vfPtsi);
            }

            for (unsigned int j=0; j<3; j++) {
                vector<float> vfPtsi;
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][0]);
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][1]);
                vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][2]);
                vfPts.push_back(vfPtsi);
            }

            vector<vector<float> > Kei = ComputeKeiC3D6(vfPts);

            vector<int> mn;
            for (unsigned j=0; j<nodes.size(); j++) {
                mn.push_back(nodes[j]*3);
            }

            for (unsigned int ni=0; ni<nodes.size(); ni++){
                for (unsigned int nj=0; nj<nodes.size(); nj++){
                    for (unsigned int m=0; m<3; m++){
                        for (unsigned int n=0; n<3; n++){
                            if ( (mn[ni]+m)>=Ksize || (mn[nj]+n)>=Ksize )   continue;
                            Ku[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                        }
                    }
                }
            }
        }
        return true;
    }
    else{
        return false;
    }
}



void FEA2::ImposeDirichletEncastre_K(int nMode, vector<vector<int> > vD, float Klarge){
    for (unsigned int i=0; i<vD.size(); i++){
        int mp0 = 3*(vD[i][0] - 1);
        int mp1 = mp0 + 1;
        int mp2 = mp0 + 2;

        if (nMode==1){
            K[mp0][mp0] = Klarge;
            K[mp1][mp1] = Klarge;
            K[mp2][mp2] = Klarge;
        }
        else if (nMode==2){
            Ku[mp0][mp0] = Klarge;
            Ku[mp1][mp1] = Klarge;
            Ku[mp2][mp2] = Klarge;
        }
        /*
        K[mp0][mp0] = Klarge;
        K[mp1][mp1] = Klarge;
        K[mp2][mp2] = Klarge;
        */
    }
}


void FEA2::ImposeDirichletEncastre_a(vector<vector<int> > vD, float Klarge){
    for (unsigned int i=0; i<vD.size(); i++){
        int mp0 = 3*(vD[i][0] - 1);
        int mp1 = mp0 + 1;
        int mp2 = mp0 + 2;

        //K[mp0][mp0] = Klarge;
        //K[mp1][mp1] = Klarge;
        //K[mp2][mp2] = Klarge;

        vva[mp0][0] = 1/Klarge;
        vva[mp1][0] = 1/Klarge;
        vva[mp2][0] = 1/Klarge;
    }
}


vector<vector<float> > FEA2::InvertMatrixEigen(vector<vector<float> > m1) {
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


vector<vector<float> > FEA2::MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2) {
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


void FEA2::Set_uf(vector<vector<float> > vPoints){
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


    /*
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
    */


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


void FEA2::ComputeDisplacement(){
    vva.clear();
    for (unsigned int i=0; i<u0.size(); i++){
        vector<float> va;
        va.push_back(uf[i]-u0[i]);
        vva.push_back(va);
    }

    ImposeDirichletEncastre_a(vvDir_t,100000000.0);
}


void FEA2::ComputeForces(){
    vvf.clear();
    vvf = MultiplyMatricesEigen(K,vva);
}


float FEA2::ComputeStrainEnergy(){
    // sE = a' · K · a = a' · F

    vector<vector<float> > vvat;
    vector<float> vvati;
    for (unsigned int i=0; i<vva.size(); i++){
        vvati.push_back(vva[i][0]);
    }
    vvat.push_back(vvati);

    /*
    cout << "vvat" << endl;
    for (unsigned int i=0; i<vvat.size(); i++){
        cout << vvat[i][0] << " ";
    }
    cout << endl << endl;

    cout << "vvf" << endl;
    for (unsigned int i=0; i<vvf.size(); i++){
        cout << vvf[i][0] << " ";
    }
    cout << endl << endl;
    */

    vector<vector<float> > vvsE = MultiplyMatricesEigen(vvat,vvf);
    sE = vvsE[0][0];


    //cout << "sE  " << vvsE[0][0] << endl;

    if (sE < 0.0)
        sE = -sE;

    CurrentSE = sE;
    return sE;
}


float FEA2::NormalizeStrainEnergy(){
    //float invNormFactor = 1/fNormFactor;
    //nsE = sE*invNormFactor;

    int nEl = Ksize/3;
    nsE = sE / nEl;

    return nsE;
}


void FEA2::UpdateForces(){
    //cout << "vvf  " << vvf.size() << endl
    //     << "     " << vMPsXYZN_t.size() << "  " << vMPsXYZN_t2.size() << endl
    //     << "     " << vMPsXYZN_ut.size() << "  " << vMPsXYZN_ut2.size() << endl
    //     << "     " << vpMPs_t.size() << "  " << vpMPs_ut.size() << endl;

    unsigned int nadd = 3*(vMPsXYZN_u.size() + vMPsXYZN_u2.size() - vMPsXYZN_t.size() - vMPsXYZN_t2.size());
    for (unsigned int i=0; i<nadd; i++){
        vector<float> vfi = vector<float>(1,0.0);
        vvf.push_back(vfi);
    }
}


void FEA2::ComputeNewDisplacement(){
    //cout << "vva2  " << vva2.size() << endl;

    vva2 = MultiplyMatricesEigen(Ku,vvf);
    vva2 = vector_resize_cols(vva2,3);
}


vector<vector<float> > FEA2::vector_resize_cols(vector<vector<float> > v1, unsigned int n){
    vector<vector<float> > v2;
    vector<float> v2i;
    for (unsigned int i=0; i<v1.size(); i++){
        for (unsigned int j=0; j<v1[i].size(); j++){
            v2i.push_back(v1[i][j]);
            if (v2i.size() == n){
                v2.push_back(v2i);
                v2i.clear();
            }
        }
    }
    return v2;
}





float FEA2::GetStrainEnergy() { return sE; }

float FEA2::GetNormalizedStrainEnergy() { return nsE; }

void FEA2::setbfea2(bool bSet) { bInFEA2 = bSet; }

void FEA2::setbfea(bool bSet) { bInFEA = bSet; }

//void FEA2::setCurrEdge(int input) { nCurrEdge = input; }
