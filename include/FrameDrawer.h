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

#ifndef FRAMEDRAWER_H
#define FRAMEDRAWER_H

#include "Tracking.h"
#include "MapPoint.h"
#include "Map.h"
#include "Statistics.h"

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include<mutex>


namespace ORB_SLAM2
{

class Tracking;
class Viewer;

class FrameDrawer
{
public:

    FrameDrawer(Map* pMap);

    // Update info from the last processed frame.
    void Update(Tracking *pTracker);

    // Draw last processed frame.
    cv::Mat DrawFrame(bool bDrawMesh);

public:

    vector<vector<cv::KeyPoint*> > vpKPs2Draw;
    vector<vector<MapPoint*> > vpMPs2Draw;
    vector<float> vpMPs2DrawWgt;

protected:

    void DrawTextInfo(cv::Mat &im, int nState, cv::Mat &imText);

    cv::Point2f DistortMapPoint(MapPoint* pMP);
    cv::Point2f UndistortPoint(cv::Point2f pt);
    bool InCircle(cv::Point2f pt, int r);

    Mat SetTransparentColor(Mat &img, vector<vector<Point> > &roi, double alpha1, double alpha2); 
    cv::Scalar SetColor(float floatValue);
    // Info of the frame to be drawn
    cv::Mat mIm;
    int N;
    vector<cv::KeyPoint> mvCurrentKeys;
    vector<float> fvMapPointSearchRadious;
    vector<bool> mvbMap, mvbVO;
    vector<bool> mvbNewMP, mvbOutlier, mvbNoClose, mvbNoSimilar;
    vector<bool> mvbKP2Draw;

    vector<cv::Point2f> vp2fNoClose;
    vector<cv::Point2f> vp2fNoSimilar;

    vector<bool> mvbFullMap;
    vector<cv::Point2f> vp2fFullMap;
    vector<int> viMapAge;

    bool mbOnlyTracking;
    int mnTracked, mnTrackedVO;
    vector<cv::KeyPoint> mvIniKeys;
    vector<int> mvIniMatches;
    int mState;

    // Frame data
    cv::Mat mRcw, mtcw;
    float fx, fy, cx, cy, k1, k2, p1, p2, k3;

    Map* mpMap;
    Statistics* pStats1;

    std::mutex mMutex;
    std::mutex mMapPointMutex;
};

} //namespace ORB_SLAM

#endif // FRAMEDRAWER_H
