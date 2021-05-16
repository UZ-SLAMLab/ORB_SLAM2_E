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

#include "FrameDrawer.h"
#include "Tracking.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include<mutex>

namespace ORB_SLAM2
{

FrameDrawer::FrameDrawer(Map* pMap):mpMap(pMap)
{
    mState=Tracking::SYSTEM_NOT_READY;
    mIm = cv::Mat(480,640,CV_8UC3, cv::Scalar(0,0,0));

    // Output text files for saving statistics
    static Statistics Stats1("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsW1.txt");
      pStats1 = &Stats1;
      pStats1->OpenFile(1);
      pStats1->ColumnHeadersT();
      pStats1->CloseFile();

}



cv::Mat FrameDrawer::DrawFrame(bool bDrawMesh)
{
    cv::Mat im;
    vector<cv::KeyPoint> vIniKeys;              // Initialization: KeyPoints in reference frame
    vector<int> vMatches;                       // Initialization: correspondeces with reference keypoints
    vector<cv::KeyPoint> vCurrentKeys;          // KeyPoints in current frame
    vector<bool> vbVO, vbMap;                   // Tracked MapPoints in current frame
    int state;                                  // Tracking state

    vector<cv::Point2f> vMap;
    vector<cv::Point2f> vNoClose;
    vector<cv::Point2f> vNoSimilar;
    vector<float> fRadious;

    vector<bool> vbFullMap, vbNewMP, vbOutlier, vbNoClose, vbNoSimilar;

    vector<int> viAge;

    int nMPsInFrame = 0;
    int nTrackedYoung = 0;
    int nTrackedAdult = 0;
    int nOutliers = 0;
    int nNoCloseORB = 0;
    int nNoSimilar = 0;

    //Copy variables within scoped mutex
    {
        unique_lock<mutex> lock(mMutex);
        state=mState;
        if(mState==Tracking::SYSTEM_NOT_READY)
            mState=Tracking::NO_IMAGES_YET;

        mIm.copyTo(im);

        if(mState==Tracking::NOT_INITIALIZED)
        {
            vCurrentKeys = mvCurrentKeys;
            vIniKeys = mvIniKeys;
            vMatches = mvIniMatches;
        }
        else if(mState==Tracking::OK)
        {
            vCurrentKeys = mvCurrentKeys;
            vbVO = mvbVO;
            vbMap = mvbMap;

            vbFullMap = mvbFullMap;
            vbNewMP = mvbNewMP;
            vbOutlier = mvbOutlier;
            vbNoClose = mvbNoClose;
            vbNoSimilar = mvbNoSimilar;

            vMap = vp2fFullMap;
            fRadious = fvMapPointSearchRadious;
            vNoClose = vp2fNoClose;
            vNoSimilar = vp2fNoSimilar;

            viAge = viMapAge;
        }
        else if(mState==Tracking::LOST)
        {
            vCurrentKeys = mvCurrentKeys;
            vMap = vp2fFullMap;
        }
    } // destroy scoped mutex -> release mutex

    if(im.channels()<3) //this should be always true
        cvtColor(im,im,CV_GRAY2BGR);

    //Draw
    if(state==Tracking::NOT_INITIALIZED) //INITIALIZING
    {
        for(unsigned int i=0; i<vMatches.size(); i++)
        {
            if(vMatches[i]>=0)
            {
                cv::line(im,vIniKeys[i].pt,vCurrentKeys[vMatches[i]].pt,
                        cv::Scalar(0,255,0));
            }
        }
    }
    else if(state==Tracking::OK) //TRACKING
    {
        mnTracked=0;
        mnTrackedVO=0;
        for(int i=0;i<N;i++)
        {
            if (vbFullMap[i])
            {
                if((vbVO[i] || vbMap[i]) && !vbOutlier[i] && !vbNoClose[i] && !vbNoSimilar[i])
                {
                    cv::Point2f p2f_mp = vMap[i];
                    cv::Point2f p2f_kp = vCurrentKeys[i].pt;

                    if(viAge[i]<3)
                    {
                        cv::circle(im,p2f_kp,2,cv::Scalar(255,190,0),-1);
                        cv::line(im,p2f_kp,p2f_mp,cv::Scalar(255,120,20),1,8,0);
                        cv::circle(im,p2f_mp,2,cv::Scalar(255,120,20),-1);
                        if (fRadious[i]>0.0) cv::circle(im,p2f_mp,fRadious[i],cv::Scalar(255,120,20),1);
                        nTrackedYoung++;
                    }
                    else
                    {
                        cv::circle(im,p2f_kp,2,cv::Scalar(0,255,0),-1);
                        cv::line(im,p2f_kp,p2f_mp,cv::Scalar(0,120,0),1,8,0);
                        cv::circle(im,p2f_mp,2,cv::Scalar(0,120,0),-1);
                        if (fRadious[i]>0.0) cv::circle(im,p2f_mp,fRadious[i],cv::Scalar(0,120,0),1);
                        nTrackedAdult++;
                    }
                    mnTracked++;
                    nMPsInFrame++;
                }

                else if(vbOutlier[i] && !vbNoClose[i] && !vbNoSimilar[i])
                {
                    cv::Point2f p2f_mp = vMap[i];
                    cv::Point2f p2f_kp = vCurrentKeys[i].pt;

                    cv::circle(im,p2f_kp,2,cv::Scalar(200,0,150),-1);
                    cv::line(im,p2f_kp,p2f_mp,cv::Scalar(150,0,100),1,8,0);
                    cv::circle(im,p2f_mp,2,cv::Scalar(150,0,100),-1);
                    if (fRadious[i]>0.0) cv::circle(im,p2f_mp,fRadious[i],cv::Scalar(150,0,100),1);
                    nOutliers++;
                    nMPsInFrame++;
                }

                else if(!vbOutlier[i] && vbNoClose[i] && !vbNoSimilar[i])
                {
                    cv::Point2f p2f_mp = vMap[i];
                    cv::circle(im,p2f_mp,2,cv::Scalar(140,0,0),-1);
                    nNoCloseORB++;
                    nMPsInFrame++;
                }

                else if (!vbOutlier[i] && !vbNoClose[i] && vbNoSimilar[i])
                {
                    cv::Point2f p2f_mp = vMap[i];
                    cv::circle(im,p2f_mp,2,cv::Scalar(140,50,50),-1);
                    nNoSimilar++;
                    nMPsInFrame++;
                }
            }
        }
        
        if (vpMPs2Draw.size() > 0 && bDrawMesh==true && mRcw.type()==5 && mtcw.type()==5){

            vector<vector<Point> > roi;

            for(unsigned int i=0; i<vpMPs2Draw.size(); i++){
                MapPoint* pMP0 = vpMPs2Draw[i][0];
                MapPoint* pMP1 = vpMPs2Draw[i][1];
                MapPoint* pMP2 = vpMPs2Draw[i][2];

                if (pMP0 && pMP1 && pMP2) {
                    cv::Point v0 = DistortMapPoint(pMP0);
                    cv::Point v1 = DistortMapPoint(pMP1);
                    cv::Point v2 = DistortMapPoint(pMP2);

                    int ncolor = 220;
                    cv::line(im,v0,v1,cv::Scalar(ncolor,ncolor,ncolor),1,8,0);
                    cv::line(im,v0,v2,cv::Scalar(ncolor,ncolor,ncolor),1,8,0);
                    cv::line(im,v1,v2,cv::Scalar(ncolor,ncolor,ncolor),1,8,0);

                    vector<Point> roii;
                    roii.push_back(v0);
                    roii.push_back(v1);
                    roii.push_back(v2);
                    roi.push_back(roii);
                }
            }
            if (vpMPs2DrawWgt.size()>0)
                im = SetTransparentColor(im, roi, 1.0, 0.5);
        }
    }

    pStats1->OpenFile(0);
    pStats1->ColumnHeadersT();
    pStats1->CloseFile();

    cv::Mat imWithInfo;
    DrawTextInfo(im,state, imWithInfo);

    //cout << "Outliers " << nOutliers << endl;

    return imWithInfo;
}


void FrameDrawer::DrawTextInfo(cv::Mat &im, int nState, cv::Mat &imText)
{
    stringstream s;
    if(nState==Tracking::NO_IMAGES_YET)
        s << " WAITING FOR IMAGES";
    else if(nState==Tracking::NOT_INITIALIZED)
        s << " TRYING TO INITIALIZE ";
    else if(nState==Tracking::OK)
    {
        if(!mbOnlyTracking)
            s << "SLAM MODE |  ";
        else
            s << "LOCALIZATION | ";
        int nKFs = mpMap->KeyFramesInMap();
        int nMPs = mpMap->MapPointsInMap();
        s << "KFs: " << nKFs << ", MPs: " << nMPs << ", Matches: " << mnTracked;
        if(mnTrackedVO>0)
            s << ", + VO matches: " << mnTrackedVO;
    }
    else if(nState==Tracking::LOST)
    {
        s << " TRACK LOST. TRYING TO RELOCALIZE ";
    }
    else if(nState==Tracking::SYSTEM_NOT_READY)
    {
        s << " LOADING ORB VOCABULARY. PLEASE WAIT...";
    }

    int baseline=0;
    cv::Size textSize = cv::getTextSize(s.str(),cv::FONT_HERSHEY_PLAIN,1,1,&baseline);

    imText = cv::Mat(im.rows+textSize.height+10,im.cols,im.type());
    im.copyTo(imText.rowRange(0,im.rows).colRange(0,im.cols));
    imText.rowRange(im.rows,imText.rows) = cv::Mat::zeros(textSize.height+10,im.cols,im.type());
    cv::putText(imText,s.str(),cv::Point(5,imText.rows-5),cv::FONT_HERSHEY_PLAIN,1,cv::Scalar(255,255,255),1,8);

}

void FrameDrawer::Update(Tracking *pTracker)
{
    unique_lock<mutex> lock(mMutex);
    if (pTracker->colorPub == 1)
        pTracker->mImColor.copyTo(mIm);
    else
        pTracker->mImGray.copyTo(mIm);
    mvCurrentKeys=pTracker->mCurrentFrame.mvKeys;
    N = mvCurrentKeys.size();
    mvbVO = vector<bool>(N,false);
    mvbMap = vector<bool>(N,false);
    mvbFullMap = vector<bool>(N,false);
    mvbNewMP = vector<bool>(N,false);
    vp2fFullMap = vector<cv::Point2f>(N,cv::Point2f(0.0,0.0));
    viMapAge = vector<int>(N,0);
    mvbOutlier = vector<bool>(N,false);
    mvbNoClose = vector<bool>(N,false);
    vp2fNoClose = vector<cv::Point2f>(N,cv::Point2f(0.0,0.0));
    mvbNoSimilar = vector<bool>(N,false);
    vp2fNoSimilar = vector<cv::Point2f>(N,cv::Point2f(0.0,0.0));
    mbOnlyTracking = pTracker->mbOnlyTracking;

    fvMapPointSearchRadious = pTracker->mCurrentFrame.fvMapPointSearchRadious;

    // Frame data
    mRcw = pTracker->mCurrentFrame.Get_mRcw();
    mtcw = pTracker->mCurrentFrame.Get_mtcw();
    fx = pTracker->mCurrentFrame.fx;
    fy = pTracker->mCurrentFrame.fy;
    cx = pTracker->mCurrentFrame.cx;
    cy = pTracker->mCurrentFrame.cy;

    cv::Mat DistCoef(4,1,CV_32F);
    DistCoef = pTracker->GetDistCoef();
    k1 = DistCoef.at<float>(0);
    k2 = DistCoef.at<float>(1);
    p1 = DistCoef.at<float>(2);
    p2 = DistCoef.at<float>(3);
    k3 = DistCoef.at<float>(4);


    if(pTracker->mLastProcessedState==Tracking::NOT_INITIALIZED)
    {
        mvIniKeys=pTracker->mInitialFrame.mvKeys;
        mvIniMatches=pTracker->mvIniMatches;
    }
    else if(pTracker->mLastProcessedState==Tracking::OK)
    {
        for(int i=0;i<N;i++)
        {
            MapPoint* pMP = pTracker->mCurrentFrame.mvpMapPoints[i];
            if(pMP)
            {
                if(!pTracker->mCurrentFrame.mvbOutlier[i])
                {
                    if(pMP->Observations()>0)
                        mvbMap[i]=true;
                    else
                        mvbVO[i]=true;
                }
            }
            MapPoint* pMP2 = pTracker->mCurrentFrame.mvpMapPointsInFrame[i];
            if(pMP2)
            {
                mvbFullMap[i] = true;
                vp2fFullMap[i] = DistortMapPoint(pMP2);
                viMapAge[i] = pMP2->MPage;
                mvbNewMP[i] = pMP2->bNewMP;
                mvbOutlier[i] = pMP2->isBad();
            }
            MapPoint* pMP3 = pTracker->mCurrentFrame.mvpMapPointsWoCloseORB[i];
            if(pMP3)
            {
                mvbNoClose[i] = true;
                vp2fNoClose[i] = DistortMapPoint(pMP3);
            }
            MapPoint* pMP4 = pTracker->mCurrentFrame.mvpMapPointsWoSimilarORB[i];
            if(pMP4)
            {
                mvbNoSimilar[i] = true;
                vp2fNoSimilar[i] = DistortMapPoint(pMP4);
            }
        }
    }
    mState=static_cast<int>(pTracker->mLastProcessedState);
}

cv::Point2f FrameDrawer::DistortMapPoint(MapPoint* pMP)
{
    unique_lock<mutex> lock(mMapPointMutex);

    // 3D in absolute coordinates
    cv::Mat P = pMP->GetWorldPos();

    // 3D in camera coordinates
    const cv::Mat Pc = mRcw*P+mtcw;
    const float PcX = Pc.at<float>(0);
    const float PcY = Pc.at<float>(1);
    const float PcZ = Pc.at<float>(2);

    //Project in image
    const float invz = 1.0/PcZ;

    float a = fx*PcX*invz + cx;
    float b = fy*PcY*invz + cy;

    float aa = (a-cx)/fx;
    float bb = (b-cy)/fy;
    float rSQ = (aa*aa+bb*bb);
    float corr = (1+rSQ*k1+rSQ*rSQ*k2);

    cv::Point2f pt;

    pt.x = (aa*corr)*fx+cx;
    pt.y = (bb*corr)*fy+cy;

    return pt;
}

cv::Point2f FrameDrawer::UndistortPoint(cv::Point2f pt)
{
    float a = (pt.x-cx)/fx;
    float b = (pt.y-cy)/fy;
    float rSQ = (a*a+b*b);
    float corr = (1+rSQ*k1+rSQ*rSQ*k2);

    cv::Point2f ptUn;

    ptUn.x = (a*corr)*fx+cx;
    ptUn.y = (b*corr)*fy+cy;

    return ptUn;
}

bool FrameDrawer::InCircle(cv::Point2f pt, int r)
{
    int r_sq = r*r;
    int pt_rx = pt.x - cx;
    int pt_ry = pt.y - cy;
    int dist = pt_rx*pt_rx + pt_ry*pt_ry;
    if (dist > r_sq)
        return false;
    else
        return true;
}


Mat FrameDrawer::SetTransparentColor(Mat &img, vector<vector<Point> > &roi, double alpha1, double alpha2) {
    Mat out;
    Mat layer = Mat::zeros(img.size(), CV_8UC3);

    for (unsigned int i=0; i<roi.size(); i++){
        float colorit = vpMPs2DrawWgt[i];
        Scalar colori = SetColor(colorit);

        vector<vector<Point> > roi1;
        roi1.push_back(roi[i]);

        fillPoly(layer,roi1,colori);
    }
    addWeighted(img, alpha1, layer, alpha2, 0, out);
    return out;
}


cv::Scalar FrameDrawer::SetColor(float floatValue){
    int cRed = floatValue * 255;
    int cBlue = 255 - floatValue * 255;
    int cGreen = 0;

    if (floatValue>=0 && floatValue<= 0.5)
       cGreen = floatValue * 512; 
    else if (floatValue > 0.5 && floatValue <= 1)
       cGreen = 255 - (floatValue - 0.5)*512;
    else 
       return -1;

    cv::Scalar output = cv::Scalar(cBlue,cGreen,cRed);
    return output;
}




} //namespace ORB_SLAM
