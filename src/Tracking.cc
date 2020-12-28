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


#include "Tracking.h"

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include"ORBmatcher.h"
#include"FrameDrawer.h"
#include"Converter.h"
#include"Map.h"
#include"Initializer.h"

#include"Optimizer.h"
#include"PnPsolver.h"

#include<iostream>

#include<mutex>


using namespace std;

namespace ORB_SLAM2
{

Tracking::Tracking(System *pSys, ORBVocabulary* pVoc, FrameDrawer *pFrameDrawer, MapDrawer *pMapDrawer, Map *pMap, KeyFrameDatabase* pKFDB, const string &strSettingPath, const int sensor):
    mState(NO_IMAGES_YET), mSensor(sensor), mbOnlyTracking(false), mbVO(false), mpORBVocabulary(pVoc),
    mpKeyFrameDB(pKFDB), mpInitializer(static_cast<Initializer*>(NULL)), mpSystem(pSys), mpViewer(NULL),
    mpFrameDrawer(pFrameDrawer), mpMapDrawer(pMapDrawer), mpMap(pMap), mnLastRelocFrameId(0)
{
    // Load camera parameters from settings file

    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    float fx = fSettings["Camera.fx"];
    float fy = fSettings["Camera.fy"];
    float cx = fSettings["Camera.cx"];
    float cy = fSettings["Camera.cy"];

    cv::Mat K = cv::Mat::eye(3,3,CV_32F);
    K.at<float>(0,0) = fx;
    K.at<float>(1,1) = fy;
    K.at<float>(0,2) = cx;
    K.at<float>(1,2) = cy;
    K.copyTo(mK);

    cv::Mat DistCoef(4,1,CV_32F);
    DistCoef.at<float>(0) = fSettings["Camera.k1"];
    DistCoef.at<float>(1) = fSettings["Camera.k2"];
    DistCoef.at<float>(2) = fSettings["Camera.p1"];
    DistCoef.at<float>(3) = fSettings["Camera.p2"];
    const float k3 = fSettings["Camera.k3"];
    if(k3!=0)
    {
        DistCoef.resize(5);
        DistCoef.at<float>(4) = k3;
    }
    DistCoef.copyTo(mDistCoef);

    mbf = fSettings["Camera.bf"];

    float fps = fSettings["Camera.fps"];
    if(fps==0)
        fps=30;

    // Max/Min Frames to insert keyframes and to check relocalisation
    mMinFrames = 0;
    mMaxFrames = fps;

    cout << endl << "Camera Parameters: " << endl;
    cout << "- fx: " << fx << endl;
    cout << "- fy: " << fy << endl;
    cout << "- cx: " << cx << endl;
    cout << "- cy: " << cy << endl;
    cout << "- k1: " << DistCoef.at<float>(0) << endl;
    cout << "- k2: " << DistCoef.at<float>(1) << endl;
    if(DistCoef.rows==5)
        cout << "- k3: " << DistCoef.at<float>(4) << endl;
    cout << "- p1: " << DistCoef.at<float>(2) << endl;
    cout << "- p2: " << DistCoef.at<float>(3) << endl;
    cout << "- fps: " << fps << endl;


    int nRGB = fSettings["Camera.RGB"];
    mbRGB = nRGB;

    if(mbRGB)
        cout << "- color order: RGB (ignored if grayscale)" << endl;
    else
        cout << "- color order: BGR (ignored if grayscale)" << endl;

    colorPub = fSettings["Camera.colorPub"];
    if (colorPub==1)
        cout << "- publishing in Color" << endl;
    else
        cout << "- publishing in Grayscale" << endl;

    // Load ORB parameters

    int nFeatures = fSettings["ORBextractor.nFeatures"];
    float fScaleFactor = fSettings["ORBextractor.scaleFactor"];
    int nLevels = fSettings["ORBextractor.nLevels"];
    int fIniThFAST = fSettings["ORBextractor.iniThFAST"];
    int fMinThFAST = fSettings["ORBextractor.minThFAST"];

    nTScaleLevels = nLevels;
    fTScaleFactor = fScaleFactor;
    mvScaleFactors.resize(nLevels);
    mvScaleFactors[0] = 1.0f;
    for (int i=1; i<nLevels; i++)
        mvScaleFactors[i] = mvScaleFactors[i-1]*fScaleFactor;

    mpORBextractorLeft = new ORBextractor(nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    if(sensor==System::STEREO)
        mpORBextractorRight = new ORBextractor(nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    if(sensor==System::MONOCULAR)
        mpIniORBextractor = new ORBextractor(2*nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    cout << endl  << "ORB Extractor Parameters: " << endl;
    cout << "- Number of Features: " << nFeatures << endl;
    cout << "- Scale Levels: " << nLevels << endl;
    cout << "- Scale Factor: " << fScaleFactor << endl;
    cout << "- Initial Fast Threshold: " << fIniThFAST << endl;
    cout << "- Minimum Fast Threshold: " << fMinThFAST << endl;

    if(sensor==System::STEREO || sensor==System::RGBD)
    {
        mThDepth = mbf*(float)fSettings["ThDepth"]/fx;
        cout << endl << "Depth Threshold (Close/Far Points): " << mThDepth << endl;
    }

    if(sensor==System::RGBD)
    {
        mDepthMapFactor = fSettings["DepthMapFactor"];
        if(fabs(mDepthMapFactor)<1e-5)
            mDepthMapFactor=1;
        else
            mDepthMapFactor = 1.0f/mDepthMapFactor;
    }

    nTestAllFrames = fSettings["RelocParam.bTestAllFrames"];
    if (nTestAllFrames==1) bTestAllFrames = true;
    nPrecisionFrames = fSettings["RelocParam.nPrecisionFrames"];
    nElType = fSettings["RelocParam.nElType"];
    if (nElType!=1 && nElType!=2) nElType = 1;
    nUseInverse = fSettings["RelocParam.bUseInverse"];
    if (nUseInverse==1) bUseInverse = true;
    cout << endl << "Reloc Settings: " << endl;
    cout << "- Test all frames: " << bTestAllFrames << endl;
    cout << "- Frames evaluated for Precision score: " << nPrecisionFrames << endl;
    cout << "- Compute Inverse in Relocalization: " << bUseInverse << endl;
    if (nElType==1) cout << "- Element Type: C3D8 - Hexahedron" << endl;
    if (nElType==2) cout << "- Element Type: C3D6 - Triangular Prism" << endl;
    cout << endl;

    // Output text files for saving statistics of PnP and NonLinearOptimization
    // One pointer for Rigid Model and one pointer for Non-Rigid Model
    static Statistics StatsReloc("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsReloc.txt");
      pStatsReloc = &StatsReloc;
      pStatsReloc->OpenFile(1);
      pStatsReloc->ColumnHeadersReloc();
      pStatsReloc->CloseFile();

    static Statistics StatsPR("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsPR.txt");
      pStatsPR = &StatsPR;
      pStatsPR->OpenFile(1);
      pStatsPR->ColumnHeadersPR();
      pStatsPR->CloseFile();

    static Statistics StatsRelocS1("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsRelocS1.txt");
      pStatsRelocS1 = &StatsRelocS1;
      pStatsRelocS1->OpenFile(1);
      pStatsRelocS1->ColumnHeadersRelocS1();
      pStatsRelocS1->CloseFile();

    static Statistics StatsRelocS2("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsRelocS2.txt");
      pStatsRelocS2 = &StatsRelocS2;
      pStatsRelocS2->OpenFile(1);
      pStatsRelocS2->CloseFile();
    static Statistics StatsRelocS3("/home/cirauqui/workspace/ORB_SLAM2/output/evaluation/StatsRelocS3.txt");
      pStatsRelocS3 = &StatsRelocS3;
      pStatsRelocS3->OpenFile(1);
      pStatsRelocS3->CloseFile();

    nRelocalizationsCounter = 0;

    // Precision and recall
    kpiTP = 0;
    kpiFP = 0;
    kpiFN = 0;
    kpiP = 0;





}

void Tracking::SetLocalMapper(LocalMapping *pLocalMapper)
{
    mpLocalMapper=pLocalMapper;
}

void Tracking::SetLoopClosing(LoopClosing *pLoopClosing)
{
    mpLoopClosing=pLoopClosing;
}

void Tracking::SetViewer(Viewer *pViewer)
{
    mpViewer=pViewer;
}


cv::Mat Tracking::GrabImageStereo(const cv::Mat &imRectLeft, const cv::Mat &imRectRight, const double &timestamp)
{
    mImGray = imRectLeft;
    mImColor = imRectLeft;
    cv::Mat imGrayRight = imRectRight;

    if(mImGray.channels()==3)
    {
        if(mbRGB)
        {
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_RGB2GRAY);
        }
        else
        {
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_BGR2GRAY);
        }
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
        {
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_RGBA2GRAY);
        }
        else
        {
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_BGRA2GRAY);
        }
    }

    mCurrentFrame = Frame(mImGray,imGrayRight,timestamp,mpORBextractorLeft,mpORBextractorRight,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);

    Track();

    return mCurrentFrame.mTcw.clone();
}


cv::Mat Tracking::GrabImageRGBD(const cv::Mat &imRGB,const cv::Mat &imD, const double &timestamp)
{
    mImGray = imRGB;
    mImColor = imRGB;
    cv::Mat imDepth = imD;

    if(mImGray.channels()==3)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
    }

    if((fabs(mDepthMapFactor-1.0f)>1e-5) || imDepth.type()!=CV_32F)
        imDepth.convertTo(imDepth,CV_32F,mDepthMapFactor);

    mCurrentFrame = Frame(mImGray,imDepth,timestamp,mpORBextractorLeft,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);

    Track();

    return mCurrentFrame.mTcw.clone();
}


cv::Mat Tracking::GrabImageMonocular(const cv::Mat &im, const double &timestamp)
{
    mImGray = im;
    mImColor = im;

    //cout << "Resolution: " << mImGray.cols << "x" << mImGray.rows << endl;

    if(mImGray.channels()==3)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
    }

    if(mState==NOT_INITIALIZED || mState==NO_IMAGES_YET)
        mCurrentFrame = Frame(mImGray,timestamp,mpIniORBextractor,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);
    else
        mCurrentFrame = Frame(mImGray,timestamp,mpORBextractorLeft,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);

    Track();

    return mCurrentFrame.mTcw.clone();
}

void Tracking::Track()
{
    //cout << " Track() " << endl;
    if(mState==NO_IMAGES_YET)
    {
        mState = NOT_INITIALIZED;
    }

    if (bMap)
    {
        bMap = BuildLoadedMap();

        if (bMap)
        {
            mState = LOST;
            cout << "OK" << endl;
        }
        else
            cout << "NOK" << endl;

        bMap = false;
    }


    mLastProcessedState=mState;

    // Get Map Mutex -> Map cannot be changed
    unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

    if(mState==NOT_INITIALIZED)
    {
        if(mSensor==System::STEREO || mSensor==System::RGBD)
            StereoInitialization();
        else
            MonocularInitialization();

        mpFrameDrawer->Update(this);

        if(mState!=OK)
            return;
    }
    else
    {
        // System is initialized. Track Frame.
        bool bOK;

        // Initial camera pose estimation using motion model or relocalization (if tracking is lost)
        if(!mbOnlyTracking)
        {
            // Local Mapping is activated. This is the normal behaviour, unless
            // you explicitly activate the "only tracking" mode.

            if(mState==OK)
            {
                // Local Mapping might have changed some MapPoints tracked in last frame
                CheckReplacedInLastFrame();

                if(mVelocity.empty() || mCurrentFrame.mnId<mnLastRelocFrameId+2)
                {
                    cout << " Track() - CurrentFrameId = RelocFrameId + 1" << endl;
                    bOK = TrackReferenceKeyFrame();
                    cout << "           TrackReferenceKeyFrame - bOK=" << bOK << endl;
                }
                else
                {
                    bOK = TrackWithMotionModel();
                    mCurrentFrame.SourceKfId = 0;
                    if(mpLocalMapper->bDataToSave)
                        mpLocalMapper->SaveDataInFrame(&mCurrentFrame);
                    if(!bOK)
                    {
                        bOK = TrackReferenceKeyFrame();
                    }

                }
            }
            else
            {
                bOK = Relocalization();
            }
        }
        else
        {
            // Localization Mode: Local Mapping is deactivated

            if(mState==LOST)
            {
                bOK = Relocalization();
            }
            else
            {
                if(!mbVO)
                {
                    // In last frame we tracked enough MapPoints in the map

                    if(!mVelocity.empty())
                    {
                        bOK = TrackWithMotionModel();
                        mCurrentFrame.SourceKfId = 0;
                        if(mpLocalMapper->bDataToSave)
                            mpLocalMapper->SaveDataInFrame(&mCurrentFrame);
                    }
                    else
                    {
                        bOK = TrackReferenceKeyFrame();
                    }
                }
                else
                {
                    // In last frame we tracked mainly "visual odometry" points.

                    // We compute two camera poses, one from motion model and one doing relocalization.
                    // If relocalization is sucessfull we choose that solution, otherwise we retain
                    // the "visual odometry" solution.

                    bool bOKMM = false;
                    bool bOKReloc = false;
                    vector<MapPoint*> vpMPsMM;
                    vector<bool> vbOutMM;
                    cv::Mat TcwMM;
                    if(!mVelocity.empty())
                    {
                        bOKMM = TrackWithMotionModel();
                        vpMPsMM = mCurrentFrame.mvpMapPoints;
                        vbOutMM = mCurrentFrame.mvbOutlier;
                        TcwMM = mCurrentFrame.mTcw.clone();
                    }
                    bOKReloc = Relocalization();

                    if(bOKMM && !bOKReloc)
                    {
                        mCurrentFrame.SetPose(TcwMM);
                        mCurrentFrame.mvpMapPoints = vpMPsMM;
                        mCurrentFrame.mvbOutlier = vbOutMM;

                        if(mbVO)
                        {
                            for(int i =0; i<mCurrentFrame.N; i++)
                            {
                                if(mCurrentFrame.mvpMapPoints[i] && !mCurrentFrame.mvbOutlier[i])
                                {
                                    mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                                }
                            }
                        }
                    }
                    else if(bOKReloc)
                    {
                        mbVO = false;
                    }

                    bOK = bOKReloc || bOKMM;
                }
            }
        }

        mCurrentFrame.mpReferenceKF = mpReferenceKF;

        // If we have an initial estimation of the camera pose and matching. Track the local map.
        if(!mbOnlyTracking)
        {
            if(bOK)
                bOK = TrackLocalMap();
        }
        else
        {
            // mbVO true means that there are few matches to MapPoints in the map. We cannot retrieve
            // a local map and therefore we do not perform TrackLocalMap(). Once the system relocalizes
            // the camera we will use the local map again.
            if(bOK && !mbVO)
                bOK = TrackLocalMap();
        }


        /// Update relocalization KPIs

        if (mnLastRelocFrameId>0){

            if (mCurrentFrame.mnId-mnLastRelocFrameId>0 && mCurrentFrame.mnId<=mnLastRelocFrameId+nPrecisionFrames && !bOK) {
                // False Positive = Track not maintained for n frames after relocalization
                pStatsPR->AddValue(mnLastRelocFrameId);
                pStatsPR->AddValue(1);
                pStatsPR->AddValue(0);
                pStatsPR->NewLine();
                kpiFP++;
                kpiTot++;
            }
            else if (mCurrentFrame.mnId==mnLastRelocFrameId+nPrecisionFrames && bOK) {
                // True Positive = Track maintained for n frames after relocalization
                pStatsPR->AddValue(mnLastRelocFrameId);
                pStatsPR->AddValue(1);
                pStatsPR->AddValue(1); //Save TP and go back to state lost
                pStatsPR->NewLine();
                kpiTP++;
                kpiTot++;
                if (bTestAllFrames) bOK = false;
            }
            else if (!bOK) {
                // False Negative = Didn't tried to relocate
                pStatsPR->AddValue(mnLastRelocFrameId);
                pStatsPR->AddValue(1);
                pStatsPR->AddValue(0);
                pStatsPR->NewLine();
                kpiFN++;
                kpiTot++;
            }

            double kpiPr = 0.0;
            int kpiTPFP = kpiTP + kpiFP;
            if (kpiTPFP>0)
                kpiPr = (float)kpiTP / (float)kpiTPFP;

            double kpiRc = 0.0;
            int kpiTPFN = kpiTP + kpiFN;
            if (kpiTPFN>0)
                kpiRc = (float)kpiTP / (float)kpiTPFN;

            kpiTot++;

            cout << "kpiTot = " << kpiTot << "\tkpiTP = " << kpiTP << "\tkpiFP = " << kpiFP << "\tkpiFN = " << kpiFN << fixed << setprecision(3) << "\tkpiPr = " << kpiPr << "\tkpiRc = " << kpiRc << endl;



        }




        if(bOK)
            mState = OK;
        else
            mState=LOST;

        // Update drawer
        mpFrameDrawer->Update(this);

        // If tracking was good, check if we insert a keyframe
        if(bOK)
        {
            // Update motion model
            if(!mLastFrame.mTcw.empty())
            {
                cv::Mat LastTwc = cv::Mat::eye(4,4,CV_32F);
                mLastFrame.GetRotationInverse().copyTo(LastTwc.rowRange(0,3).colRange(0,3));
                mLastFrame.GetCameraCenter().copyTo(LastTwc.rowRange(0,3).col(3));
                mVelocity = mCurrentFrame.mTcw*LastTwc;
            }
            else
                mVelocity = cv::Mat();

            mpMapDrawer->SetCurrentCameraPose(mCurrentFrame.mTcw);

            // Clean VO matches
            for(int i=0; i<mCurrentFrame.N; i++)
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
                if(pMP)
                    if(pMP->Observations()<1)
                    {
                        mCurrentFrame.mvbOutlier[i] = false;
                        mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                    }
            }

            // Delete temporal MapPoints
            for(list<MapPoint*>::iterator lit = mlpTemporalPoints.begin(), lend =  mlpTemporalPoints.end(); lit!=lend; lit++)
            {
                MapPoint* pMP = *lit;
                delete pMP;
            }
            mlpTemporalPoints.clear();

            // Check if we need to insert a new keyframe
            if(NeedNewKeyFrame())
                CreateNewKeyFrame();

            // We allow points with high innovation (considererd outliers by the Huber Function)
            // pass to the new keyframe, so that bundle adjustment will finally decide
            // if they are outliers or not. We don't want next frame to estimate its position
            // with those points so we discard them in the frame.
            for(int i=0; i<mCurrentFrame.N;i++)
            {
                if(mCurrentFrame.mvpMapPoints[i] && mCurrentFrame.mvbOutlier[i])
                    mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
            }
        }

        // Reset if the camera get lost soon after initialization
        if(mState==LOST)
        {
            if(mpMap->KeyFramesInMap()<=5)
            {
                cout << "Track lost soon after initialization, reseting..." << endl;
                mpSystem->Reset();
                return;
            }
        }

        if(!mCurrentFrame.mpReferenceKF)
            mCurrentFrame.mpReferenceKF = mpReferenceKF;

        mLastFrame = Frame(mCurrentFrame);
    }

    // Store frame pose information to retrieve the complete camera trajectory afterwards.
    if(!mCurrentFrame.mTcw.empty())
    {
        cv::Mat Tcr = mCurrentFrame.mTcw*mCurrentFrame.mpReferenceKF->GetPoseInverse();
        mlRelativeFramePoses.push_back(Tcr);
        mlpReferences.push_back(mpReferenceKF);
        mlFrameTimes.push_back(mCurrentFrame.mTimeStamp);
        mlbLost.push_back(mState==LOST);
    }
    else
    {
        // This can happen if tracking is lost
        mlRelativeFramePoses.push_back(mlRelativeFramePoses.back());
        mlpReferences.push_back(mlpReferences.back());
        mlFrameTimes.push_back(mlFrameTimes.back());
        mlbLost.push_back(mState==LOST);
    }

}


void Tracking::StereoInitialization()
{
    if(mCurrentFrame.N>500)
    {
        // Set Frame pose to the origin
        mCurrentFrame.SetPose(cv::Mat::eye(4,4,CV_32F));

        // Create KeyFrame
        KeyFrame* pKFini = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);

        // Insert KeyFrame in the map
        mpMap->AddKeyFrame(pKFini);

        // Create MapPoints and asscoiate to KeyFrame
        for(int i=0; i<mCurrentFrame.N;i++)
        {
            float z = mCurrentFrame.mvDepth[i];
            if(z>0)
            {
                cv::Mat x3D = mCurrentFrame.UnprojectStereo(i);
                MapPoint* pNewMP = new MapPoint(x3D,pKFini,mpMap);
                pNewMP->AddObservation(pKFini,i);
                pKFini->AddMapPoint(pNewMP,i);
                pNewMP->ComputeDistinctiveDescriptors();
                pNewMP->UpdateNormalAndDepth();
                mpMap->AddMapPoint(pNewMP);

                mCurrentFrame.mvpMapPoints[i]=pNewMP;
            }
        }

        cout << "New map created with " << mpMap->MapPointsInMap() << " points" << endl;

        mpLocalMapper->InsertKeyFrame(pKFini);

        mLastFrame = Frame(mCurrentFrame);
        mnLastKeyFrameId=mCurrentFrame.mnId;
        mpLastKeyFrame = pKFini;

        mvpLocalKeyFrames.push_back(pKFini);
        mvpLocalMapPoints=mpMap->GetAllMapPoints();
        mpReferenceKF = pKFini;
        mCurrentFrame.mpReferenceKF = pKFini;

        mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

        mpMap->mvpKeyFrameOrigins.push_back(pKFini);

        mpMapDrawer->SetCurrentCameraPose(mCurrentFrame.mTcw);

        mState=OK;
    }
}

void Tracking::MonocularInitialization()
{

    if(!mpInitializer)
    {
        // Set Reference Frame
        if(mCurrentFrame.mvKeys.size()>100)
        {
            mInitialFrame = Frame(mCurrentFrame);
            mLastFrame = Frame(mCurrentFrame);
            mvbPrevMatched.resize(mCurrentFrame.mvKeysUn.size());
            for(size_t i=0; i<mCurrentFrame.mvKeysUn.size(); i++)
                mvbPrevMatched[i]=mCurrentFrame.mvKeysUn[i].pt;         // Guardar KPs del F de referencia para emparejar contra este

            if(mpInitializer)
                delete mpInitializer;

            mpInitializer =  new Initializer(mCurrentFrame,1.0,200);

            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);

            return;
        }
    }
    else
    {
        // Try to initialize
        if((int)mCurrentFrame.mvKeys.size()<=100)		//ORB-SLAM2 init
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);
            return;
        }

        // Find correspondences
        ORBmatcher matcher(0.9,true);
        int nmatches = matcher.SearchForInitialization(mInitialFrame,mCurrentFrame,mvbPrevMatched,mvIniMatches,100);

        // Check if there are enough correspondences
        if(nmatches<30)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            return;
        }

        cv::Mat Rcw; // Current Camera Rotation
        cv::Mat tcw; // Current Camera Translation
        vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)

        if(mpInitializer->Initialize(mCurrentFrame, mvIniMatches, Rcw, tcw, mvIniP3D, vbTriangulated))
        {
            for(size_t i=0, iend=mvIniMatches.size(); i<iend;i++)
            {
                if(mvIniMatches[i]>=0 && !vbTriangulated[i])
                {
                    mvIniMatches[i]=-1;
                    nmatches--;
                }
            }

            // Set Frame Poses
            mInitialFrame.SetPose(cv::Mat::eye(4,4,CV_32F));
            cv::Mat Tcw = cv::Mat::eye(4,4,CV_32F);
            Rcw.copyTo(Tcw.rowRange(0,3).colRange(0,3));
            tcw.copyTo(Tcw.rowRange(0,3).col(3));
            mCurrentFrame.SetPose(Tcw);

            CreateInitialMapMonocular();
        }
    }
}



/*void Tracking::MonocularInitialization()
{
    if(!mpInitializer)
    {
        // Set up reference frame
        if(mCurrentFrame.mvKeys.size()>100)
        {
            mInitialFrame = Frame(mCurrentFrame);
            mLastFrame = Frame(mCurrentFrame);
            mvbPrevMatched.resize(mCurrentFrame.mvKeysUn.size());
            for(size_t i=0; i<mCurrentFrame.mvKeysUn.size(); i++)
                mvbPrevMatched[i]=mCurrentFrame.mvKeysUn[i].pt;

            if(mpInitializer)
                delete mpInitializer;

            mpInitializer =  new Initializer(mCurrentFrame,1.0,200);

            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);
        }
    }
    else
    {
        // Try to initialize
        if((int)mCurrentFrame.mvKeys.size()<=100)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);
            return;
        }

        // Find correspondences
        ORBmatcher matcher(0.9,true);
        int nmatches = matcher.SearchForInitialization(mInitialFrame,mCurrentFrame,mvbPrevMatched,mvIniMatches,100);

        // Check if there are enough correspondences
        if(nmatches<100)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            return;
        }

        cv::Mat Rcw; // Current Camera Rotation
        cv::Mat tcw; // Current Camera Translation
        vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)

        if(mpInitializer->Initialize(mCurrentFrame, mvIniMatches, Rcw, tcw, mvIniP3D, vbTriangulated))
        {
            for(size_t i=0, iend=mvIniMatches.size(); i<iend;i++)
            {
                if(mvIniMatches[i]>=0 && !vbTriangulated[i])
                {
                    mvIniMatches[i]=-1;
                    nmatches--;
                }
            }

            // Set Frame Poses
            mInitialFrame.SetPose(cv::Mat::eye(4,4,CV_32F));
            cv::Mat Tcw = cv::Mat::eye(4,4,CV_32F);
            Rcw.copyTo(Tcw.rowRange(0,3).colRange(0,3));
            tcw.copyTo(Tcw.rowRange(0,3).col(3));
            mCurrentFrame.SetPose(Tcw);

            CreateInitialMapMonocular();
        }
    }
}*/







void Tracking::CreateInitialMapMonocular()
{
    // Create KeyFrames
    KeyFrame* pKFini = new KeyFrame(mInitialFrame,mpMap,mpKeyFrameDB);
    KeyFrame* pKFcur = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);


    pKFini->ComputeBoW();
    pKFcur->ComputeBoW();

    // Insert KFs in the map
    mpMap->AddKeyFrame(pKFini);
    mpMap->AddKeyFrame(pKFcur);

    // Create MapPoints and asscoiate to keyframes
    for(size_t i=0; i<mvIniMatches.size();i++)
    {
        if(mvIniMatches[i]<0)
            continue;

        //Create MapPoint.
        cv::Mat worldPos(mvIniP3D[i]);

        MapPoint* pMP = new MapPoint(worldPos,pKFcur,mpMap);

        pKFini->AddMapPoint(pMP,i);
        pKFcur->AddMapPoint(pMP,mvIniMatches[i]);

        pMP->AddObservation(pKFini,i);
        pMP->AddObservation(pKFcur,mvIniMatches[i]);

        pMP->ComputeDistinctiveDescriptors();
        pMP->UpdateNormalAndDepth();

        //Fill Current Frame structure
        mCurrentFrame.mvpMapPoints[mvIniMatches[i]] = pMP;
        mCurrentFrame.mvbOutlier[mvIniMatches[i]] = false;

        //Add to Map
        mpMap->AddMapPoint(pMP);
    }

    // Update Connections
    pKFini->UpdateConnections();
    pKFcur->UpdateConnections();

    // Bundle Adjustment
    cout << "New Map created with " << mpMap->MapPointsInMap() << " points" << endl;

    Optimizer::GlobalBundleAdjustemnt(mpMap,20);

    // Set median depth to 1
    float medianDepth = pKFini->ComputeSceneMedianDepth(2);
    float invMedianDepth = 1.0f/medianDepth;

    if(medianDepth<0 || pKFcur->TrackedMapPoints(1)<40)		//ORB-SLAM2 init
    {
        cout << "Wrong initialization, reseting..." << endl;
        Reset();
        return;
    }

    // Scale initial baseline
    cv::Mat Tc2w = pKFcur->GetPose();
    Tc2w.col(3).rowRange(0,3) = Tc2w.col(3).rowRange(0,3)*invMedianDepth;
    pKFcur->SetPose(Tc2w);

    // Scale points
    vector<MapPoint*> vpAllMapPoints = pKFini->GetMapPointMatches();
    for(size_t iMP=0; iMP<vpAllMapPoints.size(); iMP++)
    {
        if(vpAllMapPoints[iMP])
        {
            MapPoint* pMP = vpAllMapPoints[iMP];
            pMP->SetWorldPos(pMP->GetWorldPos()*invMedianDepth);
        }
    }

    mpLocalMapper->InsertKeyFrame(pKFini);
    mpLocalMapper->InsertKeyFrame(pKFcur);

    mCurrentFrame.SetPose(pKFcur->GetPose());
    mnLastKeyFrameId=mCurrentFrame.mnId;
    mpLastKeyFrame = pKFcur;

    mvpLocalKeyFrames.push_back(pKFcur);
    mvpLocalKeyFrames.push_back(pKFini);
    mvpLocalMapPoints=mpMap->GetAllMapPoints();
    mpReferenceKF = pKFcur;
    mCurrentFrame.mpReferenceKF = pKFcur;

    mLastFrame = Frame(mCurrentFrame);

    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

    mpMapDrawer->SetCurrentCameraPose(pKFcur->GetPose());

    mpMap->mvpKeyFrameOrigins.push_back(pKFini);

    mState=OK;
}




/*void Tracking::CreateInitialMapMonocular()
{
    // Create KeyFrames
    KeyFrame* pKFini = new KeyFrame(mInitialFrame,mpMap,mpKeyFrameDB);
    KeyFrame* pKFcur = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);


    pKFini->ComputeBoW();
    pKFcur->ComputeBoW();

    // Insert KFs in the map
    mpMap->AddKeyFrame(pKFini);
    mpMap->AddKeyFrame(pKFcur);

    // Create MapPoints and asscoiate to keyframes
    for(size_t i=0; i<mvIniMatches.size();i++)
    {
        if(mvIniMatches[i]<0)
            continue;

        //Create MapPoint.
        cv::Mat worldPos(mvIniP3D[i]);

        MapPoint* pMP = new MapPoint(worldPos,pKFcur,mpMap);

        pKFini->AddMapPoint(pMP,i);
        pKFcur->AddMapPoint(pMP,mvIniMatches[i]);

        pMP->AddObservation(pKFini,i);
        pMP->AddObservation(pKFcur,mvIniMatches[i]);

        pMP->ComputeDistinctiveDescriptors();
        pMP->UpdateNormalAndDepth();

        //Fill Current Frame structure
        mCurrentFrame.mvpMapPoints[mvIniMatches[i]] = pMP;
        mCurrentFrame.mvbOutlier[mvIniMatches[i]] = false;

        //Add to Map
        mpMap->AddMapPoint(pMP);
    }

    // Update Connections
    pKFini->UpdateConnections();
    pKFcur->UpdateConnections();

    // Bundle Adjustment
    cout << "New Map created with " << mpMap->MapPointsInMap() << " points" << endl;

    Optimizer::GlobalBundleAdjustemnt(mpMap,20);

    // Set median depth to 1
    float medianDepth = pKFini->ComputeSceneMedianDepth(2);
    float invMedianDepth = 1.0f/medianDepth;

    if(medianDepth<0 || pKFcur->TrackedMapPoints(1)<100)
    {
        cout << "Wrong initialization, reseting..." << endl;
        Reset();
        return;
    }

    // Scale initial baseline
    cv::Mat Tc2w = pKFcur->GetPose();
    Tc2w.col(3).rowRange(0,3) = Tc2w.col(3).rowRange(0,3)*invMedianDepth;
    pKFcur->SetPose(Tc2w);

    // Scale points
    vector<MapPoint*> vpAllMapPoints = pKFini->GetMapPointMatches();
    for(size_t iMP=0; iMP<vpAllMapPoints.size(); iMP++)
    {
        if(vpAllMapPoints[iMP])
        {
            MapPoint* pMP = vpAllMapPoints[iMP];
            pMP->SetWorldPos(pMP->GetWorldPos()*invMedianDepth);
        }
    }

    mpLocalMapper->InsertKeyFrame(pKFini);
    mpLocalMapper->InsertKeyFrame(pKFcur);

    mCurrentFrame.SetPose(pKFcur->GetPose());
    mnLastKeyFrameId=mCurrentFrame.mnId;
    mpLastKeyFrame = pKFcur;

    mvpLocalKeyFrames.push_back(pKFcur);
    mvpLocalKeyFrames.push_back(pKFini);
    mvpLocalMapPoints=mpMap->GetAllMapPoints();
    mpReferenceKF = pKFcur;
    mCurrentFrame.mpReferenceKF = pKFcur;

    mLastFrame = Frame(mCurrentFrame);

    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

    mpMapDrawer->SetCurrentCameraPose(pKFcur->GetPose());

    mpMap->mvpKeyFrameOrigins.push_back(pKFini);

    mState=OK;
}*/








void Tracking::CheckReplacedInLastFrame()
{
    for(int i =0; i<mLastFrame.N; i++)
    {
        MapPoint* pMP = mLastFrame.mvpMapPoints[i];

        if(pMP)
        {
            MapPoint* pRep = pMP->GetReplaced();
            if(pRep)
            {
                mLastFrame.mvpMapPoints[i] = pRep;
            }
        }
    }
}


bool Tracking::TrackReferenceKeyFrame()
{
    // Compute Bag of Words vector
    mCurrentFrame.ComputeBoW();

    // We perform first an ORB matching with the reference keyframe
    // If enough matches are found we setup a PnP solver
    ORBmatcher matcher(0.7,true);
    vector<MapPoint*> vpMapPointMatches;

    int nmatches = matcher.SearchByBoW(mpReferenceKF,mCurrentFrame,vpMapPointMatches);

    if(nmatches<15)
    {
        cout << "           nmatches with TrackReferenceKeyFrame() " << nmatches << endl;
        return false;
    }

    mCurrentFrame.mvpMapPoints = vpMapPointMatches;
    mCurrentFrame.SetPose(mLastFrame.mTcw);

    Optimizer::PoseOptimization(&mCurrentFrame);

    // Discard outliers
    int nmatchesMap = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                nmatchesMap++;
        }
    }

    return nmatchesMap>=10;
}


bool Tracking::CheckReferenceKeyFrameReloc()
{
    // Compute Bag of Words vector
    mCurrentFrame.ComputeBoW();

    // We perform first an ORB matching with the reference keyframe
    // If enough matches are found we setup a PnP solver
    ORBmatcher matcher(0.7,true);
    vector<MapPoint*> vpMapPointMatches;

    int nmatches = matcher.SearchByBoW(mpReferenceKF,mCurrentFrame,vpMapPointMatches);

    if(nmatches<15)
        return false;

    mCurrentFrame.mvpMapPoints = vpMapPointMatches;
    mCurrentFrame.SetPose(mLastFrame.mTcw);

    Optimizer::PoseOptimization(&mCurrentFrame);

    // Discard outliers
    int nmatchesMap = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                nmatchesMap++;
        }
    }

    return nmatchesMap>=10;
}

void Tracking::UpdateLastFrame()
{
    // Update pose according to reference keyframe
    KeyFrame* pRef = mLastFrame.mpReferenceKF;
    cv::Mat Tlr = mlRelativeFramePoses.back();

    mLastFrame.SetPose(Tlr*pRef->GetPose());

    if(mnLastKeyFrameId==mLastFrame.mnId || mSensor==System::MONOCULAR || !mbOnlyTracking)
        return;

    // Create "visual odometry" MapPoints
    // We sort points according to their measured depth by the stereo/RGB-D sensor
    vector<pair<float,int> > vDepthIdx;
    vDepthIdx.reserve(mLastFrame.N);
    for(int i=0; i<mLastFrame.N;i++)
    {
        float z = mLastFrame.mvDepth[i];
        if(z>0)
        {
            vDepthIdx.push_back(make_pair(z,i));
        }
    }

    if(vDepthIdx.empty())
        return;

    sort(vDepthIdx.begin(),vDepthIdx.end());

    // We insert all close points (depth<mThDepth)
    // If less than 100 close points, we insert the 100 closest ones.
    int nPoints = 0;
    for(size_t j=0; j<vDepthIdx.size();j++)
    {
        int i = vDepthIdx[j].second;

        bool bCreateNew = false;

        MapPoint* pMP = mLastFrame.mvpMapPoints[i];
        if(!pMP)
            bCreateNew = true;
        else if(pMP->Observations()<1)
        {
            bCreateNew = true;
        }

        if(bCreateNew)
        {
            cv::Mat x3D = mLastFrame.UnprojectStereo(i);
            MapPoint* pNewMP = new MapPoint(x3D,mpMap,&mLastFrame,i);

            mLastFrame.mvpMapPoints[i]=pNewMP;

            mlpTemporalPoints.push_back(pNewMP);
            nPoints++;
        }
        else
        {
            nPoints++;
        }

        if(vDepthIdx[j].first>mThDepth && nPoints>100)
            break;
    }
}

bool Tracking::TrackWithMotionModel()
{
    ORBmatcher matcher(0.9,true);

    // Update last frame pose according to its reference keyframe
    // Create "visual odometry" points if in Localization Mode
    UpdateLastFrame();

    mCurrentFrame.SetPose(mVelocity*mLastFrame.mTcw);

    fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));

    // Project points seen in previous frame
    int th;
    if(mSensor!=System::STEREO)
        th=15;
    else
        th=7;
    int nmatches = matcher.SearchByProjection(mCurrentFrame,mLastFrame,th,mSensor==System::MONOCULAR);

    // If few matches, uses a wider window search
    if(nmatches<20)
    {
        fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));
        nmatches = matcher.SearchByProjection(mCurrentFrame,mLastFrame,2*th,mSensor==System::MONOCULAR);
    }

    if(nmatches<20)
        return false;

    // Optimize frame pose with all matches
    Optimizer::PoseOptimization(&mCurrentFrame);

    // Discard outliers
    int nmatchesMap = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                nmatchesMap++;
        }
    }

    if(mbOnlyTracking)
    {
        mbVO = nmatchesMap<10;
        return nmatches>20;
    }

    return nmatchesMap>=10;
}

bool Tracking::TrackLocalMap()
{
    // We have an estimation of the camera pose and some map points tracked in the frame.
    // We retrieve the local map and try to find matches to points in the local map.

    UpdateLocalMap();

    SearchLocalPoints();

    // Optimize Pose
    Optimizer::PoseOptimization(&mCurrentFrame);
    mnMatchesInliers = 0;

    // Update MapPoints Statistics
    for(int i=0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(!mCurrentFrame.mvbOutlier[i])
            {
                mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                if(!mbOnlyTracking)
                {
                    if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                        mnMatchesInliers++;
                }
                else
                    mnMatchesInliers++;
            }
            else if(mSensor==System::STEREO)
                mCurrentFrame.mvpMapPoints[i] = static_cast<MapPoint*>(NULL);

        }
    }

    // Decide if the tracking was succesful
    // More restrictive if there was a relocalization recently
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && mnMatchesInliers<50)
        return false;

    if(mnMatchesInliers<30)
        return false;
    else
        return true;
}


bool Tracking::NeedNewKeyFrame()
{
    if(mbOnlyTracking)
        return false;

    // If Local Mapping is freezed by a Loop Closure do not insert keyframes
    if(mpLocalMapper->isStopped() || mpLocalMapper->stopRequested())
        return false;

    const int nKFs = mpMap->KeyFramesInMap();

    // Do not insert keyframes if not enough frames have passed from last relocalisation
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && nKFs>mMaxFrames)
        return false;

    // Tracked MapPoints in the reference keyframe
    int nMinObs = 3;
    if(nKFs<=2)
        nMinObs=2;
    int nRefMatches = mpReferenceKF->TrackedMapPoints(nMinObs);

    // Local Mapping accept keyframes?
    bool bLocalMappingIdle = mpLocalMapper->AcceptKeyFrames();

    // Check how many "close" points are being tracked and how many could be potentially created.
    int nNonTrackedClose = 0;
    int nTrackedClose= 0;
    if(mSensor!=System::MONOCULAR)
    {
        for(int i =0; i<mCurrentFrame.N; i++)
        {
            if(mCurrentFrame.mvDepth[i]>0 && mCurrentFrame.mvDepth[i]<mThDepth)
            {
                if(mCurrentFrame.mvpMapPoints[i] && !mCurrentFrame.mvbOutlier[i])
                    nTrackedClose++;
                else
                    nNonTrackedClose++;
            }
        }
    }

    bool bNeedToInsertClose = (nTrackedClose<100) && (nNonTrackedClose>70);

    // Thresholds
    float thRefRatio = 0.75f;
    if(nKFs<2)
        thRefRatio = 0.4f;

    if(mSensor==System::MONOCULAR)
        thRefRatio = 0.9f;

    // Condition 1a: More than "MaxFrames" have passed from last keyframe insertion
    const bool c1a = mCurrentFrame.mnId>=mnLastKeyFrameId+mMaxFrames;
    // Condition 1b: More than "MinFrames" have passed and Local Mapping is idle
    const bool c1b = (mCurrentFrame.mnId>=mnLastKeyFrameId+mMinFrames && bLocalMappingIdle);
    //Condition 1c: tracking is weak
    const bool c1c =  mSensor!=System::MONOCULAR && (mnMatchesInliers<nRefMatches*0.25 || bNeedToInsertClose) ;
    // Condition 2: Few tracked points compared to reference keyframe. Lots of visual odometry compared to map matches.
    const bool c2 = ((mnMatchesInliers<nRefMatches*thRefRatio|| bNeedToInsertClose) && mnMatchesInliers>15);

    if((c1a||c1b||c1c)&&c2)
    {
        // If the mapping accepts keyframes, insert keyframe.
        // Otherwise send a signal to interrupt BA
        if(bLocalMappingIdle)
        {
            return true;
        }
        else
        {
            mpLocalMapper->InterruptBA();
            if(mSensor!=System::MONOCULAR)
            {
                if(mpLocalMapper->KeyframesInQueue()<3)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    }
    else
        return false;
}

void Tracking::CreateNewKeyFrame()
{
    if(!mpLocalMapper->SetNotStop(true))
        return;

    KeyFrame* pKF = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);

    mpReferenceKF = pKF;
    mCurrentFrame.mpReferenceKF = pKF;

    if(mSensor!=System::MONOCULAR)
    {
        mCurrentFrame.UpdatePoseMatrices();

        // We sort points by the measured depth by the stereo/RGBD sensor.
        // We create all those MapPoints whose depth < mThDepth.
        // If there are less than 100 close points we create the 100 closest.
        vector<pair<float,int> > vDepthIdx;
        vDepthIdx.reserve(mCurrentFrame.N);
        for(int i=0; i<mCurrentFrame.N; i++)
        {
            float z = mCurrentFrame.mvDepth[i];
            if(z>0)
            {
                vDepthIdx.push_back(make_pair(z,i));
            }
        }

        if(!vDepthIdx.empty())
        {
            sort(vDepthIdx.begin(),vDepthIdx.end());

            int nPoints = 0;
            for(size_t j=0; j<vDepthIdx.size();j++)
            {
                int i = vDepthIdx[j].second;

                bool bCreateNew = false;

                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
                if(!pMP)
                    bCreateNew = true;
                else if(pMP->Observations()<1)
                {
                    bCreateNew = true;
                    mCurrentFrame.mvpMapPoints[i] = static_cast<MapPoint*>(NULL);
                }

                if(bCreateNew)
                {
                    cv::Mat x3D = mCurrentFrame.UnprojectStereo(i);
                    MapPoint* pNewMP = new MapPoint(x3D,pKF,mpMap);
                    pNewMP->AddObservation(pKF,i);
                    pKF->AddMapPoint(pNewMP,i);
                    pNewMP->ComputeDistinctiveDescriptors();
                    pNewMP->UpdateNormalAndDepth();
                    mpMap->AddMapPoint(pNewMP);

                    mCurrentFrame.mvpMapPoints[i]=pNewMP;
                    nPoints++;
                }
                else
                {
                    nPoints++;
                }

                if(vDepthIdx[j].first>mThDepth && nPoints>100)
                    break;
            }
        }
    }

    mpLocalMapper->InsertKeyFrame(pKF);

    mpLocalMapper->SetNotStop(false);

    mnLastKeyFrameId = mCurrentFrame.mnId;
    mpLastKeyFrame = pKF;
}

void Tracking::SearchLocalPoints()
{
    // Do not search map points already matched
    int it = 0;
    for(vector<MapPoint*>::iterator vit=mCurrentFrame.mvpMapPoints.begin(), vend=mCurrentFrame.mvpMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP)
        {
            if(pMP->isBad())
            {
                *vit = static_cast<MapPoint*>(NULL);
            }
            else
            {
                pMP->IncreaseVisible();
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                pMP->mbTrackInView = false;

                mCurrentFrame.mvpMapPointsInFrame[it] = pMP;
                int r = RadiousByViewingCosT(pMP->mTrackViewCos);
                const int &nPredictedLevel = pMP->mnTrackScaleLevel;
                if (nPredictedLevel>=0 && nPredictedLevel<=nTScaleLevels)
                    mCurrentFrame.fvMapPointSearchRadious[it] = r*mvScaleFactors[nPredictedLevel];
            }
        }
        it++;
    }

    int nToMatch=0;

    // Project points in frame and check its visibility
    for(vector<MapPoint*>::iterator vit=mvpLocalMapPoints.begin(), vend=mvpLocalMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP->mnLastFrameSeen == mCurrentFrame.mnId)
            continue;
        if(pMP->isBad())
            continue;
        // Project (this fills MapPoint variables for matching)
        if(mCurrentFrame.isInFrustum(pMP,0.5))
        {
            pMP->IncreaseVisible();
            mCurrentFrame.mvpMapPointsInFrame[it] = pMP;
            int r = RadiousByViewingCosT(pMP->mTrackViewCos);
            const int &nPredictedLevel = pMP->mnTrackScaleLevel;
            if (nPredictedLevel>=0 && nPredictedLevel<=nTScaleLevels)
                mCurrentFrame.fvMapPointSearchRadious[it] = r*mvScaleFactors[nPredictedLevel];
            nToMatch++;
        }
    }

    if(nToMatch>0)
    {
        ORBmatcher matcher(0.8);
        int th = 1;
        if(mSensor==System::RGBD)
            th=3;
        // If the camera has been relocalised recently, perform a coarser search
        if(mCurrentFrame.mnId<mnLastRelocFrameId+2)
            th=5;
        matcher.SearchByProjection(mCurrentFrame,mvpLocalMapPoints,th);
    }
}

void Tracking::UpdateLocalMap()
{
    // This is for visualization
    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

    // Update
    UpdateLocalKeyFrames();
    UpdateLocalPoints();
}

void Tracking::UpdateLocalPoints()
{
    mvpLocalMapPoints.clear();

    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFrames.begin(), itEndKF=mvpLocalKeyFrames.end(); itKF!=itEndKF; itKF++)
    {
        KeyFrame* pKF = *itKF;
        const vector<MapPoint*> vpMPs = pKF->GetMapPointMatches();

        for(vector<MapPoint*>::const_iterator itMP=vpMPs.begin(), itEndMP=vpMPs.end(); itMP!=itEndMP; itMP++)
        {
            MapPoint* pMP = *itMP;
            if(!pMP)
                continue;
            if(pMP->mnTrackReferenceForFrame==mCurrentFrame.mnId)
                continue;
            if(!pMP->isBad())
            {
                mvpLocalMapPoints.push_back(pMP);
                pMP->mnTrackReferenceForFrame=mCurrentFrame.mnId;
            }
        }
    }
}


void Tracking::UpdateLocalKeyFrames()
{
    // Each map point vote for the keyframes in which it has been observed
    map<KeyFrame*,int> keyframeCounter;
    for(int i=0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
            if(!pMP->isBad())
            {
                const map<KeyFrame*,size_t> observations = pMP->GetObservations();
                for(map<KeyFrame*,size_t>::const_iterator it=observations.begin(), itend=observations.end(); it!=itend; it++)
                    keyframeCounter[it->first]++;
            }
            else
            {
                mCurrentFrame.mvpMapPoints[i]=NULL;
            }
        }
    }

    if(keyframeCounter.empty())
        return;

    int max=0;
    KeyFrame* pKFmax= static_cast<KeyFrame*>(NULL);

    mvpLocalKeyFrames.clear();
    mvpLocalKeyFrames.reserve(3*keyframeCounter.size());

    // All keyframes that observe a map point are included in the local map. Also check which keyframe shares most points
    for(map<KeyFrame*,int>::const_iterator it=keyframeCounter.begin(), itEnd=keyframeCounter.end(); it!=itEnd; it++)
    {
        KeyFrame* pKF = it->first;

        if(pKF->isBad())
            continue;

        if(it->second>max)
        {
            max=it->second;
            pKFmax=pKF;
        }

        mvpLocalKeyFrames.push_back(it->first);
        pKF->mnTrackReferenceForFrame = mCurrentFrame.mnId;
    }


    // Include also some not-already-included keyframes that are neighbors to already-included keyframes
    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFrames.begin(), itEndKF=mvpLocalKeyFrames.end(); itKF!=itEndKF; itKF++)
    {
        // Limit the number of keyframes
        if(mvpLocalKeyFrames.size()>80)
            break;

        KeyFrame* pKF = *itKF;

        const vector<KeyFrame*> vNeighs = pKF->GetBestCovisibilityKeyFrames(10);

        for(vector<KeyFrame*>::const_iterator itNeighKF=vNeighs.begin(), itEndNeighKF=vNeighs.end(); itNeighKF!=itEndNeighKF; itNeighKF++)
        {
            KeyFrame* pNeighKF = *itNeighKF;
            if(!pNeighKF->isBad())
            {
                if(pNeighKF->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFrames.push_back(pNeighKF);
                    pNeighKF->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        const set<KeyFrame*> spChilds = pKF->GetChilds();
        for(set<KeyFrame*>::const_iterator sit=spChilds.begin(), send=spChilds.end(); sit!=send; sit++)
        {
            KeyFrame* pChildKF = *sit;
            if(!pChildKF->isBad())
            {
                if(pChildKF->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFrames.push_back(pChildKF);
                    pChildKF->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        KeyFrame* pParent = pKF->GetParent();
        if(pParent)
        {
            if(pParent->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
            {
                mvpLocalKeyFrames.push_back(pParent);
                pParent->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                break;
            }
        }

    }

    if(pKFmax)
    {
        mpReferenceKF = pKFmax;
        mCurrentFrame.mpReferenceKF = mpReferenceKF;
    }
}

bool Tracking::Relocalization()
{
    cout << endl;
    cout << "-R- START RELOC PROCESS - Frame " << mCurrentFrame.mnId << endl;
    nRelocalizationsCounter++;
    pStatsReloc->NewIteration(nRelocalizationsCounter);
    if ( ((nRelocalizationsCounter-1)==0) || (((nRelocalizationsCounter-1)%10)==0) )
        pStatsReloc->ColumnHeadersReloc();

    float reloct1 = 0.0;
    float reloct11, reloct12;
    double relocf = cv::getTickFrequency();

    if (bDebugMode) cout << "-R- Computing frame BoW and looking for Reloc candidates" << endl;

    // Compute Bag of Words Vector
    mCurrentFrame.ComputeBoW();

    // Track Lost: Query KF Database for KF candidates for relocalization when tracking is lost
    vector<KeyFrame*> vpCandidateKFs = mpKeyFrameDB->DetectRelocalizationCandidates(&mCurrentFrame);

    if(vpCandidateKFs.empty())
    {
        cout << "-R- DBoW can't find any candidate KFs" << endl;
        pStatsReloc->AddValue(0);
        pStatsRelocS1->AddValue(mCurrentFrame.mnId);
        pStatsRelocS1->AddValue(0);
        pStatsPR->AddValue(mCurrentFrame.mnId);
        pStatsPR->AddValue(0);
        pStatsPR->AddValue(0);
        pStatsPR->NewLine();
        return false;
    }

    const int nKFs = vpCandidateKFs.size();
    if (bDebugMode) cout << "-R- " << nKFs << " candidate KFs found" << endl;
    pStatsReloc->NewLine();
    pStatsReloc->AddValue(nKFs);
    pStatsRelocS1->AddValue(nKFs);

    // We perform first an ORB matching with each candidate
    // If enough matches are found we setup a PnP solver
    ORBmatcher matcher(0.75,true);

    vector<PnPsolver*> vpPnPsolvers;
    vpPnPsolvers.resize(nKFs);

    vector<vector<MapPoint*> > vvpMapPointMatches;
    vector<vector<cv::KeyPoint*> > vvpKeyPointMatches;
    vvpMapPointMatches.resize(nKFs);
    vvpKeyPointMatches.resize(nKFs);

    vector<bool> vbDiscarded;
    vbDiscarded.resize(nKFs);

    int nCandidates = 0;    //(>4)
    int nCandidates15 = 0;  //(>15)
    vector<int> number_of_matches;

    if (bDebugMode) cout << "-R- Looking for ORB matches using DBoW" << endl;
    for(int i=0; i<nKFs; i++)
    {
        KeyFrame* pKF = vpCandidateKFs[i];
        if(pKF->isBad())
            vbDiscarded[i] = true;
        else
        {
            int nmatches = matcher.SearchByBoW(pKF,mCurrentFrame,vvpMapPointMatches[i]);
            number_of_matches.push_back(nmatches);
            if (bDebugMode) cout << "    KF " << i+1 << "of" << vpCandidateKFs.size() << " (KFid " << pKF->mnId << "):  ";
            if(nmatches<4)  // prev 15
            {
                vbDiscarded[i] = true;
                if (bDebugMode) cout << "NOK (" << nmatches << " matches, not enough)" << endl;
                continue;
            }
            else
            {
                if (bDebugMode) cout << "OK (" << nmatches << " matches, PnPsolver added)" << endl;
                PnPsolver* pSolver = new PnPsolver(mCurrentFrame,vvpMapPointMatches[i]);
                pSolver->SetRansacParameters(0.99,10,300,4,0.5,5.991);
                vpPnPsolvers[i] = pSolver;
                nCandidates++;
                if (nmatches>=15)
                    nCandidates15++;
            }
        }
    }

    if (nCandidates==0)
        cout << "-R- Unable to find enough matches among the KFs." << endl;
    else {
        cout << "-R- " << nCandidates << " KFs with enough ORB matches. 5 RANSAC its with each KF." << endl;
        kpiP++;
    }

    pStatsReloc->AddValue(nCandidates15);
    pStatsReloc->AddValue(nCandidates);

    pStatsRelocS1->AddValue(mCurrentFrame.mnId);
    pStatsRelocS1->AddValue(nCandidates15);
    pStatsRelocS1->AddValue(nCandidates);
    for (unsigned int nn=0; nn<number_of_matches.size(); nn++)
        pStatsRelocS1->AddValue(number_of_matches[nn]);
    pStatsRelocS1->NewLine();


    // Alternatively perform some iterations of P4P RANSAC
    // Until we found a camera pose supported by enough inliers
    bool bMatch = false;
    ORBmatcher matcher2(0.9,true);

    int nCount = 0;
    int nPosesR = 0;
    int nPosesNR = 0;

    while(nCandidates>0 && !bMatch)
    {
        for(int i=0; i<nKFs; i++)
        {
            if(vbDiscarded[i])
                continue;

            nCount++;
            if (bDebugMode) cout << "    Candidate " << nCount << endl;

            // Perform 5 Ransac Iterations
            vector<bool> vbInliers;
            int nInliers;
            bool bNoMore;

            // Choose relocation method: 1(original) 2(TFG)
            int iterationModel = 2;

            PnPsolver* pSolver = vpPnPsolvers[i];
            cv::Mat Tcw;
            vector<MapPoint*> vpMatchesByProjection;
            vpMatchesByProjection = vector<MapPoint*>(vvpMapPointMatches[i].size(), static_cast<MapPoint*>(NULL));

            if (iterationModel==1)
            {
                reloct11 = cv::getTickCount();
                Tcw = pSolver->iterate(5,bNoMore,vbInliers,nInliers);
                reloct12 = cv::getTickCount();
                reloct1 = (reloct12-reloct11)/relocf;
                nPosesR++;
                pStatsReloc->AddValue(nInliers);
            }
            else if (iterationModel==2)
            {
                reloct11 = cv::getTickCount();
                Tcw = pSolver->iterate(mCurrentFrame,vpMatchesByProjection,5,bNoMore,vbInliers,nInliers,mpMap,pStatsRelocS2);
                reloct12 = cv::getTickCount();
                reloct1 = (reloct12-reloct11)/relocf;
                nPosesNR++;
                pStatsReloc->AddValue(nInliers);
            }
            if (bDebugMode) cout << "        PnP- Iteration time " << reloct1 << endl;
            pStatsReloc->AddValueFl(reloct1);

            // If Ransac reachs max. iterations discard keyframe
            if(bNoMore)
            {
                vbDiscarded[i]=true;
                nCandidates--;
            }

            // If a Camera Pose is computed, optimize
            if(!Tcw.empty())
            {
                cout << "        Camera pose computed, optimizing." << endl;

                cv::Mat mTcwBackup, mTcwR, mTcwNR;
                Tcw.copyTo(mTcwBackup);
                Tcw.copyTo(mCurrentFrame.mTcw);

                set<MapPoint*> sFound;
                /*const int np = vbInliers.size();*/
                int nAdded = 0;
                int nAddedProj = 0;

                for(unsigned int j=0; j<vbInliers.size(); j++)     // Set inliers as MPs in mCurrentFrame
                {
                    if(vbInliers[j])
                    {
                        mCurrentFrame.mvpMapPoints[j]=vvpMapPointMatches[i][j];
                        sFound.insert(vvpMapPointMatches[i][j]);
                        if (mCurrentFrame.mvpMapPoints[j])
                            nAdded++;
                    }
                    else
                        mCurrentFrame.mvpMapPoints[j]=NULL;
                }
                for (unsigned int j=0; j<vpMatchesByProjection.size(); j++)
                {
                    if (vpMatchesByProjection[j])
                    {
                        //nAddedProj++;
                        mCurrentFrame.mvpMapPoints[j] = vpMatchesByProjection[j];
                        sFound.insert(vpMatchesByProjection[j]);
                        if (mCurrentFrame.mvpMapPoints[j])
                            nAddedProj++;
                    }
                }

                int nTotalAdded = 0;
                for (unsigned int j=0; j<mCurrentFrame.mvpMapPoints.size(); j++)
                {
                    if (mCurrentFrame.mvpMapPoints[j])
                        nTotalAdded++;
                }

                if (bDebugMode) cout << "        " << nAdded << " inliers added." << endl;
                cout << "        Inliers added: DBoW(" << nAdded << ") Proj(" << nAddedProj << ")" << endl;
                pStatsReloc->AddValue(nAdded);
                int nGood = 0;

                float reloct2_R = 0.0;
                float reloct2_NR = 0.0;
                float reloct21, reloct22, reloct23;
                float reloct3_R = 0.0;
                float reloct3_NR = 0.0;
                float reloct31, reloct32, reloct33;
                float reloct4_R = 0.0;
                float reloct4_NR = 0.0;
                float reloct41, reloct42, reloct43;

                RestoreRigidityFlag();
                int nFinalInliers_R = 0;
                int nFinalInliers_NR = 0;

                reloct21 = cv::getTickCount();

                // First we attempt Rigid Relocation
                int nGoodR = Optimizer::PoseOptimization(&mCurrentFrame);
                mCurrentFrame.mTcw.copyTo(mTcwR);
                SetRigidityFlag(true);

                reloct22 = cv::getTickCount();
                reloct2_R = (reloct22-reloct21)/relocf;

                // Restore camPose from Backup. Attempt NonRigid Relocation. Save Non-Rigid pose optimization
                mTcwBackup.copyTo(mCurrentFrame.mTcw);
                int nGoodNR = Optimizer::PoseOptimizationNR(&mCurrentFrame,mpMap,mpFrameDrawer,mpMapDrawer,nElType,bDebugMode);
                mCurrentFrame.mTcw.copyTo(mTcwNR);
                SetRigidityFlag(false);
                mCurrentFrame.mTcw.copyTo(mTcwBackup);

                reloct23 = cv::getTickCount();
                reloct2_NR = (reloct23-reloct22)/relocf;
                if (bDebugMode) cout << "             S1- Time(R-NR)(" << reloct2_R << "-" << reloct2_NR << ") \t R(" << nGoodR << ")  NR(" << nGoodNR << ") \t ";
                if((nGoodR<10) && (nGoodNR<10))
                {
                    if (bDebugMode) cout << "Not enough." << endl;
                    //continue;
                }
                else if (nGoodR>=10 && nGoodNR<10)
                {
                    if (bDebugMode) cout << "Enough rigid." << endl;
                    mTcwR.copyTo(mCurrentFrame.mTcw);
                    nGood = nGoodR;
                }
                else if (nGoodNR>=10)
                {
                    if (bDebugMode) cout << "Enough non-rigid." << endl;
                    nGood = nGoodNR;
                }

                nFinalInliers_R = nGoodR;
                nFinalInliers_NR = nGoodNR;

                pStatsReloc->AddValue(nGoodR);
                pStatsReloc->AddValue(nGoodNR);
                pStatsReloc->AddValueFl(reloct2_NR);

                for(int io =0 ; io<mCurrentFrame.N; io++)
                    if(mCurrentFrame.mvbOutlier[io])
                        mCurrentFrame.mvpMapPoints[io]=static_cast<MapPoint*>(NULL);

                // If few inliers, search by projection in a coarse window and optimize again
                if(nGood<50)
                {
                    if (bDebugMode) cout << "                 Few inliers (" << nGood << "<50)." << endl;
                    if (bDebugMode) cout << "             S2- Search in a coarse window (projection)." << endl;
                    int nadditional = matcher2.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFound,10,100);
                    cout << "                 " << nadditional+nGood << " inliers pre-optimization." << endl;

                    if (nadditional+nGood<50)
                        cout << "                 not enough (" << nadditional+nGood << "<50)." << endl;

                    if(nadditional+nGood>=50)
                    {
                        RestoreRigidityFlag();

                        reloct31 = getTickCount();

                        nGoodR = Optimizer::PoseOptimization(&mCurrentFrame);
                        mCurrentFrame.mTcw.copyTo(mTcwR);
                        SetRigidityFlag(true);

                        reloct32 = getTickCount();
                        reloct3_R = (reloct32-reloct31)/relocf;

                        mTcwBackup.copyTo(mCurrentFrame.mTcw);
                        nGoodNR = Optimizer::PoseOptimizationNR(&mCurrentFrame,mpMap,mpFrameDrawer,mpMapDrawer,nElType,bDebugMode);
                        mCurrentFrame.mTcw.copyTo(mTcwNR);
                        SetRigidityFlag(false);
                        mCurrentFrame.mTcw.copyTo(mTcwBackup);

                        reloct33 = getTickCount();
                        reloct3_NR = (reloct33-reloct32)/relocf;
                        pStatsReloc->AddValue(nGoodR);
                        pStatsReloc->AddValue(nGoodNR);
                        pStatsReloc->AddValueFl(reloct3_NR);

                        if (bDebugMode) cout << "                 Time(R-NR)(" << reloct3_R << "-" << reloct3_NR << ") \t R(" << nGoodR << ")  NR(" << nGoodNR << ") \t ";
                        if((nGoodR<10) && (nGoodNR<10))
                        {
                            if (bDebugMode) cout << "Not enough." << endl;
                            //continue;
                        }
                        else if (nGoodR>=10 && nGoodNR<10)
                        {
                            if (bDebugMode) cout << "Enough rigid." << endl;
                            mTcwR.copyTo(mCurrentFrame.mTcw);
                            nGood = nGoodR;
                        }
                        else if (nGoodNR>=10)
                        {
                            if (bDebugMode) cout << "Enough non-rigid." << endl;
                            nGood = nGoodNR;
                        }

                        nFinalInliers_R = nGoodR;
                        nFinalInliers_NR = nGoodNR;

                        // If many inliers but still not enough, search by projection again in a narrower window
                        // the camera has been already optimized with many points
                        if(nGood>30 && nGood<50)
                        {
                            if (bDebugMode) cout << "                 Still not enough (30<" << nGood << "<50)." << endl;
                            if (bDebugMode) cout << "             S3- Search in a narrowed window (projection)." << endl;
                            sFound.clear();
                            for(int ip =0; ip<mCurrentFrame.N; ip++)
                                if(mCurrentFrame.mvpMapPoints[ip])
                                    sFound.insert(mCurrentFrame.mvpMapPoints[ip]);
                            nadditional =matcher2.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFound,3,64);
                            if (bDebugMode) cout << "        " << nadditional+nGood << " inliers pre-optimization." << endl;

                            if (nadditional+nGood<50)
                                if (bDebugMode) cout << "                 not enough (" << nadditional+nGood << "<50)." << endl;

                            // Final optimization
                            if(nGood+nadditional>=50)
                            {
                                RestoreRigidityFlag();

                                reloct41 = getTickCount();

                                nGoodR = Optimizer::PoseOptimization(&mCurrentFrame);
                                SetRigidityFlag(true);

                                reloct42 = getTickCount();
                                reloct4_R = (reloct42-reloct41)/relocf;

                                mTcwBackup.copyTo(mCurrentFrame.mTcw);
                                nGoodNR = Optimizer::PoseOptimizationNR(&mCurrentFrame,mpMap,mpFrameDrawer,mpMapDrawer,nElType,bDebugMode);
                                mCurrentFrame.mTcw.copyTo(mTcwNR);
                                SetRigidityFlag(false);

                                reloct43 = getTickCount();
                                reloct4_NR = (reloct43-reloct42)/relocf;
                                if (bDebugMode) cout << "                 Time(R-NR)(" << reloct4_R << "-" << reloct4_NR << ") \t R(" << nGoodR << ")  NR(" << nGoodNR << ") \t ";
                                if((nGoodR<50) && (nGoodNR<50))
                                {
                                    if (bDebugMode) cout << "Not enough." << endl;
                                    //continue;
                                }
                                else if (nGoodR>=50 && nGoodNR<50)
                                {
                                    if (bDebugMode) cout << "Enough rigid." << endl;
                                    mTcwR.copyTo(mCurrentFrame.mTcw);
                                    nGood = nGoodR;
                                }
                                else if (nGoodNR>=50)
                                {
                                    if (bDebugMode) cout << "Enough non-rigid." << endl;
                                    nGood = nGoodNR;
                                }

                                nFinalInliers_R = nGoodR;
                                nFinalInliers_NR = nGoodNR;

                                pStatsReloc->AddValue(nGoodR);
                                pStatsReloc->AddValue(nGoodNR);
                                pStatsReloc->AddValueFl(reloct4_NR);

                                for(int io =0; io<mCurrentFrame.N; io++)
                                    if(mCurrentFrame.mvbOutlier[io])
                                        mCurrentFrame.mvpMapPoints[io]=NULL;
                            }
                        }
                    }
                }

                if (nGood>=50)
                    pStatsRelocS3->AddValue(1);
                else
                    pStatsRelocS3->AddValue(0);

                pStatsRelocS3->AddValue(nTotalAdded);
                pStatsRelocS3->AddValue(nFinalInliers_R);
                pStatsRelocS3->AddValue(nFinalInliers_NR);

                pStatsRelocS3->NewLine();

                // If the pose is supported by enough inliers stop ransacs and continue
                if(nGood>=50)
                {
                    bMatch = true;
                    cout << "                 Enough inliers (" << nGood << ")." << endl;
                    pStatsReloc->AddValue(nGood);
                    break;
                }
                else
                    nCandidates--;
            }

            pStatsReloc->EmptyCols(3);
        }
    }

    if(!bMatch)
    {
        cout << "-R- Relocation failure!" << endl;
        pStatsReloc->AddValue(111);
        pStatsReloc->NewLine();
        pStatsPR->AddValue(mCurrentFrame.mnId);
        pStatsPR->AddValue(0);
        pStatsPR->AddValue(0);
        pStatsPR->NewLine();
        return false;
    }
    else
    {
        mnLastRelocFrameId = mCurrentFrame.mnId;
        cout << "-R- Relocation successful!" << endl;
        pStatsReloc->AddValue(222);
        pStatsReloc->NewLine();
        //pStatsPR->AddValue(mCurrentFrame.mnId);
        //pStatsPR->AddValue(1);
        return true;
        //if (bTestAllFrames) return false;   // Don't allow relocalization for statistic data gathering.
        //else return true;    // Normal mode
    }

}

void Tracking::Reset()
{

    cout << "System Reseting" << endl;
    if(mpViewer)
    {
        mpViewer->RequestStop();
        while(!mpViewer->isStopped())
            usleep(3000);
    }

    // Reset Local Mapping
    cout << "Reseting Local Mapper...";
    mpLocalMapper->RequestReset();
    cout << " done" << endl;

    // Reset Loop Closing
    cout << "Reseting Loop Closing...";
    mpLoopClosing->RequestReset();
    cout << " done" << endl;

    // Clear BoW Database
    cout << "Reseting Database...";
    mpKeyFrameDB->clear();
    cout << " done" << endl;

    // Clear Map (this erase MapPoints and KeyFrames)
    mpMap->clear();

    KeyFrame::nNextId = 0;
    Frame::nNextId = 0;
    mState = NO_IMAGES_YET;

    if(mpInitializer)
    {
        delete mpInitializer;
        mpInitializer = static_cast<Initializer*>(NULL);
    }

    mlRelativeFramePoses.clear();
    mlpReferences.clear();
    mlFrameTimes.clear();
    mlbLost.clear();

    if(mpViewer)
        mpViewer->Release();
}

void Tracking::ChangeCalibration(const string &strSettingPath)
{
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    float fx = fSettings["Camera.fx"];
    float fy = fSettings["Camera.fy"];
    float cx = fSettings["Camera.cx"];
    float cy = fSettings["Camera.cy"];

    cv::Mat K = cv::Mat::eye(3,3,CV_32F);
    K.at<float>(0,0) = fx;
    K.at<float>(1,1) = fy;
    K.at<float>(0,2) = cx;
    K.at<float>(1,2) = cy;
    K.copyTo(mK);

    cv::Mat DistCoef(4,1,CV_32F);
    DistCoef.at<float>(0) = fSettings["Camera.k1"];
    DistCoef.at<float>(1) = fSettings["Camera.k2"];
    DistCoef.at<float>(2) = fSettings["Camera.p1"];
    DistCoef.at<float>(3) = fSettings["Camera.p2"];
    const float k3 = fSettings["Camera.k3"];
    if(k3!=0)
    {
        DistCoef.resize(5);
        DistCoef.at<float>(4) = k3;
    }
    DistCoef.copyTo(mDistCoef);

    mbf = fSettings["Camera.bf"];

    Frame::mbInitialComputations = true;
}

void Tracking::InformOnlyTracking(const bool &flag)
{
    mbOnlyTracking = flag;
}

float Tracking::RadiousByViewingCosT(const float &viewCos)
{
    if (viewCos>0.998)
        return 3.0;
    else
        return 4.5;
}

cv::Mat Tracking::GetDistCoef()
{
    return mDistCoef;
}

void Tracking::SetRigidityFlag(bool rigid)
{
    for (unsigned ii=0; ii<mCurrentFrame.mvpMapPoints.size(); ii++)
    {
        MapPoint* pMPf = mCurrentFrame.mvpMapPoints[ii];
        if (pMPf)
        {
            if (rigid==true)
                pMPf->bFlagRigid = true;
            else
                pMPf->bFlagNonRigid = true;
        }
    }
}

void Tracking::RestoreRigidityFlag()
{
    for (unsigned int ii=0; ii<mCurrentFrame.mvpMapPoints.size(); ii++)
    {
        MapPoint* pMPf = mCurrentFrame.mvpMapPoints[ii];
        if (pMPf)
        {
            pMPf->bFlagRigid = false;
            pMPf->bFlagNonRigid = false;
        }
    }
}


bool Tracking::LoadMap(vector<KeyFrame*> &output_KFs, vector<MapPoint*> &output_MPs)
{
    cout << endl << endl << "Lading map:" << endl;

    //const char * filename_F_char = "../../../evaluation/map/Frame.bin";
    const char * filename_KFs_char = "../../../evaluation/map/KeyFrames.bin";
    const char * filename_MPs_char = "../../../evaluation/map/MapPoints.bin";

    //FILE *pFile_F;
    FILE *pFile_KFs;
    FILE *pFile_MPs;

    /*pFile_F = fopen ( filename_F_char , "r+b" );
    if (pFile_F==NULL)
    {
        cout << "Unable to open Frame file" << endl;
        return false;
    }
    fseek(pFile_F,0,SEEK_END);
    size_t lSize_F = ftell(pFile_F);
    fclose(pFile_F);*/

    pFile_KFs = fopen ( filename_KFs_char , "r+b" );
    if (pFile_KFs==NULL)
    {
        cout << "Unable to open KeyFrames file" << endl;
        return false;
    }
    fseek(pFile_KFs,0,SEEK_END);
    size_t lSize_KFs = ftell(pFile_KFs);
    fclose(pFile_KFs);

    pFile_MPs = fopen ( filename_MPs_char , "r+b" );
    if (pFile_MPs==NULL)
    {
        cout << "Unable to open MapPoints file" << endl;
        return false;
    }
    fseek(pFile_MPs,0,SEEK_END);
    size_t lSize_MPs = ftell(pFile_MPs);
    fclose(pFile_MPs);

    //long size_F = sizeof(Frame);
    long size_KF = sizeof(KeyFrame);
    long size_MP = sizeof(MapPoint);

    //unsigned int nFs = lSize_F/size_F;
    unsigned int nKFs = lSize_KFs/size_KF;
    unsigned int nMPs = lSize_MPs/size_MP;




    /*pFile_F = fopen ( filename_F_char , "r+b" );
    cout << "  " << nFs << " Frames to load ... ";

    Frame buffer_F = Frame();
    int ressult_F = fread(&buffer_F,sizeof(Frame),1,pFile_F);

    fclose(pFile_F);

    cout << "OK(" << ressult_F << ")" << endl;*/



    pFile_KFs = fopen ( filename_KFs_char , "r+b" );
    cout << "  " << nKFs << " KeyFrames to load ... ";

    vector<KeyFrame*> vpKFs;
    int nOkKFs = 0;

    for (unsigned int i=0; i<nKFs; i++)
    {
        KeyFrame* pKF;// = new KeyFrame();
        pKF = (KeyFrame*) malloc(sizeof(KeyFrame));
        fseek(pFile_KFs,size_KF*i,0);
        int ressult_KF = fread(pKF,size_KF,1,pFile_KFs);
        vpKFs.push_back(pKF);
        nOkKFs = nOkKFs + ressult_KF;
    }

    fclose(pFile_KFs);
    cout << " OK(" << nOkKFs << ")" << endl;



    pFile_MPs = fopen ( filename_MPs_char , "r+b" );
    cout << "  " << nMPs << " MapPoints to load ... ";

    vector<MapPoint*> vpMPs;
    int nOkMPs = 0;

    for (unsigned int i=0; i<nMPs; i++)
    {
        MapPoint* pMP = new MapPoint();
        pMP = (MapPoint*) malloc(sizeof(MapPoint));
        fseek(pFile_MPs,size_MP*i,0);
        int ressult_MP = fread(pMP,size_MP,1,pFile_MPs);
        vpMPs.push_back(pMP);
        nOkMPs = nOkMPs + ressult_MP;
    }

    fclose(pFile_MPs);
    cout << " OK(" << nOkMPs << ")" << endl << endl;

    output_KFs = vpKFs;
    output_MPs = vpMPs;

    return true;
}


bool Tracking::BuildLoadedMap()
{
    int countKFs = 0;
    int countMPs = 0;

    cout << "Building Map..." << endl;



    unsigned int nMaxKfId = 0;
    vector<int> vnKfIds;

    for (unsigned int i=0; i<vpKFs_Loaded.size(); i++)
    {
        KeyFrame *pKF = vpKFs_Loaded[i];
        if (!pKF)
            continue;

        pKF->UpdateData(mpMap, mpKeyFrameDB, mpORBVocabulary);

        vnKfIds.push_back(pKF->mnId);
        if (pKF->mnId > nMaxKfId)
            nMaxKfId = pKF->mnId;

        mpMap->AddKeyFrame(pKF);
        countKFs++;
    }

    // Inverse index for KFs, position of each kf by its Id
    vector<int> vnKfIndex = vector<int>(nMaxKfId,0);
    for (unsigned int i=0; i<vnKfIds.size(); i++)
    {
        int VecPos = vnKfIds[i];
        vnKfIndex[VecPos] = i;
    }

    cout << "  KeyFrames Updated" << endl;



    vector<int> vnMP_Id;
    vector<int> vnMP_KFId;
    unsigned int nMaxMpId = 0;
    vector<vector<int> > vvnObservers;

    for (unsigned int i=0; i<vpMPs_Loaded.size(); i++)
    {
        MapPoint *pMP = vpMPs_Loaded[i];
        if (!pMP)
            continue;

        vector<int> vnObservers;

        pMP->UpdateData(mpMap, vnMP_Id, vnMP_KFId, vnObservers);

        vvnObservers.push_back(vnObservers);
        if (pMP->mnId > nMaxMpId)
            nMaxMpId = pMP->mnId;

        mpMap->AddMapPoint(vpMPs_Loaded[i]);
        countMPs++;
    }

    // Inverse index for MPs, position of each mp by its Id
    vector<int> vnMpIndex = vector<int>(nMaxMpId,0);
    for (unsigned int i=0; i<vnMP_Id.size(); i++)
    {
        int VecPos = vnMP_Id[i];
        vnMpIndex[VecPos] = i;
    }

    cout << "  MapPoints Updated" << endl;



    {
        unique_lock<mutex> lock(mMutexLoadMap);

        for (unsigned int i=0; i<vpMPs_Loaded.size(); i++)
        {
            //cout << "000" << endl;
            MapPoint *pMP = vpMPs_Loaded[i];
            if (!pMP)
                continue;
            //pMP->ClearObservations();
            //cout << pMP->nObsCounter << endl;
            for (unsigned int j=0; j<pMP->nObsCounter; j++)
            {
                //unsigned int todo = pMP->vvvnObservations[j][0];
                unsigned int nKfId = pMP->vvvnObservations[j][1];
                //unsigned int nMpIdx = pMP->vvvnObservations[j][2];

                if (nKfId>=vnKfIndex.size())
                    continue;

                unsigned int nKfIdx = vnKfIndex[nKfId];

                if (nKfIdx==0)
                    continue;

                if (nKfIdx>=vpKFs_Loaded.size())
                    continue;

                size_t size1= sizeof(*vpKFs_Loaded[nKfIdx]);
                if (size1 != sizeof(KeyFrame))
                    continue;

                KeyFrame *pKF = vpKFs_Loaded[nKfIdx];

                if (!pKF)
                    continue;

                /*if (todo == 1)
                {
                    cout << "pre" << endl;
                    pMP->AddObservationLoadMap(pKF,nMpIdx);
                    cout << "post" << endl;
                }
                else if (todo == 2)
                    pMP->EraseObservationLoadMap(pKF,nMpIdx);
                else
                    continue;*/
            }
        }
    }


    cout << "vpMPs_Loaded" << endl;

    for (unsigned int i=0; i<vpKFs_Loaded.size(); i++)
    {
        KeyFrame *pKF = vpKFs_Loaded[i];
        if (!pKF)
            continue;
        cout << pKF->mnId << endl;


        for (unsigned int j=0; j<pKF->nMPsCounter-1; j++)
        {
            cout << " ant " << endl;
            int todo = pKF->vvvnMPsObserved[j][0];
            int nMpId = pKF->vvvnMPsObserved[j][1];
            //int nMpIdx = pKF->vvvnMPsObserved[j][2];
            //cout << " post " << todo << " " << nMpIdx << endl;



            cout << "ant2" << endl;
            int nMpIdx = vnMpIndex[nMpId];
            if (nMpIdx==0)
                continue;
            MapPoint *pMP = vpMPs_Loaded[nMpIdx];
            cout << "post2" << endl;

            if (!pMP)
                continue;

            cout << "MP id " << pMP->mnId << endl;



            if (todo == 1)
            {
                pKF->AddMapPointLoadMap(pMP,nMpIdx);
                cout << "a" << endl;
            }
            else if (todo == 2)
            {
                pKF->EraseMapPointMatchLoadMap(nMpIdx);
                cout << "b" << endl;
            }
            else
            {
                cout << "c" << endl;
                continue;
            }
        }
    }
    cout << "vpKFs_Loaded" << endl;

    for (unsigned int i=0; i<vpMPs_Loaded.size(); i++)
    {
        MapPoint *pMP = vpMPs_Loaded[i];

        pMP->ComputeDistinctiveDescriptors();

        pMP->UpdateNormalAndDepth();

        mpMap->AddMapPoint(pMP);
    }
    cout << "vpMPs_Loaded" << endl;






    /*for (unsigned int i=0; i<vvnObservers.size(); i++)
    {
        MapPoint *pMP = vpMPs_Loaded[i];

        int nKfRefId = vnMP_KFId[i];
        int nKfRefIdx = vnKfIndex[nKfRefId];

        KeyFrame *pKF = vpKFs_Loaded[nKfRefIdx];

        //pMP->AddObservation(pKF,idx1);
        //pKF->AddMapPoint(pMP,idx1);

        for (unsigned int j=0; j<vvnObservers[i].size(); j++)
        {
            int nKfId = vvnObservers[i][j];
            if (nKfId == nKfRefId)
                continue;

            int nKfIdx = vnKfIndex[nKfId];

            KeyFrame *pKF2 = vpKFs_Loaded[nKfIdx];

            //pMP->AddObservation(pKF2,idx2);
            //pKF2->AddMapPoint(pMP,idx2);
        }

        pMP->ComputeDistinctiveDescriptors();

        pMP->UpdateNormalAndDepth();

        mpMap->AddMapPoint(pMP);
    }*/

    cout << "  Relations updated" << endl;

    if (countKFs>0 || countMPs>0)
        return true;
    else
        return false;
}






} //namespace ORB_SLAM
