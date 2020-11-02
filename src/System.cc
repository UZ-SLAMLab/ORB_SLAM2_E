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



#include "System.h"
#include "Converter.h"
#include <thread>
#include <pangolin/pangolin.h>
#include <iomanip>
//#include<ros/ros.h>
//#include<ros/package.h>

#include "Statistics.h"

#include "fstream"

namespace ORB_SLAM2
{

System::System(const string &strVocFile, const string &strSettingsFile, const eSensor sensor,
               const bool bUseViewer):mSensor(sensor), mpViewer(static_cast<Viewer*>(NULL)), mbReset(false),mbActivateLocalizationMode(false),
        mbDeactivateLocalizationMode(false)
{
    // Output welcome message
    cout << endl <<
    "ORB-SLAM2 Copyright (C) 2014-2016 Raul Mur-Artal, University of Zaragoza." << endl <<
    "This program comes with ABSOLUTELY NO WARRANTY;" << endl  <<
    "This is free software, and you are welcome to redistribute it" << endl <<
    "under certain conditions. See LICENSE.txt." << endl << endl;

    cout << "Input sensor was set to: ";

    if(mSensor==MONOCULAR)
        cout << "Monocular" << endl;
    else if(mSensor==STEREO)
        cout << "Stereo" << endl;
    else if(mSensor==RGBD)
        cout << "RGB-D" << endl;

    //Check settings file
    cv::FileStorage fsSettings(strSettingsFile.c_str(), cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
       cerr << "Failed to open settings file at: " << strSettingsFile << endl;
       exit(-1);
    }

    //Load ORB Vocabulary
    cout << endl << "Loading ORB Vocabulary. This could take a while..." << endl << endl;

    mpVocabulary = new ORBVocabulary();
    bool bVocLoad = mpVocabulary->loadFromTextFile(strVocFile);
    if(!bVocLoad)
    {
        cerr << "Wrong path to vocabulary. " << endl;
        cerr << "Falied to open at: " << strVocFile << endl;
        exit(-1);
    }
    //cout << "Vocabulary loaded!" << endl << endl;

    //Create KeyFrame Database
    mpKeyFrameDatabase = new KeyFrameDatabase(*mpVocabulary);

    //Create the Map
    mpMap = new Map();

    //cout << "Map ready" << endl;

    //Create Drawers. These are used by the Viewer
    mpFrameDrawer = new FrameDrawer(mpMap);
    mpMapDrawer = new MapDrawer(mpMap, strSettingsFile);

    //cout << "Drawers ready" << endl;

    //Initialize the Tracking thread
    //(it will live in the main thread of execution, the one that called this constructor)
    mpTracker = new Tracking(this, mpVocabulary, mpFrameDrawer, mpMapDrawer,
                             mpMap, mpKeyFrameDatabase, strSettingsFile, mSensor);

    //cout << "Tracker ready" << endl;

    //Initialize the Local Mapping thread and launch
    mpLocalMapper = new LocalMapping(mpMap, mSensor==MONOCULAR);
    mptLocalMapping = new thread(&ORB_SLAM2::LocalMapping::Run,mpLocalMapper);

    //cout << "Local mapper ready" << endl;

    //Initialize the Loop Closing thread and launch
    mpLoopCloser = new LoopClosing(mpMap, mpKeyFrameDatabase, mpVocabulary, mSensor!=MONOCULAR);
    mptLoopClosing = new thread(&ORB_SLAM2::LoopClosing::Run, mpLoopCloser);

    //cout << "Loop closer ready" << endl;

    //Initialize the Viewer thread and launch
    if(bUseViewer)
    {
        mpViewer = new Viewer(this, mpFrameDrawer,mpMapDrawer,mpTracker,strSettingsFile);
        mptViewer = new thread(&Viewer::Run, mpViewer);
        mpTracker->SetViewer(mpViewer);
        //cout << "Viewer ready" << endl;
    }

    //Set pointers between threads
    mpTracker->SetLocalMapper(mpLocalMapper);
    mpTracker->SetLoopClosing(mpLoopCloser);

    mpLocalMapper->SetTracker(mpTracker);
    mpLocalMapper->SetLoopCloser(mpLoopCloser);

    mpLoopCloser->SetTracker(mpTracker);
    mpLoopCloser->SetLocalMapper(mpLocalMapper);

    //cout << "Pointers between threads ready" << endl << endl;
}

cv::Mat System::TrackStereo(const cv::Mat &imLeft, const cv::Mat &imRight, const double &timestamp)
{
    if(mSensor!=STEREO)
    {
        cerr << "ERROR: you called TrackStereo but input sensor was not set to STEREO." << endl;
        exit(-1);
    }

    // Check mode change
    {
        unique_lock<mutex> lock(mMutexMode);
        if(mbActivateLocalizationMode)
        {
            mpLocalMapper->RequestStop();

            // Wait until Local Mapping has effectively stopped
            while(!mpLocalMapper->isStopped())
            {
                usleep(1000);
            }

            mpTracker->InformOnlyTracking(true);
            mbActivateLocalizationMode = false;
        }
        if(mbDeactivateLocalizationMode)
        {
            mpTracker->InformOnlyTracking(false);
            mpLocalMapper->Release();
            mbDeactivateLocalizationMode = false;
        }
    }

    // Check reset
    {
    unique_lock<mutex> lock(mMutexReset);
    if(mbReset)
    {
        mpTracker->Reset();
        mbReset = false;
    }
    }

    cv::Mat Tcw = mpTracker->GrabImageStereo(imLeft,imRight,timestamp);

    unique_lock<mutex> lock2(mMutexState);
    mTrackingState = mpTracker->mState;
    mTrackedMapPoints = mpTracker->mCurrentFrame.mvpMapPoints;
    mTrackedKeyPointsUn = mpTracker->mCurrentFrame.mvKeysUn;
    return Tcw;
}

cv::Mat System::TrackRGBD(const cv::Mat &im, const cv::Mat &depthmap, const double &timestamp)
{
    if(mSensor!=RGBD)
    {
        cerr << "ERROR: you called TrackRGBD but input sensor was not set to RGBD." << endl;
        exit(-1);
    }

    // Check mode change
    {
        unique_lock<mutex> lock(mMutexMode);
        if(mbActivateLocalizationMode)
        {
            mpLocalMapper->RequestStop();

            // Wait until Local Mapping has effectively stopped
            while(!mpLocalMapper->isStopped())
            {
                usleep(1000);
            }

            mpTracker->InformOnlyTracking(true);
            mbActivateLocalizationMode = false;
        }
        if(mbDeactivateLocalizationMode)
        {
            mpTracker->InformOnlyTracking(false);
            mpLocalMapper->Release();
            mbDeactivateLocalizationMode = false;
        }
    }

    // Check reset
    {
    unique_lock<mutex> lock(mMutexReset);
    if(mbReset)
    {
        mpTracker->Reset();
        mbReset = false;
    }
    }

    cv::Mat Tcw = mpTracker->GrabImageRGBD(im,depthmap,timestamp);

    unique_lock<mutex> lock2(mMutexState);
    mTrackingState = mpTracker->mState;
    mTrackedMapPoints = mpTracker->mCurrentFrame.mvpMapPoints;
    mTrackedKeyPointsUn = mpTracker->mCurrentFrame.mvKeysUn;
    return Tcw;
}

cv::Mat System::TrackMonocular(const cv::Mat &im, const double &timestamp)
{
    if(mSensor!=MONOCULAR)
    {
        cerr << "ERROR: you called TrackMonocular but input sensor was not set to Monocular." << endl;
        exit(-1);
    }

    // Check mode change
    {
        unique_lock<mutex> lock(mMutexMode);
        if(mbActivateLocalizationMode)
        {
            mpLocalMapper->RequestStop();

            // Wait until Local Mapping has effectively stopped
            while(!mpLocalMapper->isStopped())
            {
                usleep(1000);
            }

            mpTracker->InformOnlyTracking(true);
            mbActivateLocalizationMode = false;
        }
        if(mbDeactivateLocalizationMode)
        {
            mpTracker->InformOnlyTracking(false);
            mpLocalMapper->Release();
            mbDeactivateLocalizationMode = false;
        }
    }

    // Check reset
    {
    unique_lock<mutex> lock(mMutexReset);
    if(mbReset)
    {
        mpTracker->Reset();
        mbReset = false;
    }
    }

    cv::Mat Tcw = mpTracker->GrabImageMonocular(im,timestamp);

    unique_lock<mutex> lock2(mMutexState);
    mTrackingState = mpTracker->mState;
    mTrackedMapPoints = mpTracker->mCurrentFrame.mvpMapPoints;
    mTrackedKeyPointsUn = mpTracker->mCurrentFrame.mvKeysUn;

    return Tcw;
}

void System::ActivateLocalizationMode()
{
    unique_lock<mutex> lock(mMutexMode);
    mbActivateLocalizationMode = true;
}

void System::DeactivateLocalizationMode()
{
    unique_lock<mutex> lock(mMutexMode);
    mbDeactivateLocalizationMode = true;
}

bool System::MapChanged()
{
    static int n=0;
    int curn = mpMap->GetLastBigChangeIdx();
    if(n<curn)
    {
        n=curn;
        return true;
    }
    else
        return false;
}

void System::Reset()
{
    unique_lock<mutex> lock(mMutexReset);
    mbReset = true;
}

void System::Shutdown()
{
    mpLocalMapper->RequestFinish();
    mpLoopCloser->RequestFinish();
    if(mpViewer)
        mpViewer->RequestFinish();

    // Wait until all thread have effectively stopped
    while(!mpLocalMapper->isFinished() || !mpLoopCloser->isFinished()  ||
          !mpViewer->isFinished()      || mpLoopCloser->isRunningGBA())
    {
        usleep(5000);
    }

    pangolin::BindToContext("ORB-SLAM2: Map Viewer");
}

void System::SaveTrajectoryTUM(const string &filename)
{
    cout << endl << "Saving camera trajectory to " << filename << " ..." << endl;
    if(mSensor==MONOCULAR)
    {
        cerr << "ERROR: SaveTrajectoryTUM cannot be used for monocular." << endl;
        return;
    }

    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    sort(vpKFs.begin(),vpKFs.end(),KeyFrame::lId);

    // Transform all keyframes so that the first keyframe is at the origin.
    // After a loop closure the first keyframe might not be at the origin.
    cv::Mat Two = vpKFs[0]->GetPoseInverse();

    ofstream f;
    f.open(filename.c_str());
    f << fixed;

    // Frame pose is stored relative to its reference keyframe (which is optimized by BA and pose graph).
    // We need to get first the keyframe pose and then concatenate the relative transformation.
    // Frames not localized (tracking failure) are not saved.

    // For each frame we have a reference keyframe (lRit), the timestamp (lT) and a flag
    // which is true when tracking failed (lbL).
    list<ORB_SLAM2::KeyFrame*>::iterator lRit = mpTracker->mlpReferences.begin();
    list<double>::iterator lT = mpTracker->mlFrameTimes.begin();
    list<bool>::iterator lbL = mpTracker->mlbLost.begin();
    for(list<cv::Mat>::iterator lit=mpTracker->mlRelativeFramePoses.begin(),
        lend=mpTracker->mlRelativeFramePoses.end();lit!=lend;lit++, lRit++, lT++, lbL++)
    {
        if(*lbL)
            continue;

        KeyFrame* pKF = *lRit;

        cv::Mat Trw = cv::Mat::eye(4,4,CV_32F);

        // If the reference keyframe was culled, traverse the spanning tree to get a suitable keyframe.
        while(pKF->isBad())
        {
            Trw = Trw*pKF->mTcp;
            pKF = pKF->GetParent();
        }

        Trw = Trw*pKF->GetPose()*Two;

        cv::Mat Tcw = (*lit)*Trw;
        cv::Mat Rwc = Tcw.rowRange(0,3).colRange(0,3).t();
        cv::Mat twc = -Rwc*Tcw.rowRange(0,3).col(3);

        vector<float> q = Converter::toQuaternion(Rwc);

        f << setprecision(6) << *lT << " " <<  setprecision(9) << twc.at<float>(0) << " " << twc.at<float>(1) << " " << twc.at<float>(2) << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
    }
    f.close();
    cout << endl << "trajectory saved!" << endl;
}


void System::SaveKeyFrameTrajectoryTUM(const string &filename)
{
    cout << endl << "Saving keyframe trajectory to " << filename << " ..." << endl;

    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    sort(vpKFs.begin(),vpKFs.end(),KeyFrame::lId);

    // Transform all keyframes so that the first keyframe is at the origin.
    // After a loop closure the first keyframe might not be at the origin.
    //cv::Mat Two = vpKFs[0]->GetPoseInverse();

    ofstream f;
    f.open(filename.c_str());
    f << fixed;

    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];

       // pKF->SetPose(pKF->GetPose()*Two);

        if(pKF->isBad())
            continue;

        cv::Mat R = pKF->GetRotation().t();
        vector<float> q = Converter::toQuaternion(R);
        cv::Mat t = pKF->GetCameraCenter();
        f << setprecision(6) << pKF->mTimeStamp << setprecision(7) << " " << t.at<float>(0) << " " << t.at<float>(1) << " " << t.at<float>(2)
          << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;

    }

    f.close();
    cout << endl << "trajectory saved!" << endl;
}

void System::SaveTrajectoryKITTI(const string &filename)
{
    cout << endl << "Saving camera trajectory to " << filename << " ..." << endl;
    if(mSensor==MONOCULAR)
    {
        cerr << "ERROR: SaveTrajectoryKITTI cannot be used for monocular." << endl;
        return;
    }

    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    sort(vpKFs.begin(),vpKFs.end(),KeyFrame::lId);

    // Transform all keyframes so that the first keyframe is at the origin.
    // After a loop closure the first keyframe might not be at the origin.
    cv::Mat Two = vpKFs[0]->GetPoseInverse();

    ofstream f;
    f.open(filename.c_str());
    f << fixed;

    // Frame pose is stored relative to its reference keyframe (which is optimized by BA and pose graph).
    // We need to get first the keyframe pose and then concatenate the relative transformation.
    // Frames not localized (tracking failure) are not saved.

    // For each frame we have a reference keyframe (lRit), the timestamp (lT) and a flag
    // which is true when tracking failed (lbL).
    list<ORB_SLAM2::KeyFrame*>::iterator lRit = mpTracker->mlpReferences.begin();
    list<double>::iterator lT = mpTracker->mlFrameTimes.begin();
    for(list<cv::Mat>::iterator lit=mpTracker->mlRelativeFramePoses.begin(), lend=mpTracker->mlRelativeFramePoses.end();lit!=lend;lit++, lRit++, lT++)
    {
        ORB_SLAM2::KeyFrame* pKF = *lRit;

        cv::Mat Trw = cv::Mat::eye(4,4,CV_32F);

        while(pKF->isBad())
        {
          //  cout << "bad parent" << endl;
            Trw = Trw*pKF->mTcp;
            pKF = pKF->GetParent();
        }

        Trw = Trw*pKF->GetPose()*Two;

        cv::Mat Tcw = (*lit)*Trw;
        cv::Mat Rwc = Tcw.rowRange(0,3).colRange(0,3).t();
        cv::Mat twc = -Rwc*Tcw.rowRange(0,3).col(3);

        f << setprecision(9) << Rwc.at<float>(0,0) << " " << Rwc.at<float>(0,1)  << " " << Rwc.at<float>(0,2) << " "  << twc.at<float>(0) << " " <<
             Rwc.at<float>(1,0) << " " << Rwc.at<float>(1,1)  << " " << Rwc.at<float>(1,2) << " "  << twc.at<float>(1) << " " <<
             Rwc.at<float>(2,0) << " " << Rwc.at<float>(2,1)  << " " << Rwc.at<float>(2,2) << " "  << twc.at<float>(2) << endl;
    }
    f.close();
    cout << endl << "trajectory saved!" << endl;
}

void System::SaveMap(const string &filename_MPs, const string &filename_KFs)
{
    double freq = cv::getTickFrequency();


    cout << endl << "Saving Map: " << endl;


    cout << "  KeyFrames ... ";

    const char * filename_KFs_char = filename_KFs.c_str();

    float t11 = cv::getTickCount();

    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    if (vpKFs.size()==0)
        return;

    FILE * pFile_KFs;
    pFile_KFs = fopen(filename_KFs_char,"w+b");
    for (unsigned int i=0; i<vpKFs.size(); i++)
    {
        //vpKFs[i]->ReadyToSave();
        fwrite(vpKFs[i],sizeof(KeyFrame),1,pFile_KFs);
    }
    fclose(pFile_KFs);

    float t12 = cv::getTickCount();
    float t1 = (t12-t11)/freq;

    cout << "OK (" << vpKFs.size() << " KeyFrames, " << sizeof(KeyFrame) << " BytesUd, " << t1 << " seconds)" << endl;


    cout << "  MapPoints ... ";

    const char * filename_MPs_char = filename_MPs.c_str();

    float t21 = cv::getTickCount();

    vector<MapPoint*> vpMPs = mpMap->GetAllMapPoints();

    FILE *pFile_MPs;
    pFile_MPs = fopen(filename_MPs_char,"w+b");
    for (unsigned int i=0; i<vpMPs.size(); i++)
    {
        vpMPs[i]->RecordObservers();
        fwrite(vpMPs[i],sizeof(MapPoint),1,pFile_MPs);
    }
    fclose(pFile_MPs);

    float t22 = cv::getTickCount();
    float t2 = (t22-t21)/freq;

    cout << "OK (" << vpMPs.size() << " MapPoints, " << sizeof(MapPoint) << " BytesUd, " << t2 << " seconds)" << endl;


    cout << "  Save MapPoint 3D coordinates for mesh generation ... ";

    /*const char * filename_MPs_3D_char = "../../../evaluation/map/MPs_3D.txt";

    vector<MapPoint*> vpMPs_3D = mpMap->GetAllMapPoints();

    FILE *pFile_MPs_3D;
    pFile_MPs_3D = fopen(filename_MPs_3D_char,"w+b");
    for (unsigned int i=0; i<vpMPs.size(); i++)
    {
        cv::Mat pos = vpMPs[i]->GetWorldPos();
        char comma = ";";

        for (unsigned int j=0; j<3; j++)
        {
            float pos.at<char>(49,39);
            fwrite(pos,sizeof(float),1,pFile_MPs_3D);
            fwrite(comma,sizeof(char),1,pFile_MPs_3D);
    }
    fclose(pFile_MPs_3D);*/

    vector<MapPoint*> vpMPs_3D = mpMap->GetAllMapPoints();

    ofstream MPs_3Ds;
    MPs_3Ds.open("../../../evaluation/map/MPs_3D.txt");

    for (unsigned int i=0; i<vpMPs.size(); i++)
    {
        if (vpMPs[i]->isBad())
            continue;

        cv::Mat pos = vpMPs[i]->GetWorldPos();
        cv::Mat normal = vpMPs[i]->GetNormal();

        MPs_3Ds << "vector<float> MP" << i << ";";

        Vec3f vPos = pos;
        float fPos1 = vPos.val[0];
        float fPos2 = vPos.val[1];
        float fPos3 = vPos.val[2];

        Vec3f vNor = normal;
        float fNor1 = vNor.val[0];
        float fNor2 = vNor.val[1];
        float fNor3 = vNor.val[2];
        //cout << fPos1 << "-" << fPos2 << "-" << fPos3 << ":" << fNor1 << "-" << fNor2 << "-" << fNor3 << endl;

        // Save 3D position
        MPs_3Ds << "\tMP" << i << ".push_back(" << fPos1 << ");";
        MPs_3Ds << "\tMP" << i << ".push_back(" << fPos2 << ");";
        MPs_3Ds << "\tMP" << i << ".push_back(" << fPos3 << ");";

        // Save normal vector
        MPs_3Ds << "\tMP" << i << ".push_back(" << fNor1 << ");";
        MPs_3Ds << "\tMP" << i << ".push_back(" << fNor2 << ");";
        MPs_3Ds << "\tMP" << i << ".push_back(" << fNor3 << ");";


        /*for (unsigned int j=0; j<3; j++)
        {
            float pos_i = pos.at<char>(j,0);

            MPs_3Ds << "\tMP" << i << ".push_back(" << pos_i << ");";

        }*/

        MPs_3Ds << "\tvMPs.push_back(MP" << i << ");";
        MPs_3Ds << endl;

    }

    MPs_3Ds.close();

    cout << "OK (" << vpMPs_3D.size() << " MapPoints)" << endl;

    cout << "  LastFrame ... ";

    string filename_F = "../../../evaluation/map/Frame.bin";
    const char * filename_F_char = filename_F.c_str();

    Frame* pF = &mpTracker->mCurrentFrame;

    FILE * pFile_F;
    pFile_F = fopen(filename_F_char,"w+b");
    fwrite(pF,sizeof(Frame),1,pFile_F);
    fclose(pFile_F);

    cout << "OK (" << 1 << " Frame, " << sizeof(Frame) << " Bytes" << ")" << endl;


}


/*void System::LoadMap(const string &filename_MPs, const string &filename_KFs)
{
    cout << endl << endl << "Lading map:" << endl;

    const char * filename_F_char = "../../../evaluation/map/Frame.bin";
    const char * filename_KFs_char = filename_KFs.c_str();
    const char * filename_MPs_char = filename_MPs.c_str();

    FILE * pFile_F;
    FILE *pFile_KFs;
    FILE *pFile_MPs;

    pFile_F = fopen ( filename_F_char , "r+b" );
    if (pFile_F==NULL)
    {
        cout << "Unable to open Frame file" << endl;
        return;
    }
    fseek(pFile_F,0,SEEK_END);
    size_t lSize_F = ftell(pFile_F);
    fclose(pFile_F);

    pFile_KFs = fopen ( filename_KFs_char , "r+b" );
    if (pFile_KFs==NULL)
    {
        cout << "Unable to open KeyFrames file" << endl;
        return;
    }
    fseek(pFile_KFs,0,SEEK_END);
    size_t lSize_KFs = ftell(pFile_KFs);
    fclose(pFile_KFs);

    pFile_MPs = fopen ( filename_MPs_char , "r+b" );
    if (pFile_MPs==NULL)
    {
        cout << "Unable to open MapPoints file" << endl;
        return;
    }
    fseek(pFile_MPs,0,SEEK_END);
    size_t lSize_MPs = ftell(pFile_MPs);
    fclose(pFile_MPs);

    long size_F = sizeof(Frame);
    long size_KF = sizeof(KeyFrame);
    long size_MP = sizeof(MapPoint);

    unsigned int nFs = lSize_F/size_F;
    unsigned int nKFs = lSize_KFs/size_KF;
    unsigned int nMPs = lSize_MPs/size_MP;




    pFile_F = fopen ( filename_F_char , "r+b" );
    cout << "  " << nFs << " Frames to load ... ";

    Frame buffer_F = Frame();
    int ressult_F = fread(&buffer_F,sizeof(Frame),1,pFile_F);

    fclose(pFile_F);

    cout << "OK(" << ressult_F << ")" << endl;



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
    cout << " OK(" << nOkMPs << ")" << endl;
}*/



int System::GetTrackingState()
{
    unique_lock<mutex> lock(mMutexState);
    return mTrackingState;
}

vector<MapPoint*> System::GetTrackedMapPoints()
{
    unique_lock<mutex> lock(mMutexState);
    return mTrackedMapPoints;
}

vector<cv::KeyPoint> System::GetTrackedKeyPointsUn()
{
    unique_lock<mutex> lock(mMutexState);
    return mTrackedKeyPointsUn;
}

} //namespace ORB_SLAM
