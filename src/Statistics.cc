/**
* This file is part of a modified version of ORB-SLAM, intended to carry
* on a performance evaluation of the software with Medical Endoscopy Sequences.
*
* Copyright (C) 2016 Iñigo Cirauqui-Viloria <cirauqui21 at gmail dot com> (University of Zaragoza)
*
* For more information of the original software see <http://webdiis.unizar.es/~raulmur/orbslam/>
*
* This version of ORB-SLAM is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
*
* This version of ORB-SLAM is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this version of ORB-SLAM. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <limits>

#include <opencv2/highgui/highgui.hpp>

#include "Statistics.h"


using namespace std;
using namespace cv;

namespace ORB_SLAM2
{


Statistics::Statistics(string file)
{
    filename = file;
}


Statistics::~Statistics()
{
}


void Statistics::OpenFile(bool newFile)
{
    if (newFile) {
        o.open(filename.c_str());
    }
    if (!newFile) {
        o.open(filename.c_str(), ios_base::app);
    }
}


void Statistics::CloseFile()
{
    o.close();
}


void Statistics::ColumnHeadersT()  // Prints info about the tracking
{
    o << "Frame\t"      // Frame ID
      << "nORBs\t"      // Features in current Frame
      << "MPinF\t"      // MPs in current image
      << "Tracked\t"    // Tracked MPs
      << "noClose\t"    // Reprojected MPs wo close feature in Frame
      << "nNoEq\t"      // Reprojected MPs wo similar ORB
      << "noBA\t";      // Bad Point by BA

    string text;
    text = "Frame";   Heads.push_back(text);
    text = "nORBs";   Heads.push_back(text);
    text = "MPinF";   Heads.push_back(text);
    text = "Tracked"; Heads.push_back(text);
    text = "noClose"; Heads.push_back(text);
    text = "noORBeq"; Heads.push_back(text);
    text = "BA";      Heads.push_back(text);
}


void Statistics::ColumnHeadersLM()  // Saves info about the stress test of the MPs in the 3 KFs since they are born
{
    o << "Frame\t"  // Frame ID
      << "KF\t"     // KeyFrame ID
      << "Life0\t"  // MPs created in last KF
      << "Life1\t"  // MPs alife for 2 KF
      << "Life2\t"  // MPs alife for 3 KF
      << "Life3\t"  // MPs alife for 4 KF, they become permanent
      << "Cond1\t"  // Not matched in 25% frames
      << "Cond2\t"  // Not observed in 3 KFs
      << "BA\t"     // MPs rejected by Bundle adjustment
      << "Adult\t"  // New adult MP
      << "fpkf1c1\t"
      << "fpkf1ba\t"
      << "fpkf2c1\t"
      << "fpkf2c2\t"
      << "fpkf2ba\t"
      << "fpkf3c1\t"
      << "fpkf3c2\t"
      << "fpkf3ba\t";

    string text;
    text = "Frame"; Heads.push_back(text);
    text = "KF";    Heads.push_back(text);
    text = "kf1c1"; Heads.push_back(text);
    text = "kf1ba"; Heads.push_back(text);
    text = "kf2c1"; Heads.push_back(text);
    text = "kf2c2"; Heads.push_back(text);
    text = "kf2ba"; Heads.push_back(text);
    text = "kf3c1"; Heads.push_back(text);
    text = "kf3c2"; Heads.push_back(text);
    text = "kf3ba"; Heads.push_back(text);
    text = "Life0"; Heads.push_back(text);
    text = "Life1"; Heads.push_back(text);
    text = "Life2"; Heads.push_back(text);
    text = "Life3"; Heads.push_back(text);
}


void Statistics::ColumnHeadersMPRej()
{
    o << "Frame\t"      // Frame ID
      << "KF\t"         // KeyFrame ID
      << "fORBs\t"      // ORBs available
      << "EpipOK\t"     // Epipolar geometry ok
      << "Depth\t"      // Negative depth in any of the 2 KFs
      << "Paral\t"      // Not enought parallax
      << "ReprE\t"      // Too much reprojection error
      << "Scale\t"      // Scale inconsistent
      << "BA\t"         // BadPoint by BA
      << "newMP\t";     // New MPs.

    string text;
    text = "Frame";  Heads.push_back(text);
    text = "KF";     Heads.push_back(text);
    text = "ORBs";   Heads.push_back(text);
    text = "epipOK"; Heads.push_back(text);
    text = "Depth";  Heads.push_back(text);
    text = "Paral";  Heads.push_back(text);
    text = "ReprE";  Heads.push_back(text);
    text = "Scale";  Heads.push_back(text);
}


void Statistics::ColumnHeadersReloc()
{
    /*
    o << "kfcan\t"      // Frame ID
      << "acc15\t"         // KeyFrame ID
      << "acc4 \t"      // ORBs available
      << "nInl \t"     // Epipolar geometry ok
      << "tPnP \t"      // Negative depth in any of the 2 KFs
      << "nAdd \t"      // Not enought parallax
      << "1okR \t"      // Too much reprojection error
      << "1okNR\t"      // Scale inconsistent
      << "tS1  \t"
      << "1okR \t"      // Too much reprojection error
      << "1okNR\t"      // Scale inconsistent
      << "tS2  \t"
      << "1okR \t"      // Too much reprojection error
      << "1okNR\t"      // Scale inconsistent
      << "tS3  \t"
      << "final\t"         // BadPoint by BA
      << "ress \t"     // New MPs.
      << endl;
      */


    o << "KF_candidates\t"      // Candidate KeyFrames
      << "KF_15matches \t"      // KeyFrames with 15 matches or more
      << "KF_04matches \t"      // KeyFrames with 4 matches or more
      << ".            \t"      // KeyFrames with 4 matches or more
      << "Inliers_PnP_R\t"      // Inliers original PnP
      << "Time_PnP_R   \t"      // Time original PnP
      << "Inliers_PnP_D\t"      // Inliers new PnP
      << "Time_PnP_D   \t"      // Time new PnP
      << "InliersAddBoW\t"      // Inliers added BoW
      << "InliersAddMap\t"      // Inliers added Map Projection
      << "nGoodR\t"             // Final Inliers original NL Optimization
      << "timeR\t"              // Time original NL Optimization
      << "nGoodD\t"             // Final Inliers new NL Optimization
      << "timeD\t"              // Time new NL Optimization
      << endl;


    string text;
    text = "Frame";  Heads.push_back(text);
    text = "KF";     Heads.push_back(text);
    text = "ORBs";   Heads.push_back(text);
    text = "epipOK"; Heads.push_back(text);
    text = "Depth";  Heads.push_back(text);
    text = "Paral";  Heads.push_back(text);
    text = "ReprE";  Heads.push_back(text);
    text = "Scale";  Heads.push_back(text);
}

void Statistics::ColumnHeadersRelocS1()
{
    o << "CurF\t"      // Frame ID
      << "KF15\t"         // KeyFrame ID
      << "KF0y\t"      // ORBs available
      << "KFnn\t"     // Epipolar geometry ok
      << "matc\t"      // Negative depth in any of the 2 KFs
      << endl;

    string text;
    text = "CurF";  Heads.push_back(text);
    text = "KF15";  Heads.push_back(text);
    text = "KF0y";  Heads.push_back(text);
    text = "KFnn";  Heads.push_back(text);
    text = "matc";  Heads.push_back(text);
}


void Statistics::ColumnHeadersPR()  // Prints info about the tracking
{
    o << "Fr\t"      // Frame ID
      << "TFP\t"      // KeyFrame ID
      << "TP\t"
      << endl;

    string text;
    text = "Fr";     Heads.push_back(text);
    text = "TFP";     Heads.push_back(text);
    text = "TP";   Heads.push_back(text);
}


void Statistics::ColumnHeadersPnP()
{
    o << "Frame\t"      // Frame ID
      << "KF\t"         // KeyFrame ID
      << "fORBs\t"      // ORBs available
      << "EpipOK\t"     // Epipolar geometry ok
      << "Depth\t"      // Negative depth in any of the 2 KFs
      << "Paral\t"      // Not enought parallax
      << "ReprE\t"      // Too much reprojection error
      << "Scale\t"      // Scale inconsistent
      << "BA\t"         // BadPoint by BA
      << "newMP\t";     // New MPs.

    string text;
    text = "Frame";  Heads.push_back(text);
    text = "KF";     Heads.push_back(text);
    text = "ORBs";   Heads.push_back(text);
    text = "epipOK"; Heads.push_back(text);
    text = "Depth";  Heads.push_back(text);
    text = "Paral";  Heads.push_back(text);
    text = "ReprE";  Heads.push_back(text);
    text = "Scale";  Heads.push_back(text);
}


void Statistics::ColumnHeadersNLOpt()
{
    o << "Frame\t"      // Frame ID
      << "KF\t"         // KeyFrame ID
      << "fORBs\t"      // ORBs available
      << "EpipOK\t"     // Epipolar geometry ok
      << "Depth\t"      // Negative depth in any of the 2 KFs
      << "Paral\t"      // Not enought parallax
      << "ReprE\t"      // Too much reprojection error
      << "Scale\t"      // Scale inconsistent
      << "BA\t"         // BadPoint by BA
      << "newMP\t";     // New MPs.

    string text;
    text = "Frame";  Heads.push_back(text);
    text = "KF";     Heads.push_back(text);
    text = "ORBs";   Heads.push_back(text);
    text = "epipOK"; Heads.push_back(text);
    text = "Depth";  Heads.push_back(text);
    text = "Paral";  Heads.push_back(text);
    text = "ReprE";  Heads.push_back(text);
    text = "Scale";  Heads.push_back(text);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f)  // 6 inputs
{
    o  << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g)  // 7 inputs
{
    o  << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f); input.push_back(g);

    UpdateSMA_Tracking(input);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g, int h)  // 8 inputs
{
    o << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f); input.push_back(g); input.push_back(h);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i)  // 9 inputs
{
    o << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t" << i << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f); input.push_back(g); input.push_back(h); input.push_back(i);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j)  // 10 inputs
{
    o << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t" << i << "\t" << j << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f); input.push_back(g); input.push_back(h); input.push_back(i); input.push_back(j);

    UpdateSMA_NCMP(input);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j, int k)  // 11 inputs
{
    o << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t" << i << "\t" << j << "\t" << k << "\t";

    vector<int> input;
    input.push_back(a); input.push_back(b); input.push_back(c); input.push_back(d); input.push_back(e);
    input.push_back(f); input.push_back(g); input.push_back(h); input.push_back(i); input.push_back(j);
    input.push_back(k);
}


void Statistics::NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j, int k, int l, int m, int n, int p, int q, int r, int s)  // 18 inputs
{
    o << endl << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t" << i << "\t"
              << j << "\t" << k << "\t" << l << "\t" << m << "\t" << n << "\t" << p << "\t" << q << "\t" << r << "\t" << s << "\t";

    vector<int> input;
    input.push_back(a);//1
    input.push_back(b);//2
    input.push_back(k);//11
    input.push_back(l);//12
    input.push_back(m);//13
    input.push_back(n);//14
    input.push_back(p);//15
    input.push_back(q);//16
    input.push_back(r);//17
    input.push_back(s);//18
    input.push_back(c);//3
    input.push_back(d);//4
    input.push_back(e);//5
    input.push_back(f);//6

    UpdateSMA_Mapping(input);
}


void Statistics::NewIteration(int nIt)
{
    OpenFile(0);
    o << endl
      << endl
      << "It " << nIt << endl;
    ColumnHeadersReloc();
    CloseFile();
}


void Statistics::AddValue(int a)
{
    OpenFile(0);
    o << a << "\t";
    CloseFile();
}


void Statistics::AddValueS1(int a)
{
    OpenFile(0);
    EmptyCols(1);
    o << a << endl;
    CloseFile();
}


void Statistics::AddValueFl(float a)
{
    OpenFile(0);
    o << setprecision(5) << a << "\t";
    CloseFile();
}


void Statistics::AddText(string a)
{
    OpenFile(0);
    o << a << "\t";
    CloseFile();
}


void Statistics::EmptyCols(int a)
{
    OpenFile(0);
    for (int aa=0; aa<a; aa++)
        o << " \t";
    CloseFile();
}


void Statistics::NewLine()
{
    OpenFile(0);
    o << endl;
    CloseFile();
}


/*
void Statistics::BarGraph()
{
    i.open(filename.c_str());
    i.close();

    for (int i=0; i<10; i++){
        for (int j=0; j<5; j++){
            if (j==0)
                A[i][j] = i+1;
            else if (j!=0)
                A[i][j] = rand() % 100;
        }
    }
    for (int i=0; i<10; i++){
        for (int j=0; j<5; j++){
            cout << " " << A[i][j];
            if (j==4)
                cout << endl;
        }
    }

    namedWindow("Histogram",CV_WINDOW_AUTOSIZE);

    int heigh = 500;
    int width = 200;
    int space = 5;
    int jmax = 5 - 1;
    int barWidth = ((width/(5-1))-space-space/jmax);

    Mat Hist(Size(width,heigh),CV_8UC3);
    for (int i=0; i<10; i++) {
        Point2f pt0, pt1, pt2;
        pt0.x = 1;
	pt0.y = 1;
        rectangle (Hist,pt0,500*pt0,cv::Scalar(0,0,0),-1);
        for (int j=0; j<5; j++) {
                pt1.x = barWidth*(j-1) + space*j;
                pt1.y = 500;
                pt2.x = barWidth*(j) + space*j;
                pt2.y = pt1.y - ((A[i][j])*5);
                rectangle (Hist,pt1,pt2,cv::Scalar(rand()&255,rand()&255,rand()&255),-1);
        }
        imshow("Histogram",Hist);
        //waitKey(0);
    }
}
*/


void Statistics::StartChrono()
{
    freq = getTickFrequency();
    t1 = getTickCount();
}

void Statistics::StopChrono()
{
    t2 = getTickCount();
    time = (t2-t1)/freq;

    cout << endl << " Time:\t" << time << endl;
}


void Statistics::UpdateSMA_Tracking(vector<int> input)
{
    if (SMAvalues.empty()) {
        for (unsigned int i=0; i<input.size(); i++) {
            vector<int> value;
            value.push_back(input[i]);
            SMAvalues.push_back(value);
        }
    }
    else {
        for (unsigned int i=0; i<SMAvalues.size(); i++) {
            if (SMAvalues[i].size()<20)
                SMAvalues[i].push_back(input[i]);
            else {
                vector<int> temp;
                for (unsigned int j=1; j<SMAvalues[i].size(); j++)
                    temp.push_back(SMAvalues[i][j]);
                temp.push_back(input[i]);
                SMAvalues[i] = temp;
            }
        }
    }

    SMApercentajes.clear();
    for (unsigned int i=0; i<SMAvalues.size(); i++) {
        int mean;
        int sum = 0;
        int last = SMAvalues[i].size()-1;
        if (i==0) mean = SMAvalues[i][last];
        else {
            for (unsigned int j=0; j<SMAvalues[i].size(); j++)
                sum += SMAvalues[i][j];
            mean = sum / 20;
        }

        if (i<3)
            SMApercentajes.push_back(mean);
        else {
            double perc = 100 * mean / SMApercentajes[2];
            SMApercentajes.push_back(perc);
        }
    }
}


void Statistics::UpdateSMA_Mapping(vector<int> input)
{
    if (SMAvalues.empty()) {
        for (unsigned int i=0; i<input.size(); i++) {
            vector<int> value;
            value.push_back(input[i]);
            SMAvalues.push_back(value);
        }
    }
    else {
        for (unsigned int i=0; i<SMAvalues.size(); i++) {
            if (SMAvalues[i].size()<20)
                SMAvalues[i].push_back(input[i]);
            else {
                vector<int> temp;
                for (unsigned int j=1; j<SMAvalues[i].size(); j++)
                    temp.push_back(SMAvalues[i][j]);
                temp.push_back(input[i]);
                SMAvalues[i] = temp;
            }
        }
    }

    vector<int>means;
    for (unsigned int i=0; i<SMAvalues.size(); i++) {
        int sum = 0;
        int mean;
        int last = SMAvalues[i].size()-1;
        if (i<2) mean = SMAvalues[i][last];
        else {
            for (unsigned int j=0; j<SMAvalues[i].size(); j++)
                sum += SMAvalues[i][j];
            mean = sum / 20;
        }
        means.push_back(mean);
    }

    double perc;

    SMApercentajes.clear();
    for (unsigned int i=0; i<means.size(); i++) {
        if (i<2) {
            perc = means[i];
            SMApercentajes.push_back(perc);
        }
        else if (i<4) {
            int cat = means[2] + means[3];
            if (cat==0) perc = 0; else perc = 100 * means[i] / cat;
            SMApercentajes.push_back(perc);
        }
        else if (i<7) {
            int cat = means[4] + means[5] + means[6];
            if (cat==0) perc = 0; else perc = 100 * means[i] / cat;
            SMApercentajes.push_back(perc);
        }
        else if (i<10) {
            int cat = means[7] + means[8] + means[9];
            if (cat==0) perc = 0; else perc = 100 * means[i] / cat;
            SMApercentajes.push_back(perc);
        }
        else if (i<11) {
            perc = means[i];
            SMApercentajes.push_back(perc);
        }
        else if (i<13) {
            int cat = means[10];
            if (cat==0) perc = 0; else perc = 100 * means[i] / cat;
            SMApercentajes.push_back(perc);
        }
        else {
            int cat = means[10];
            int remains = means[12] - means[7] - means[8] - means[9];
            if (cat==0) perc = 0; else perc = 100 * remains / cat;
            SMApercentajes.push_back(perc);
        }
    }
}


void Statistics::UpdateSMA_NCMP(vector<int> input)
{
    if (SMAvalues.empty()) {
        for (unsigned int i=0; i<input.size(); i++) {
            vector<int> value;
            value.push_back(input[i]);
            SMAvalues.push_back(value);
        }
    }
    else {
        for (unsigned int i=0; i<SMAvalues.size(); i++) {
            if (SMAvalues[i].size()<20)
                SMAvalues[i].push_back(input[i]);
            else {
                vector<int> temp;
                for (unsigned int j=1; j<SMAvalues[i].size(); j++)
                    temp.push_back(SMAvalues[i][j]);
                temp.push_back(input[i]);
                SMAvalues[i] = temp;
            }
        }
    }

    vector<int>means;
    for (unsigned int i=0; i<SMAvalues.size(); i++) {
        int sum = 0;
        int mean;
        int last = SMAvalues[i].size()-1;
        if (i<2) mean = SMAvalues[i][last];
        else {
            for (unsigned int j=0; j<SMAvalues[i].size(); j++)
                sum += SMAvalues[i][j];
            mean = sum / 20;
        }
        means.push_back(mean);
    }

    double perc;

    SMApercentajes.clear();
    for (unsigned int i=0; i<means.size(); i++) {
        if (i<3) {
            perc = means[i];
            SMApercentajes.push_back(perc);
        }
        else if (i<4) {
            if (means[2]==0) perc = 0; else perc = 100 * means[i] / means[2];
            SMApercentajes.push_back(perc);
        }
        else {
            int tot = means[4] + means[5] + means[6] + means[7];
            if (tot==0) perc = 0; else perc = 100 * means[i] / tot;
            SMApercentajes.push_back(perc);
        }
    }
}


vector<string> Statistics::GetHeaders()
{
    vector<string> output = Heads;
    return output;
}


vector<int> Statistics::GetPercentajes()
{
    vector<int> output;
    output = SMApercentajes;
    return output;
}


} //namespace ORB_SLAM2

