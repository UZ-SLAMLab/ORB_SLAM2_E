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

#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <limits>

#include <opencv2/highgui/highgui.hpp>

namespace ORB_SLAM2
{

using namespace std;
using namespace cv;


class Statistics {
 public:
    Statistics(string file);
    ~Statistics();

    void OpenFile(bool newFile);
    void CloseFile();

    void ColumnHeadersT();
    void ColumnHeadersPR(); //Precision and Recall
    void ColumnHeadersLM();
    void ColumnHeadersMPRej();
    void ColumnHeadersReloc();
    void ColumnHeadersRelocS1();
    void ColumnHeadersPnP();
    void ColumnHeadersNLOpt();

    void FillTable();

    void NewRow(int a, int b, int c, int d, int e, int f);
    void NewRow(int a, int b, int c, int d, int e, int f, int g);
    void NewRow(int a, int b, int c, int d, int e, int f, int g, int h);
    void NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i);
    void NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j);
    void NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j, int k);
    void NewRow(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j, int k, int l, int m, int n, int p, int q, int r, int s);

    void NewIteration(int nIt);
    void AddValue(int a);
    void AddValueS1(int a);
    void AddValueFl(float a);
    void AddText(string a);
    void EmptyCols(int a);
    void NewLine();

    void BarGraph();

    void StartChrono();
    void StopChrono();

    vector<string> GetHeaders();
    vector<int> GetPercentajes();

 private:
    string filename;

    ofstream o;

    fstream i;

    vector<vector<int> > MapStat;

    int A[10][5];

    float t1,t2,freq,time;

    vector<string> Heads;

    vector<vector<int> > SMAvalues;
    vector<int> SMApercentajes;

    void UpdateSMA_Tracking(vector<int> input);
    void UpdateSMA_Mapping(vector<int> input);
    void UpdateSMA_NCMP(vector<int> input);

};

} //namespace ORB_SLAM2

#endif // STATISTICS_H
