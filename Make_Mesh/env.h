/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#define __ENVIRONMENT__

#ifndef __SEARCH__
#include "search.h"
#endif

#ifndef __DISPLAY__
#include "display.h"
#endif

#include "Eigen/Core"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include<iostream>
#include<string>

using namespace Eigen;
using namespace std;

class TEnvironment {
 public:
  TEnvironment(); 
  ~TEnvironment(); 
  void SetTarget(); // set the target figure
  void SetInit();   // initial setting
  void DoIt();      // main procedure
  void WriteData();

  void Display_top(); // draw top tile shapes 
  void Write_Tile();  // write top tilig shapes to a file

  int fN; // the number of vertices
  VectorXd fW; // the coordinates of the target shape
  TSearch *tSearch; // a class for the search procedure

  clock_t fTimeStart, fTimeEnd; // for execution tme
  char *fInstanceName;          // name of the target figure
  char *fOutputName = NULL;     // output file name
  char *fFileNameWeight = NULL;  // file name for weight setting

  TDisplay tDisplay;  // a class of drawing the tile shapes
  int fFlagDisplay;   // 0: not display, 1: display tile figures
  double fAlpha;      // parameter
  double fBeta;      // parameter
};

#endif

  
