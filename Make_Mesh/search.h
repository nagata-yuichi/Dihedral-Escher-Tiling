/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_E__
#define __SEARCH_E__

#ifndef __DISPLAY__
#include "display.h"
#endif

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <ncurses.h>

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;

class TSearch {
 public:
  TSearch( int N );
  ~TSearch(); 
  void SetParameter();

  void SetInit( VectorXd w );
  void DoIt();
  void SetRegion();
  bool Check_InnerResion( double px, double py );
  void SetInitialInnerPoints();
  void ModiInnerPoints();
  void SetK(); 
  bool Check_Internal_link( int s1, int s2 );
  bool Intersection( int a1, int a2, int b1, int b2 );
  void SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd );
  int iix( int i );
  int iiy( int i );


  int fN;  // ターゲット図形の点の数
  VectorXd fW;  // ターゲット図形の座標ベクトル
  MatrixXd fG;      // Eでは使わない（display用）

  int fN_in;      // the number of the inner points 
  int fNN;        // the number of all points   
  VectorXd fW_in; // the positions of the Inner points of W
  VectorXd fWW;   // the positions of all points of W

  double fLengthX, fLengthY;
  double fX_min, fY_min;
  MatrixXi fInnerRegion; 
  double fL, fL_min, fL_max;
  double fAlpha;
  double fBeta;
  MatrixXi fKKi;
  MatrixXi fNei;      
  VectorXi fNei_Size; 
  VectorXd fV_short = VectorXd(2);
  VectorXd fV_long = VectorXd(2);
  TDisplay tDisplay;


  int inline ix( int i ){ // index iの点のx座標の位置
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i;
  }
  int inline iy( int i ){ // index iの点のy座標の位置
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i+1;
  }

};


#endif

  
