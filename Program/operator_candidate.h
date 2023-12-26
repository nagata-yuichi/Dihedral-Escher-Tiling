/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_CANDIDATE__
#define __OPERATOR_CANDIDATE__

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;


class TOperator_Candidate
{
 public:
  TOperator_Candidate( int N1, int N2 );
  ~TOperator_Candidate();
  void SetParameter();
  void SetInit( VectorXd w1,  VectorXd w2 ); 

  void Set_eval_JS();
  double Diff_segment_rotation( VectorXd segment1, VectorXd segment2, int len );
  double Diff_segment_parallel( VectorXd segment1, VectorXd segment2, int len );
  double Diff_segment_rotation_GR( VectorXd segment1, VectorXd segment2, int len );

  int fN1, fN2;
  int fIH, fNa, fNv;
  //  int *fn;
  int fReverse1, fReverse2;
  VectorXd fW1, fW2; 
  
  double ***fEval_J12_edge;
  double ***fEval_J1_edge;
  double ***fEval_J2_edge;
  double **fEval_S1_edge;
  double **fEval_S2_edge;
  double ***fEval_J12_GR_edge;
  double ***fEval_J1_GR_edge;
  double ***fEval_J2_GR_edge;


  int inline ix1( int i ){ // index iの点のx座標の位置
    if( i < 0 )
      i += fN1;
    else if ( i >= fN1 )
      i -= fN1;
    return 2*i;
  }
  int inline iy1( int i ){ // index iの点のy座標の位置
    if( i < 0 )
      i += fN1;
    else if ( i >= fN1 )
      i -= fN1;
    return 2*i+1;
  }

  int inline ix2( int i ){ // index iの点のx座標の位置
    if( i < 0 )
      i += fN2;
    else if ( i >= fN2 )
      i -= fN2;
    return 2*i;
  }
  int inline iy2( int i ){ // index iの点のy座標の位置
    if( i < 0 )
      i += fN2;
    else if ( i >= fN2 )
      i -= fN2;
    return 2*i+1;
  }
};

#endif
