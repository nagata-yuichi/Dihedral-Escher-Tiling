/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#define __DISPLAY__

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include  "xw.h"
#include  <X11/Xlib.h>
#include <unistd.h>

using namespace Eigen;
using namespace std;

class TDisplay {
 public:
  TDisplay(); 
  ~TDisplay(); 
  void Disp_bit( VectorXd& u , VectorXd& w, int *fi, int nv );
  void Disp_AdjacentRelation( int n, VectorXd w, MatrixXi K ); 
  void Disp_AdjacentRelation2( int n, VectorXd w, MatrixXi K ); 

  void Disp_goal( VectorXd w , MatrixXd G );
  void Disp_tile( VectorXd& w, VectorXd& u );
  void Disp_all( VectorXd* u , VectorXd& w, int** fi, int* nv );
  void Disp_all_sup( VectorXd* u , VectorXd& w, int** fi, int* nv );

  int fBoxL; 
};

#endif
