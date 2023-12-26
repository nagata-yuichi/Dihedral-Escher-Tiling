/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#define __ENVIRONMENT__

#ifndef __DISPLAY__
#include "display.h"
#endif

#include "Eigen/Core"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>


using namespace Eigen;
using namespace std;

class TEnvironment {
 public:
  TEnvironment(); 
  ~TEnvironment(); 
  void SetInit();   // Initial setup
  void ReadTile();  // Read tile shapes from the specified file
  void DoIt();      // Main procedure
  void Display_top(); // Display tile shpaes and tiling results
  void SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd ); // Sort data
  
  // int fN, fNa, fNv;        // n, n_v
  int fN1, fN1_in, fNN1;
  int fN2, fN2_in, fNN2;  
  VectorXd fW1;        // w
  VectorXd fWW1;       // \tilde{w}
  MatrixXi fKK1i;  
  VectorXd fW2;        // w
  VectorXd fWW2;       // \tilde{w}
  MatrixXi fKK2i;  
  
  TDisplay *tDisplay; // // a class of drawing the tile shapes
  char *fResultFileName; // input file name 

  int fNum_top;       // the number of top solutions stored in the EST
  double fTime;       // execution time 
  double *fEval_top; // the evaluating values of the top solutions 
  int *fIndex_top;   // [i] -> index of the i-th best solution
  VectorXd *fU_top;  // the top solutions (tile shape)
  VectorXd *fUa_top;  // the top solutions (tile shape)
  VectorXd *fUU1_top; // the top solutions (mesh representation)
  VectorXd *fUU2_top; // the top solutions (mesh representation) 
  int **ffn_top;     // [s] -> fi[] of the tile shape fU_top[s]
  int *fNv_top;      // [s] -> fNv of the tile shape fU_top[s]
  int *fIH_top;      // [s] -> IH type of the tile shape fU_top[s]
  int *fN_top;       //
  int *fNa_top;       // 
  int *fSt0_top;     // [s] -> st0 of the tile shape fU_top[s]
  int *fSt1_top;     // [s] -> st1 of the tile shape fU_top[s]
  int *fSt2_top;     // [s] -> st2 of the tile shape fU_top[s]
  int *fReverse_top;

  double fScale_T;
};

  
#endif
