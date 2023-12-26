/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#define __DISPLAY__

#ifndef  __OPERATOR_ADD__
#include "operator_additional.h"
#endif

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

  void SetInit( int N1, int N1_in, VectorXd ww1, MatrixXi kk1i, int N2, int N2_in, VectorXd ww2, MatrixXi kk2i ); // Set data of the goal mesh
  void Set_range1(); // Get the maximum and minimum values of the coordinates of the goal shape
  void Set_range2(); // Get the maximum and minimum values of the coordinates of the goal shape

  void Goal_mesh1(); // Display the goal mesh
  void Goal_mesh2(); // Display the goal mesh
  void Tile_mesh1( VectorXd uu ); // Display the tile mesh
  void Tile_mesh2( VectorXd uu ); // Display the tile mesh

  void Dihedral( VectorXd u, VectorXd ua, int st0, int nv, int ffn[], VectorXd uu1, VectorXd uu2 ); // デフォルト引数で統一したい
  void Dihedral2( VectorXd u, VectorXd ua, int st0, VectorXd uu1, VectorXd uu2 );

  void Tiling_IH4( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH5( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH6( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH1( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH2( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH3( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result
  void Tiling_IH27( VectorXd& u, VectorXd& ua, int *fn, int nv ); // Display the tiling result

  void Tiling_IH4( VectorXd u, VectorXd ua, int *fn, int nv, VectorXd uu1, VectorXd uu2 );
  void Tiling_IH5( VectorXd u, VectorXd ua, int *fn, int nv, VectorXd uu1, VectorXd uu2 );
  void Tiling_IH6( VectorXd u, VectorXd ua, int *fn, int nv, VectorXd uu1, VectorXd uu2 );
  

  void All_tiles( VectorXd* uu1_top, VectorXd* uu2_top );

  TOperator_ADD *tOperator;
  
  int fNumDisplay;
  int fOpen_flag_tiling;  
  double fWidth, fHeight;
  int fxs1, fxe1, fys1, fye1, fxs2, fxe2, fys2, fye2;
  int fN, fN1, fN2;
  int fN1_in, fN2_in;
  int fNN1, fNN2;
  VectorXd fW1, fW2;
  VectorXd fW1_in, fW2_in;
  VectorXd fWW1, fWW2;
  MatrixXi fKK1i, fKK2i;

  double fScale_T;


 int inline ix( int i ){ // positoin of the x-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i;
  }

  int inline iy( int i ){ // positoin of the y-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i+1;
  }
  
  int inline ix1( int i ){ // positoin of the x-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN1;
    else if ( i >= fN1 )
      i -= fN1;
    return 2*i;
  }

  int inline iy1( int i ){ // positoin of the y-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN1;
    else if ( i >= fN1 )
      i -= fN1;
    return 2*i+1;
  }

  int inline ix2( int i ){ // positoin of the x-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN2;
    else if ( i >= fN2 )
      i -= fN2;
    return 2*i;
  }

  int inline iy2( int i ){ // positoin of the y-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN2;
    else if ( i >= fN2 )
      i -= fN2;
    return 2*i+1;
  }

  int inline iix1( int i ){ // positoin of the x-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN1;
    else if ( i >= fNN1 )
      i -= fNN1;
    return 2*i;
  }

  int inline iiy1( int i ){ // positoin of the y-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN1;
    else if ( i >= fNN1 )
      i -= fNN1;
    return 2*i+1;
  }

   int inline iix2( int i ){ // positoin of the x-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN2;
    else if ( i >= fNN2 )
      i -= fNN2;
    return 2*i;
  }

  int inline iiy2( int i ){ // positoin of the y-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN2;
    else if ( i >= fNN2 )
      i -= fNN2;
    return 2*i+1;
  }
};

#endif
