/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_ADD__
#define __OPERATOR_ADD__

#ifndef  __OPERATOR_IR__
#include "operator_IR.h"
#endif


class TOperator_ADD : public TOperator_IR {
 public:
  TOperator_ADD( int N1, int N1_in, int N2, int N2_in );
  ~TOperator_ADD(); 
  void SetParameter();
  void SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2,  VectorXd w2_in, MatrixXi kk2i );
  void SetCell1_1_ling();
  void SetCell2_1_ling();
  void SetNei1_Sub();
  void SetNei2_Sub();

  void Find_Triangle1( double xp, double yp, VectorXd& uup, int& index_triangle, double& alpha, double& beta );
  void Find_Triangle2( double xp, double yp, VectorXd& uup, int& index_triangle, double& alpha, double& beta );
  bool Is_inside_Up1( double x, double y, VectorXd& uup );
  bool Is_inside_Up2( double x, double y, VectorXd& uup );
      
  double Cal_Area1( VectorXd& uu );
  double Cal_Area2( VectorXd& uu );

  // MatrixXi fNei1, fNei2; // TOperator_Iで定義
  // VectorXi fNei1_Size, fNei2_Size; 
  MatrixXi fTriangleList1,fTriangleList2; // triangle list 
  int fTriangle1_num, fTriangle2_num; 
  MatrixXi fEdgeNei1_1, fEdgeNei1_2; // edgeを含むtriangleの残りの頂点
  MatrixXi fEdgeNei2_1, fEdgeNei2_2; // -1の場合はなし
  int fNumOfCell1, fNumOfCell2;
  MatrixXi *fCell1_Edge, *fCell2_Edge;
  VectorXi fCell1_Size, fCell2_Size;

};
#endif

