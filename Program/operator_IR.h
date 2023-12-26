/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_IR__
#define __OPERATOR_IR__

#ifndef  __OPERATOR_I__
#include "operator_I.h"
#endif

class TOperator_IR : public TOperator_I {
 public:
  TOperator_IR( int N1, int N1_in, int N2, int N2_in );
  ~TOperator_IR(); 
  void SetParameter() override;
  void SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2, VectorXd w2_in, MatrixXi kk2i ) override;
  void Cal_u_SGD() override;
  void Cal_u_SGD( int st0, int st1, int st2 ) override;
  double SGD_main( VectorXd& u1, VectorXd& u2, VectorXd& u, VectorXd& ua, VectorXd& uu1, VectorXd& uu2 );
  double Cal_eval( VectorXd u1, VectorXd u2 );
  void Check_eval( VectorXd u1,  VectorXd u2, double eval ) override;
  void Scale_Rotate_UU1_Procrustes( VectorXd u1, VectorXd& uu1 );
  void Scale_Rotate_UU2_Procrustes( VectorXd u2, VectorXd& uu2 );
  void SetK1() override;
  void SetK2() override;
  void Set_RotationMatrices1( VectorXd uu );
  void Set_RotationMatrices2( VectorXd uu );
  void Scale_Rotate_UU1( VectorXd& uu1 );
  void Scale_Rotate_UU2( VectorXd& uu1 );

  VectorXd fWW1, fWW2;     
  VectorXd fWW1_d, fWW2_d; 
  
  VectorXd fR1_cos, fR2_cos; 
  VectorXd fR1_sin, fR2_sin; 
  MatrixXd *fRot1, *fRot2; 
  MatrixXd *fRotD1, *fRotD2;

  MatrixXd fGG1, fGG2;   // sparseにすべきかも
  double fEval1_wGGw, fEval2_wGGw;
  SparseMatrix<double> fAA1;
  SparseMatrix<double> fAA2;
  VectorXd fTTWW1, fTTWW2;
  VectorXd fTTWW1_d, fTTWW2_d;
};

#endif

