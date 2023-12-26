/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_I__
#define __OPERATOR_I__

#ifndef  __OPERATOR_BASE__
#include "operator_base.h"
#endif

class TOperator_I : public TOperator_base {
 public:
  TOperator_I( int N1, int N1_in, int N2, int N2_in );
  ~TOperator_I(); 
  virtual void SetParameter(); // set default parameters
  virtual void SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2, VectorXd w2_in, MatrixXi kk2i ); // set information of the goal mesh
  void Set_values( int na, int reverse, double eval_best );
  virtual void Cal_u_SGD();
  virtual void Cal_u_SGD( int st0, int st1, int st2 );  
  virtual double SGD_main( VectorXd& u1, VectorXd& u2, VectorXd& u, VectorXd& ua );
  virtual void Check_eval( VectorXd u1,  VectorXd u2, double eval );
  double Iterative_Euclid_main( VectorXd& u1, VectorXd& u2, VectorXd& u, VectorXd& ua ); // III
  double Cal_eval_gzai( VectorXd& gzai, VectorXd& u1, VectorXd& u2 ); // III
  void Construct_B1_B2( int s );
  void Set_B1j( int st );
  void Set_B2j( int st );
  virtual void SetK1(); // Set fKK, fG, fH, fG2_inv
  virtual void SetK2(); // Set fKK, fG, fH, fG2_inv
  void Scale_Rotate_U1_Procrustes( VectorXd& u );
  void Scale_Rotate_U2_Procrustes( VectorXd& u );
  void Trans_UU1_center( VectorXd& uu );
  void Trans_UU2_center( VectorXd& uu );
  bool Check_Intersection( VectorXd& u1, VectorXd& u2 );

  //  virtual void Update_UU_top( double eval, VectorXd& u, VectorXd& ua, VectorXd& uu1, VectorXd& uu2,int st0, int st1, int st2 );
  //  virtual double Return_Eval_LS();
  //  virtual void Update_UU_LS( double eval, VectorXd& u, VectorXd& ua, VectorXd& uu1, VectorXd& uu2,int st0, int st1, int st2 );
  //  void Store_tiles_for_analize( double eval, VectorXd& gzai ); // analize
  
  int fMMd, fMMs;  // dens,sparse部分のサイズ
  MatrixXd fB1d; 
  SparseMatrix<double> fB1s; 
  MatrixXd fB2d; 
  SparseMatrix<double> fB2s;
  MatrixXd fB1dj; 
  SparseMatrix<double> fB1sj; 
  MatrixXd fB2dj; 
  SparseMatrix<double> fB2sj; 
  
  SparseMatrix<double> fG1, fG2;   // (G_I)
  int fMk;                   // diminsion of the parameters (=fBdj.cols()+fBsj.cols())
  int fMMk;                  // diminsion of the parameters (B1, B2で共通)

  double fAlpha; // A parameter for the weight setting (alpha_i for the inner points)

  MatrixXd fKK1, fKK2;   // Laplacian matrix with weights of the goal mesh (K)
  MatrixXi fKK1i, fKK2i; // Laplacian matrix or adjacency matrix 
  MatrixXd fH1, fH2;   // (- G_2^{-1} * G_1)
  MatrixXd fG2_inv1, fG2_inv2; // (G_2^{-1})

  int fNa; // fNはoperation_base.hで定義
  int fN1, fN2;         // the number of points of the goal polygon (n)
  int fN1_in, fN2_in;      // the number of the inner points (n')
  int fNN1, fNN2;        // the number of all points (n+n')
  VectorXd fW1, fW2;    // the coordinates of the goal polygon (w)
  VectorXd fW1_d, fW2_d;  // (w_c)
  VectorXd fW1_in, fW2_in; // the coordinates of the inner points of the goal mesh (w')

  MatrixXd fA1, fA2; // Asymmetric Laplacian matrix with weights of the goal mesh 
  MatrixXd fDout1, fDin1;
  MatrixXd fDout2, fDin2;

  MatrixXi fNei1, fNei2;      // fNei(i,*) is a set of the neighbors of point i (N(i))
  VectorXi fNei1_Size, fNei2_Size; // fNei_Size(i) is the number of the neighbors of point i 
  double fEval2_wGw, fEval1_wGw;   // fW.transpose()*fG*fW

  vector<T> fTripletList1;
  vector<T> fTripletList2;

  MatrixXd fG1_d; 
  MatrixXd fG2_d; 

  int fReverse1, fReverse2;

  double fEval;
  double fEval_best; 
  VectorXd fU;
  VectorXd fUa;
  VectorXd fUU1;
  VectorXd fUU2;
  int fSt0, fSt1, fSt2;
	    
	    


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

