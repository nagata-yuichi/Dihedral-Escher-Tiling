/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_I__
#include "operator_I.h"
#endif


TOperator_I::TOperator_I( int N1, int N1_in, int N2, int N2_in ) : TOperator_base( N1+N2 )
{
  fN1 = N1;
  fN1_in = N1_in; 
  fNN1 = fN1 + fN1_in;
  fW1 = VectorXd::Zero(2*fN1);
  fW1_d = VectorXd::Zero(2*fN1);
  fH1 = MatrixXd::Zero(2*fN1_in, 2*fN1);
  fG1 = SparseMatrix<double>(2*fN1,2*fN1);
  fG1.reserve(8*fN1);

  fW1_in = VectorXd::Zero(2*fN1_in);
  fKK1i = MatrixXi::Zero(fNN1,fNN1);
  fKK1 = MatrixXd::Zero(fNN1,fNN1); 
  fA1 = MatrixXd::Zero(fNN1,fNN1);
  fG2_inv1 = MatrixXd::Zero(2*fN1_in, 2*fN1_in);
  fG1_d = MatrixXd::Zero(2*fN1,2*fN1);
  
  fB1d = MatrixXd(2*fN1,2*fN1);
  fB1s = SparseMatrix<double>(2*fN1,2*fN1);
  fB1s.reserve(8*fN1);
  fB1dj = MatrixXd(2*fN1,2*fN1);
  fB1sj = SparseMatrix<double>(2*fN1,2*fN1);
  fB1sj.reserve(8*fN1);
  fNei1 = MatrixXi::Zero(fNN1,20);
  fNei1_Size = VectorXi::Zero(fNN1);
  fTripletList1.reserve(8*fN1);

  fN2 = N2;
  fN2_in = N2_in; 
  fNN2 = fN2 + fN2_in;
  fW2 = VectorXd::Zero(2*fN2);
  fW2_d = VectorXd::Zero(2*fN2);
  fH2 = MatrixXd::Zero(2*fN2_in, 2*fN2);
  fG2 = SparseMatrix<double>(2*fN2,2*fN2);
  fG2.reserve(8*fN2);

  fW2_in = VectorXd::Zero(2*fN2_in);
  fKK2i = MatrixXi::Zero(fNN2,fNN2);
  fKK2 = MatrixXd::Zero(fNN2,fNN2); 
  fA2 = MatrixXd::Zero(fNN2,fNN2);
  fG2_inv2 = MatrixXd::Zero(2*fN2_in, 2*fN2_in);
  fG2_d = MatrixXd::Zero(2*fN2,2*fN2);

  fB2d = MatrixXd(2*fN2,2*fN2);
  fB2s = SparseMatrix<double>(2*fN2,2*fN2);
  fB2s.reserve(8*fN2);
  fB2dj = MatrixXd(2*fN2,2*fN2);
  fB2sj = SparseMatrix<double>(2*fN2,2*fN2);
  fB2sj.reserve(8*fN2);
  fNei2 = MatrixXi::Zero(fNN2,20);
  fNei2_Size = VectorXi::Zero(fNN2);
  fTripletList2.reserve(8*fN2);


  fU = VectorXd::Zero(2*(fN1+fN2));
  fUa = VectorXd::Zero(2*(fN1+fN2));
  fUU1 = VectorXd::Zero(2*fNN1);
  fUU2 = VectorXd::Zero(2*fNN2);  
} 

TOperator_I::~TOperator_I() 
{

} 

void TOperator_I::SetParameter() 
{
  fAlpha = 1.0; // default (set as an input parameter)
}

void TOperator_I::SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2, VectorXd w2_in, MatrixXi kk2i ) 
{
  fW1 = w1;
  fW2 = w2;
  for( int i = 0; i < fN1; ++i ){
    fW1_d(ix1(i)) = w1(iy1(i));
    fW1_d(iy1(i)) = -w1(ix1(i));
  }
  for( int i = 0; i < fN2; ++i ){
    fW2_d(ix2(i)) = w2(iy2(i));
    fW2_d(iy2(i)) = -w2(ix2(i));
  }
  
  fW1_in = w1_in;
  fW2_in = w2_in;

  fKK1i = kk1i;  
  fKK2i = kk2i;

  this->SetK1();
  this->SetK2();
}

void TOperator_I::Set_values( int na, int reverse, double eval_best )
{
  fNa = na;
  fN = fN1+fN2-2*na-2;
  fReverse2 = reverse;
  fEval_best = eval_best;
}

void TOperator_I::Cal_u_SGD() 
{
  double eval; 
  VectorXd u1(2*fN1);
  VectorXd u2(2*fN2);
  VectorXd u(2*fN);     // points of the template
  VectorXd ua(2*fNa+4); // points of the ua (including both ends)

  VectorXd u1_in(2*fN1_in); 
  VectorXd uu1(2*fNN1); 
  VectorXd u2_in(2*fN2_in);
  VectorXd uu2(2*fNN2);

 
  ///////////////////////////////////////////////////////
  
  for( int st0 = 0; st0 < fN; ++st0 ){
    this->Construct_B1_B2( st0 ); 
    
    for( int st1 = 0; st1 < fN1; ++st1 ){ 
      this->Set_B1j( st1 );
      for( int st2 = 0; st2 < fN2 ; ++st2 ){  
	this->Set_B2j( st2 );

	eval = this->SGD_main( u1, u2, u, ua );

	fEval = 999999999.9;
	if( eval < fEval_best ){
	  if( Check_Intersection( u, ua ) == false ){
	    this->Scale_Rotate_U1_Procrustes( u1 ); 
	    this->Scale_Rotate_U2_Procrustes( u2 );

	    u1_in = fH1*(u1-fW1) + fW1_in; 
	    uu1.head(2*fN1) = u1;        
	    uu1.tail(2*fN1_in) = u1_in;  
	    u2_in = fH2*(u2-fW2) + fW2_in; 
	    uu2.head(2*fN2) = u2;        
	    uu2.tail(2*fN2_in) = u2_in;  

	    this->Trans_UU1_center( uu1 );
	    this->Trans_UU2_center( uu2 );

	    fEval = eval;
	    fU.head(2*fN) = u.head(2*fN);
	    fUa.head(2*fNa+4) = ua.head(2*fNa+4);
	    fUU1 = uu1;
	    fUU2 = uu2;
	    fSt0 = st0;
	    fSt1 = st1;
	    fSt2 = st2;
	  }
	}
      }
    }
  }
}

void TOperator_I::Cal_u_SGD( int st0, int st1, int st2 ) 
{
  double eval; 
  VectorXd u1(2*fN1);
  VectorXd u2(2*fN2);
  VectorXd u(2*fN);     // points of the template
  VectorXd ua(2*fNa+4); // points of the ua (including both ends)

  VectorXd u1_in(2*fN1_in); 
  VectorXd uu1(2*fNN1); 
  VectorXd u2_in(2*fN2_in);
  VectorXd uu2(2*fNN2);

 
  ///////////////////////////////////////////////////////

  this->Construct_B1_B2( st0 ); 
  this->Set_B1j( st1 );
  this->Set_B2j( st2 );

  eval = this->SGD_main( u1, u2, u, ua );

  fEval = 999999999.9;
  if( eval < fEval_best ){
    if( Check_Intersection( u, ua ) == false ){ 
      // this->Scale_Rotate_U1_Procrustes( u1 ); 
      // this->Scale_Rotate_U2_Procrustes( u2 );


      // とりあえず。後で整理したい。
      double tr, rt, cos1_th, sin1_th, cos2_th, sin2_th;
      VectorXd ur1(2*fN1);
      VectorXd ur2(2*fN2);
      double ugu;

      ugu = u1.transpose()*fG1*u1;
      tr = u1.transpose()*fG1*fW1; 
      rt = u1.transpose()*fG1*fW1_d;
      cos1_th = tr/ugu; // sqrt(tr*tr + rt*rt);
      sin1_th = rt/ugu; // sqrt(tr*tr + rt*rt);
      for( int k = 0; k < fN1; ++k ){
	ur1(ix1(k)) = cos1_th*u1(ix1(k)) - sin1_th*u1(iy1(k)); 
	ur1(iy1(k)) = sin1_th*u1(ix1(k)) + cos1_th*u1(iy1(k)); 
      }
      u1 = ur1;

      ugu = u2.transpose()*fG2*u2;
      tr = u2.transpose()*fG2*fW2; 
      rt = u2.transpose()*fG2*fW2_d;
      cos2_th = tr/ugu; // sqrt(tr*tr + rt*rt);
      sin2_th = rt/ugu; // sqrt(tr*tr + rt*rt);
      for( int k = 0; k < fN2; ++k ){
	ur2(ix2(k)) = cos2_th*u2(ix2(k)) - sin2_th*u2(iy2(k)); 
	ur2(iy2(k)) = sin2_th*u2(ix2(k)) + cos2_th*u2(iy2(k)); 
      }
      u2 = ur2;


      u1_in = fH1*(u1-fW1) + fW1_in; 
      uu1.head(2*fN1) = u1;        
      uu1.tail(2*fN1_in) = u1_in;  
      u2_in = fH2*(u2-fW2) + fW2_in; 
      uu2.head(2*fN2) = u2;        
      uu2.tail(2*fN2_in) = u2_in;


      // mesh representation of the tiles correspondint to u and ua
      double d1, d2;
      VectorXd uu1_d(2*fNN1), uu2_d(2*fNN2);
      d1 = cos1_th*cos1_th + sin1_th*sin1_th;
      d2 = cos2_th*cos2_th + sin2_th*sin2_th;

      for( int k = 0; k < fNN1; ++k ){
	uu1_d(iix1(k)) = 1.0/d1*(cos1_th*uu1(iix1(k)) + sin1_th*uu1(iiy1(k))); 
	uu1_d(iiy1(k)) = 1.0/d1*(-sin1_th*uu1(iix1(k)) + cos1_th*uu1(iiy1(k))); 
      }
      for( int k = 0; k < fNN2; ++k ){
	uu2_d(iix2(k)) = 1.0/d2*(cos2_th*uu2(iix2(k)) + sin2_th*uu2(iiy2(k))); 
	uu2_d(iiy2(k)) = 1.0/d2*(-sin2_th*uu2(iix2(k)) + cos2_th*uu2(iiy2(k))); 
      }
      if( fReverse2 == 1 )
	for( int k = 0; k < fNN2; ++k )
	  uu2_d(iix2(k)) *= -1.0;

      // uu1 = uu1_d;
      // uu2 = uu2_d;
      
      fEval = eval;
      fU.head(2*fN) = u.head(2*fN);
      fUa.head(2*fNa+4) = ua.head(2*fNa+4);
      fUU1 = uu1_d;
      fUU2 = uu2_d;
      fSt0 = st0;
      fSt1 = st1;
      fSt2 = st2;
    }
  }
}

double TOperator_I::SGD_main( VectorXd& u1, VectorXd& u2, VectorXd& u, VectorXd& ua )
{
  VectorXd gzai(fMMd+fMMs);
  VectorXd ag1(fMMd+fMMs); // A*gzai
  VectorXd cg1(fMMd+fMMs); // C*gzai 
  VectorXd ag2(fMMd+fMMs); // A*gzai
  VectorXd cg2(fMMd+fMMs); // C*gzai
  VectorXd av1(fMMd+fMMs); // A*gzai
  VectorXd cv1(fMMd+fMMs); // C*gzai 
  VectorXd av2(fMMd+fMMs); // A*gzai
  VectorXd cv2(fMMd+fMMs); // C*gzai

  VectorXd gzai_curr(fMMd+fMMs);
  VectorXd ag1_curr(fMMd+fMMs);
  VectorXd cg1_curr(fMMd+fMMs); 
  VectorXd ag2_curr(fMMd+fMMs); 
  VectorXd cg2_curr(fMMd+fMMs);
  double gcg1_curr, gag1_curr, gcg2_curr, gag2_curr;
  double eval, eval1, eval2, eval_curr;
  double alpha, alpha_best, alpha_init, alpha_min, alpha_max;

  VectorXd grad(fMMd+fMMs);  
  VectorXd grad1(fMMd+fMMs);  
  VectorXd grad2(fMMd+fMMs);  
  VectorXd grad_curr(fMMd+fMMs);
  VectorXd grad_pre(fMMd+fMMs);
  VectorXd grad_v(fMMd+fMMs);
  VectorXd tmp1(2*fN1);  
  VectorXd tmp2(2*fN2);  
  VectorXd p1(fMMd+fMMs);
  VectorXd p1_d(fMMd+fMMs);
  VectorXd p2(fMMd+fMMs);
  VectorXd p2_d(fMMd+fMMs);
  
  VectorXd ud1(2*fN1);
  VectorXd ud2(2*fN2);
  double ratio = (double)fN2/(double)fN1;

  static VectorXd gzai_previous = VectorXd::Random(2*(fN1+fN2));
  static int MMd_previous = 2*(fN1+fN2);
  static int MMs_previous = 2*(fN1+fN2);

  
  if( fMMd == MMd_previous && fMMs == MMs_previous ) 
    gzai_curr = gzai_previous.head(fMMd+fMMs); // option (useful in LS)
  else
    gzai_curr = VectorXd::Random(fMMd+fMMs);
  gzai_curr = gzai_curr.normalized();

  // p1, p1_d 
  tmp1 = fG1*fW1; 
  p1.head(fMMd) = fB1dj.transpose()*tmp1;
  p1.tail(fMMs) = fB1sj.transpose()*tmp1;
  tmp1 = fG1*fW1_d; 
  p1_d.head(fMMd) = fB1dj.transpose()*tmp1;
  p1_d.tail(fMMs) = fB1sj.transpose()*tmp1;

  // p2, p2_d 
  tmp2 = fG2*fW2; 
  p2.head(fMMd) = fB2dj.transpose()*tmp2;
  p2.tail(fMMs) = fB2sj.transpose()*tmp2;
  tmp2 = fG2*fW2_d;
		      
  p2_d.head(fMMd) = fB2dj.transpose()*tmp2;
  p2_d.tail(fMMs) = fB2sj.transpose()*tmp2;

  // cg1
  u1 = fB1dj*gzai_curr.head(fMMd); 
  u1 += fB1sj*gzai_curr.tail(fMMs);
  tmp1 = fG1*u1;
  cg1_curr.head(fMMd) = fB1dj.transpose()*tmp1;
  cg1_curr.tail(fMMs) = fB1sj.transpose()*tmp1;
  // ag1
  ag1_curr = p1.dot(gzai_curr)*p1 + p1_d.dot(gzai_curr)*p1_d;
  // gcg1, gag1
  gcg1_curr = cg1_curr.dot(gzai_curr);
  gag1_curr = ag1_curr.dot(gzai_curr);

  // cg2
  u2 = fB2dj*gzai_curr.head(fMMd); 
  u2 += fB2sj*gzai_curr.tail(fMMs);
  tmp2 = fG2*u2;
  cg2_curr.head(fMMd) = fB2dj.transpose()*tmp2;
  cg2_curr.tail(fMMs) = fB2sj.transpose()*tmp2;
  // ag2
  ag2_curr = p2.dot(gzai_curr)*p2 + p2_d.dot(gzai_curr)*p2_d;
  // gcg2, gag2
  gcg2_curr = cg2_curr.dot(gzai_curr);
  gag2_curr = ag2_curr.dot(gzai_curr);

  // grad
  grad1 = gcg1_curr*ag1_curr - gag1_curr*cg1_curr;
  grad1 *= -1.0/(gcg1_curr*gcg1_curr*fEval1_wGw);
  grad2 = gcg2_curr*ag2_curr - gag2_curr*cg2_curr;
  grad2 *= -1.0/(gcg2_curr*gcg2_curr*fEval2_wGw);
  grad_curr = grad1 + ratio*grad2;
  grad_v = -grad_curr;
  

  eval_curr = 0.0;
  eval_curr += 1.0 - gag1_curr/(gcg1_curr*fEval1_wGw);
  eval_curr += ratio*(1.0 - gag2_curr/(gcg2_curr*fEval2_wGw));

  alpha_init = 1.0;
  for( int iter = 0; iter < 1000; ++iter ){
    // cv1
    ud1 = fB1dj*grad_v.head(fMMd); 
    ud1 += fB1sj*grad_v.tail(fMMs);
    tmp1 = fG1*ud1;
    cv1.head(fMMd) = fB1dj.transpose()*tmp1;
    cv1.tail(fMMs) = fB1sj.transpose()*tmp1;
    // av1
    av1 = p1.dot(grad_v)*p1 + p1_d.dot(grad_v)*p1_d;
    double gcv1 = cv1.dot(gzai_curr);
    double gav1 = av1.dot(gzai_curr);
    double vcv1 = cv1.dot(grad_v);
    double vav1 = av1.dot(grad_v);

    // cv2
    ud2 = fB2dj*grad_v.head(fMMd); 
    ud2 += fB2sj*grad_v.tail(fMMs);
    tmp2 = fG2*ud2;
    cv2.head(fMMd) = fB2dj.transpose()*tmp2;
    cv2.tail(fMMs) = fB2sj.transpose()*tmp2;
    // av2
    av2 = p2.dot(grad_v)*p2 + p2_d.dot(grad_v)*p2_d;
    double gcv2 = cv2.dot(gzai_curr);
    double gav2 = av2.dot(gzai_curr);
    double vcv2 = cv2.dot(grad_v);
    double vav2 = av2.dot(grad_v);

    int flag_local; // 局所解含むか: 0:no, 1:yes 
    flag_local = 0;
    alpha_min = 0.0;
    alpha_max = alpha_init;

    for( int t = 0; t < 5; ++t ){
      double eval_best = 99999999.9;
      int k_best;
      for( int k = 0; k <= 10; ++k ){
	alpha = alpha_min + 0.1*(double)k*(alpha_max-alpha_min);
	double gag1 = gag1_curr + 2*alpha*gav1 + alpha*alpha*vav1;
	double gcg1 = gcg1_curr + 2*alpha*gcv1 + alpha*alpha*vcv1;
	double gag2 = gag2_curr + 2*alpha*gav2 + alpha*alpha*vav2;
	double gcg2 = gcg2_curr + 2*alpha*gcv2 + alpha*alpha*vcv2;

	eval1 = 1.0 - gag1/(gcg1*fEval1_wGw);
	eval2 = 1.0 - gag2/(gcg2*fEval2_wGw);
	eval = eval1 + ratio*eval2;
	assert( eval1 > 0.0 );
	assert( eval2 > 0.0 );
	if( eval < eval_best ){
	  eval_best = eval;
	  alpha_best = alpha;
	  k_best = k;
	}
	// printf( "t = %d k = %d alpha = %lf eval = %lf \n", t, k, alpha, eval );
      }

      if( flag_local == 0 ){
	if( k_best == 0 )
	  alpha_max *= 0.1;
	else if( k_best == 10 )
	  alpha_max *= 10.0;
	else{
	  double alpha_min_temp = alpha_min + 0.1*(double)(k_best-1)*(alpha_max-alpha_min);
	  double alpha_max_temp = alpha_min + 0.1*(double)(k_best+1)*(alpha_max-alpha_min);
	  alpha_min = alpha_min_temp;
	  alpha_max = alpha_max_temp;
	  flag_local = 1;
	}
      }
      else if( flag_local == 1 ){
	if( k_best == 0 )
	  alpha_max = alpha_min + 0.1*(alpha_max-alpha_min);
	else if( k_best == 10 )
	  alpha_min = alpha_min + 0.9*(alpha_max-alpha_min);
	else{
	  double alpha_min_temp = alpha_min + 0.1*(double)(k_best-1)*(alpha_max-alpha_min);
	  double alpha_max_temp = alpha_min + 0.1*(double)(k_best+1)*(alpha_max-alpha_min);
	  alpha_min = alpha_min_temp;
	  alpha_max = alpha_max_temp;
	}
      }
    }

    alpha = alpha_best;
    alpha_init = 0.5*alpha_best + 0.5*alpha_init;
    double gag1 = gag1_curr + 2*alpha*gav1 + alpha*alpha*vav1;
    double gcg1 = gcg1_curr + 2*alpha*gcv1 + alpha*alpha*vcv1;
    double gag2 = gag2_curr + 2*alpha*gav2 + alpha*alpha*vav2;
    double gcg2 = gcg2_curr + 2*alpha*gcv2 + alpha*alpha*vcv2;

    eval1 = 1.0 - gag1/(gcg1*fEval1_wGw);
    eval2 = 1.0 - gag2/(gcg2*fEval2_wGw);
    eval = eval1 + ratio*eval2;
    assert( eval1 > 0.0 );
    assert( eval2 > 0.0 );
    
    // printf( "i=%d (%lf): %lf -> %lf \n", iter, alpha, eval_curr, eval );
    if( eval_curr - eval < 0.00001 )
      break;
    
    if( eval < eval_curr ){
      eval_curr = eval;
      gzai_curr += alpha*grad_v;
      double norm_before = gzai_curr.norm();
      gzai_curr = gzai_curr.normalized();
      double norm_after = gzai_curr.norm();
      double r = norm_after/norm_before;

      gcg1_curr = r*r*gcg1;
      gag1_curr = r*r*gag1;
      gcg2_curr = r*r*gcg2;
      gag2_curr = r*r*gag2;

      cg1_curr += alpha*cv1;
      ag1_curr += alpha*av1;
      cg2_curr += alpha*cv2;
      ag2_curr += alpha*av2;
      cg1_curr *= r;
      ag1_curr *= r;
      cg2_curr *= r;
      ag2_curr *= r;

      grad_pre = grad_curr;
      grad1 = gcg1_curr*ag1_curr - gag1_curr*cg1_curr;
      grad1 *= -1.0/(gcg1_curr*gcg1_curr*fEval1_wGw);
      grad2 = gcg2_curr*ag2_curr - gag2_curr*cg2_curr;
      grad2 *= -1.0/(gcg2_curr*gcg2_curr*fEval2_wGw);
      grad_curr = grad1 + ratio*grad2;

      // double alpha_c = pow(grad_curr.norm(),2)/pow(grad_pre.norm(),2);
      // double alpha_c = grad_curr.dot(grad_curr-grad_pre)/pow(grad_pre.norm(),2);
      double alpha_c = grad_curr.dot(grad_curr-grad_pre)/grad_v.dot(grad_curr-grad_pre);

      assert( (grad_curr-grad_pre).norm() > 0.0000001 ); 

      if( flag_local == 1 ){ 
	grad_v *= alpha_c;
	grad_v -= grad_curr;
      }
      else 
	grad_v = -grad_curr;  // これはあったほうが良い

      grad_v -= grad_v.dot(gzai_curr)*gzai_curr; // あってもなくても同じかも

      // this->Store_tiles_for_analize( eval_curr, gzai_curr ); // analize
    }
  }
  gzai = gzai_curr;
  eval = eval_curr;
  assert( eval < 9999999999.9 );

  gzai = gzai.normalized(); 
  // u1, u2はW1,W2のインデックスに対応している
  u1 = fB1dj*gzai.head(fMMd); 
  u1 += fB1sj*gzai.tail(fMMs);
  u2 = fB2dj*gzai.head(fMMd); 
  u2 += fB2sj*gzai.tail(fMMs);

  u = fBd*gzai.head(fMd) + fBs*gzai.segment(fMd,fMs);

  VectorXd tmp(2*fN1);
  tmp = fB1d*gzai.head(fMMd) + fB1s*gzai.segment(fMMd,fMMs);
  ua.head(2*fNa+2) = tmp.tail(2*fNa+2);
  ua(2*(fNa+1)) = tmp(0);
  ua(2*(fNa+1)+1) = tmp(1);

  // this->Check_eval( u1, u2, eval ); // for check

  gzai_previous.head(fMMd+fMMs) = gzai;
  MMd_previous = fMMd;
  MMs_previous = fMMs;


  

  return eval;
}

void TOperator_I::Check_eval( VectorXd u1, VectorXd u2, double eval )
{
  VectorXd uu1(2*fNN1), uu2(2*fNN2);
  VectorXd u1_in(2*fN1_in), u2_in(2*fN2_in);
  VectorXd ww1(2*fN1+2*fN1_in), ww2(2*fN2+2*fN2_in);
  VectorXd v(2); 
  double sum1, sum2, sum;
  double ratio = (double)fN2/(double)fN1;

  ww1.head(2*fN1) = fW1;
  ww2.head(2*fN2) = fW2;
  ww1.tail(2*fN1_in) = fW1_in;
  ww2.tail(2*fN2_in) = fW2_in;

  this->Scale_Rotate_U1_Procrustes( u1 ); 
  this->Scale_Rotate_U2_Procrustes( u2 );
  
  u1_in = fH1*(u1-fW1) + fW1_in; 
  uu1.head(2*fN1) = u1;        
  uu1.tail(2*fN1_in) = u1_in;  
  u2_in = fH2*(u2-fW2) + fW2_in; 
  uu2.head(2*fN2) = u2;        
  uu2.tail(2*fN2_in) = u2_in;  
  
  sum1 = 0.0;
  for( int i = 0; i < fNN1; ++i ){
    int size = fNei1_Size(i);

    for( int s = 0; s < size; ++s ){
      int j = fNei1(i,s); 
      double aij = fA1(i,j);
      assert( i != j ); 
      v(0) = (uu1(iix1(j)) - uu1(iix1(i))) - (ww1(iix1(j)) - ww1(iix1(i)));
      v(1) = (uu1(iiy1(j)) - uu1(iiy1(i))) - (ww1(iiy1(j)) - ww1(iiy1(i)));
      sum1 += aij*(v(0)*v(0)+v(1)*v(1));
    }
  }
  
  sum2 = 0.0;
  for( int i = 0; i < fNN2; ++i ){
    int size = fNei2_Size(i);
    
    for( int s = 0; s < size; ++s ){
      int j = fNei2(i,s); 
      double aij = fA2(i,j);
      assert( i != j ); 
      v(0) = (uu2(iix2(j)) - uu2(iix2(i))) - (ww2(iix2(j)) - ww2(iix2(i)));
      v(1) = (uu2(iiy2(j)) - uu2(iiy2(i))) - (ww2(iiy2(j)) - ww2(iiy2(i)));
      sum2 += aij*(v(0)*v(0)+v(1)*v(1));
    }
  }
	
  // v(0) = uu1(iix1(0)) - ww1(iix1(0));
  // v(1) = uu1(iiy1(0)) - ww1(iiy1(0));
  // sum1 += v(0)*v(0)+v(1)*v(1);
  // v(0) = uu2(iix2(0)) - ww2(iix2(0));
  // v(1) = uu2(iiy2(0)) - ww2(iiy2(0));
  // sum2 += v(0)*v(0)+v(1)*v(1);

  sum = sum1/fEval1_wGw + ratio*sum2/fEval2_wGw;
  
  printf( "%lf %lf -> %lf\n", sum , eval, eval-sum );
  assert( fabs(eval-sum) < 0.000001 ); 
}


void TOperator_I::Construct_B1_B2( int s )
{
  int i, j;
  int i_after; 
  double val;
  VectorXd tmp(2*fN2);
  VectorXd b1(2*fN); // ith row がB1
  VectorXd b2(2*fN); // ith row がB2
  int s1, s2;

  fMMd = fMd; 
  fMMs = fMs+2*fNa;
  fMMk = fMMd+fMMs;

  fB1d.resize(2*fN1,fMMd); 
  fB1s.resize(2*fN1,fMMs);
  fB2d.resize(2*fN2,fMMd); 
  fB2s.resize(2*fN2,fMMs);

  fB1dj.resize(2*fN1,fMMd); 
  fB1sj.resize(2*fN1,fMMs);
  fB2dj.resize(2*fN2,fMMd); 
  fB2sj.resize(2*fN2,fMMs);
  

  s1 = s;
  s2 = s+fN1-fNa-1;

  for( int k = 0; k < 2*fN; ++k ){
    b1(k) = -1;
    b2(k) = -1;
  }
  for( int k = 0; k < fN1-fNa; ++k ){
    b1(ix(s1+k)) = ix1(k);
    b1(iy(s1+k)) = iy1(k);
  }

  for( int k = 0; k < fN2-fNa; ++k ){
    b2(ix(s2+k)) = ix2(k);
    b2(iy(s2+k)) = iy2(k);
  }

  fB1d = MatrixXd::Zero(2*fN1,fMMd);
  for( int k = 0; k < fN1-fNa; ++k ){
    fB1d.row(ix1(k)) = fBd.row(ix(s1+k));
    fB1d.row(iy1(k)) = fBd.row(iy(s1+k));
  } 
  fB2d = MatrixXd::Zero(2*fN2,fMMd);
  for( int k = 0; k < fN2-fNa; ++k ){
    fB2d.row(ix2(k)) = fBd.row(ix(s2+k));
    fB2d.row(iy2(k)) = fBd.row(iy(s2+k));
  } 


  fTripletList1.clear();
  fTripletList2.clear();
  for ( int k = 0 ; k < fMs; ++k ) { 
    for ( SparseMatrix<double>::InnerIterator it(fBs,k); it; ++it ){
      val = it.value();
      i = (int)it.row();
      j = (int)it.col();
      assert( j == k ); 

      if( b1(i) != -1 ) 
	fTripletList1.push_back(T(b1(i),j,val));
      if( b2(i) != -1 ) 
	fTripletList2.push_back(T(b2(i),j,val));
    }
  }
  for ( int k = 0 ; k < fNa; ++k ){  
    fTripletList1.push_back(T(2*(fN1-fNa+k),fMs+2*k,1.0));
    fTripletList1.push_back(T(2*(fN1-fNa+k)+1,fMs+2*k+1,1.0));
    fTripletList2.push_back(T(2*(fN2-1-k),fMs+2*k,1.0));
    fTripletList2.push_back(T(2*(fN2-1-k)+1,fMs+2*k+1,1.0));
  }
  fB1s.setFromTriplets(fTripletList1.begin(), fTripletList1.end());
  fB2s.setFromTriplets(fTripletList2.begin(), fTripletList2.end());
}


void TOperator_I::Set_B1j( int st )
{
  if( fReverse1 == 0 ){
    for( int i = 0; i < fN1; ++i ){
      fB1dj.row(ix1(i)) = fB1d.row(ix1(i-st));
      fB1dj.row(iy1(i)) = fB1d.row(iy1(i-st));
    }

    fTripletList1.clear();
    for ( int k = 0 ; k < fB1s.outerSize(); ++k ) { 
      for ( SparseMatrix<double>::InnerIterator it(fB1s,k); it; ++it ){
	double val = it.value();
	int i = (int)it.row();
	int j = (int)it.col();
	fTripletList1.push_back(T((i+2*st)%(2*fN1),j,val));
      }
    }
    fB1sj.setFromTriplets(fTripletList1.begin(), fTripletList1.end());

  }
  else if( fReverse1 == 1 ){
    for( int i = 0; i < fN1; ++i ){
      fB1dj.row(ix1(i)) = - fB1d.row(ix1(-i-st+fN1));
      fB1dj.row(iy1(i)) = fB1d.row(iy1(-i-st+fN1));
    }

    fTripletList1.clear();
    for ( int k = 0 ; k < fB1s.outerSize(); ++k ) { 
      for ( SparseMatrix<double>::InnerIterator it(fB1s,k); it; ++it ){
	double val = it.value();
	int i = (int)it.row();
	int j = (int)it.col();
	int ix = i/2;
	if( i % 2 == 0 ){
	  val *= -1.0;
	  fTripletList1.push_back(T((-ix*2-2*st + 4*fN1)%(2*fN1),j,val));
	}
	else{
	  fTripletList1.push_back(T((-ix*2-2*st+1 + 4*fN1)%(2*fN1),j,val));
	}
      }
    }
    fB1sj.setFromTriplets(fTripletList1.begin(), fTripletList1.end());
  }
}

void TOperator_I::Set_B2j( int st )
{
  if( fReverse2 == 0 ){
    for( int i = 0; i < fN2; ++i ){
      fB2dj.row(ix2(i)) = fB2d.row(ix2(i-st));
      fB2dj.row(iy2(i)) = fB2d.row(iy2(i-st));
    }

    fTripletList2.clear();
    for ( int k = 0 ; k < fB2s.outerSize(); ++k ) { 
      for ( SparseMatrix<double>::InnerIterator it(fB2s,k); it; ++it ){
	double val = it.value();
	int i = (int)it.row();
	int j = (int)it.col();
	fTripletList2.push_back(T((i+2*st)%(2*fN2),j,val));
      }
    }
    fB2sj.setFromTriplets(fTripletList2.begin(), fTripletList2.end());

  }
  else if( fReverse2 == 1 ){
    for( int i = 0; i < fN2; ++i ){
      fB2dj.row(ix2(i)) = - fB2d.row(ix2(-i-st+fN2));
      fB2dj.row(iy2(i)) = fB2d.row(iy2(-i-st+fN2));
    }

    fTripletList2.clear();
    for ( int k = 0 ; k < fB2s.outerSize(); ++k ) { 
      for ( SparseMatrix<double>::InnerIterator it(fB2s,k); it; ++it ){
	double val = it.value();
	int i = (int)it.row();
	int j = (int)it.col();
	int ix = i/2;
	if( i % 2 == 0 ){
	  val *= -1.0;
	  fTripletList2.push_back(T((-ix*2-2*st + 4*fN2)%(2*fN2),j,val));
	}
	else{
	  fTripletList2.push_back(T((-ix*2-2*st+1 + 4*fN2)%(2*fN2),j,val));
	}
      }
    }
    fB2sj.setFromTriplets(fTripletList2.begin(), fTripletList2.end());
  }
}



void TOperator_I::SetK1() 
{
  /////// set fKK, etc  //////
  MatrixXd K0 = MatrixXd::Zero(fN1,fN1); 
  MatrixXd K1 = MatrixXd::Zero(fN1,fN1_in); 
  MatrixXd K2 = MatrixXd::Zero(fN1_in,fN1_in); 
  MatrixXd K2_inv = MatrixXd::Zero(fN1_in,fN1_in);
  MatrixXd H = MatrixXd::Zero(fN1_in,fN1); 

  for( int i = 0; i < fNN1; ++i ){
    int num = 0;
    for( int j = 0; j < fNN1; ++j ){
      if( fKK1i(i,j) != 0 && i != j )
	fNei1(i,num++) = j;
      assert( num < 20 );
    }
    fNei1_Size(i) = num;
  }

  double alpha;
  fA1 = MatrixXd::Zero(fNN1,fNN1); 
  for( int i = 0; i < fNN1; ++i ){
    for( int s = 0; s < fNei1_Size(i); ++ s ){
      int j = fNei1(i,s);
      if( i < fN1 ) 
	alpha = 1.0;
      else 
	alpha = fAlpha;

      fA1(i,j) = alpha;
    }
  }

  fDout1 = MatrixXd::Zero(fNN1,fNN1); 
  fDin1 = MatrixXd::Zero(fNN1,fNN1); 
  for( int i = 0; i < fNN1; ++i ){
    for( int s = 0; s < fNei1_Size(i); ++ s ){
      int j = fNei1(i,s);
      fDout1(i,i) += fA1(i,j);
      fDin1(i,i) += fA1(j,i);
    }
  }

  fKK1 = fDout1 + fDin1 - 2.0*fA1;
  for( int i = 0; i < fNN1; ++i ){
    for( int j = i+1; j < fNN1; ++j ){
      double tmp = 0.5*(fKK1(i,j)+fKK1(j,i));
      fKK1(i,j) = tmp;
      fKK1(j,i) = tmp;
    }
  }
  // fKK1(0,0) += 1.0; // 正則化 // III
  
  K0 = fKK1.topLeftCorner(fN1,fN1);
  K1 = fKK1.topRightCorner(fN1,fN1_in);
  K2 = fKK1.bottomRightCorner(fN1_in,fN1_in);
  K2_inv = K2.inverse();

  // check
  MatrixXd TT = MatrixXd::Zero(fN1_in,fN1_in); 
  TT = K2*K2_inv; 
  for( int i = 0; i < fN1_in; ++i ){
    for( int j = 0; j < fN1_in; ++j ){
      if( i == j ) assert( fabs(TT(i,j)-1.0) < 0.000001 );
      else assert( fabs(TT(i,j)) < 0.000001 );
    }
  }

  MatrixXd K = MatrixXd::Zero(fN1,fN1);
  K = K0 - K1*K2_inv*K1.transpose();

  // check
  for( int i = 0; i < fN1; ++i )
    for( int j = 0; j < fN1; ++j )
      assert( fabs(K(i,j)-K(j,i)) < 0.000001 );

  /////// set fG and fH //////
  fG1_d = MatrixXd::Zero(2*fN1,2*fN1);
  for( int i = 0; i < fN1; ++i ){
    for( int j = 0; j < fN1; ++j ){
      fG1_d(ix1(i),ix(j)) = K(i,j);
      fG1_d(iy1(i),iy(j)) = K(i,j);
    }
  }
  
  double val;
  fTripletList.clear();
  for( int j = 0; j < 2*fN1; ++j ){
    for( int i = 0; i < 2*fN1; ++i ){
      val = fG1_d(i,j); 
      if( fabs(val) > 0.000000001 ) // このくらい小さくしないと誤差が出る
	fTripletList.push_back(T(i,j,val)); //代入の順番によって速度が変化
    }
  }
  fG1.setFromTriplets(fTripletList.begin(), fTripletList.end());
  fEval1_wGw = fW1.transpose()*fG1*fW1;
  // fEval1_wGw -= fW1[0]*fW1[0]+fW1[1]*fW1[1];  // 正則化 // III

  // check
  SelfAdjointEigenSolver<MatrixXd> es(fG1_d,false);
  int  ne = es.eigenvalues().size();
  // std::cout << es.eigenvalues() << std::endl;
  assert( es.eigenvalues()(0) >= -0.0000001 );

  // H
  H = - K2_inv*K1.transpose();
  for( int i = 0; i < fN1_in; ++i ){
    for( int j = 0; j < fN1; ++j ){
      fH1(2*i,ix1(j)) = H(i,j); // %いらないのでは
      fH1(2*i+1,iy1(j)) = H(i,j);
    }
  }

  fG2_inv1 = MatrixXd::Zero(2*fN1_in, 2*fN1_in);
  for( int i = 0; i < fN1_in; ++i ){
    for( int j = 0; j < fN1_in; ++j ){
      fG2_inv1(2*i,2*j) = K2_inv(i,j);
      fG2_inv1(2*i+1,2*j+1) = K2_inv(i,j);
    }
  }
}

void TOperator_I::SetK2() 
{
  /////// set fKK, etc  //////
  MatrixXd K0 = MatrixXd::Zero(fN2,fN2); 
  MatrixXd K1 = MatrixXd::Zero(fN2,fN2_in); 
  MatrixXd K2 = MatrixXd::Zero(fN2_in,fN2_in); 
  MatrixXd K2_inv = MatrixXd::Zero(fN2_in,fN2_in);
  MatrixXd H = MatrixXd::Zero(fN2_in,fN2); 

  for( int i = 0; i < fNN2; ++i ){
    int num = 0;
    for( int j = 0; j < fNN2; ++j ){
      if( fKK2i(i,j) != 0 && i != j )
	fNei2(i,num++) = j;
      assert( num < 20 );
    }
    fNei2_Size(i) = num;
  }

  double alpha;
  fA2 = MatrixXd::Zero(fNN2,fNN2); 
  for( int i = 0; i < fNN2; ++i ){
    for( int s = 0; s < fNei2_Size(i); ++ s ){
      int j = fNei2(i,s);
      if( i < fN2 ) 
	alpha = 1.0;
      else 
	alpha = fAlpha;

      fA2(i,j) = alpha;
    }
  }

  fDout2 = MatrixXd::Zero(fNN2,fNN2); 
  fDin2 = MatrixXd::Zero(fNN2,fNN2); 
  for( int i = 0; i < fNN2; ++i ){
    for( int s = 0; s < fNei2_Size(i); ++ s ){
      int j = fNei2(i,s);
      fDout2(i,i) += fA2(i,j);
      fDin2(i,i) += fA2(j,i);
    }
  }

  fKK2 = fDout2 + fDin2 - 2.0*fA2;
  for( int i = 0; i < fNN2; ++i ){
    for( int j = i+1; j < fNN2; ++j ){
      double tmp = 0.5*(fKK2(i,j)+fKK2(j,i));
      fKK2(i,j) = tmp;
      fKK2(j,i) = tmp;
    }
  }
  // fKK2(0,0) += 1.0; // 正則化
  
  K0 = fKK2.topLeftCorner(fN2,fN2);
  K1 = fKK2.topRightCorner(fN2,fN2_in);
  K2 = fKK2.bottomRightCorner(fN2_in,fN2_in);
  K2_inv = K2.inverse();

  // check
  MatrixXd TT = MatrixXd::Zero(fN2_in,fN2_in); 
  TT = K2*K2_inv; 
  for( int i = 0; i < fN2_in; ++i ){
    for( int j = 0; j < fN2_in; ++j ){
      if( i == j ) assert( fabs(TT(i,j)-1.0) < 0.000001 );
      else assert( fabs(TT(i,j)) < 0.000001 );
    }
  }

  MatrixXd K = MatrixXd::Zero(fN2,fN2);
  K = K0 - K1*K2_inv*K1.transpose();
  // check
  for( int i = 0; i < fN2; ++i )
    for( int j = 0; j < fN2; ++j )
      assert( fabs(K(i,j)-K(j,i)) < 0.000001 );

  /////// set fG and fH //////
  fG2_d = MatrixXd::Zero(2*fN2,2*fN2);
  for( int i = 0; i < fN2; ++i ){
    for( int j = 0; j < fN2; ++j ){
      fG2_d(ix2(i),ix(j)) = K(i,j);
      fG2_d(iy2(i),iy(j)) = K(i,j);
    }
  }

  double val;
  fTripletList.clear();
  for( int j = 0; j < 2*fN2; ++j ){
    for( int i = 0; i < 2*fN2; ++i ){
      val = fG2_d(i,j); 
      if( fabs(val) > 0.000000001 ) // このくらい小さくしないと誤差が出る
	fTripletList.push_back(T(i,j,val)); //代入の順番によって速度が変化
    }
  }
  fG2.setFromTriplets(fTripletList.begin(), fTripletList.end());
  fEval2_wGw = fW2.transpose()*fG2*fW2;

  // check
  SelfAdjointEigenSolver<MatrixXd> es(fG2_d,false);
  int  ne = es.eigenvalues().size();
  // std::cout << es.eigenvalues() << std::endl;
  assert( es.eigenvalues()(0) >= -0.0000001 );

  // H
  H = - K2_inv*K1.transpose();
  for( int i = 0; i < fN2_in; ++i ){
    for( int j = 0; j < fN2; ++j ){
      fH2(2*i,ix2(j)) = H(i,j);
      fH2(2*i+1,iy2(j)) = H(i,j);
    }
  }

  fG2_inv2 = MatrixXd::Zero(2*fN2_in, 2*fN2_in);
  for( int i = 0; i < fN2_in; ++i ){
    for( int j = 0; j < fN2_in; ++j ){
      fG2_inv2(2*i,2*j) = K2_inv(i,j);
      fG2_inv2(2*i+1,2*j+1) = K2_inv(i,j);
    }
  }
}

void TOperator_I::Scale_Rotate_U1_Procrustes( VectorXd& u )
{
  double tr, rt, cos_th, sin_th;
  VectorXd ur(2*fN1);
  double ugu;
  assert( u.rows() == 2*fN1 );

  ugu = u.transpose()*fG1*u;
  tr = u.transpose()*fG1*fW1; 
  rt = u.transpose()*fG1*fW1_d;
  cos_th = tr/ugu; // sqrt(tr*tr + rt*rt);
  sin_th = rt/ugu; // sqrt(tr*tr + rt*rt);
  for( int k = 0; k < fN1; ++k ){
    ur(ix1(k)) = cos_th*u(ix1(k)) - sin_th*u(iy1(k)); 
    ur(iy1(k)) = sin_th*u(ix1(k)) + cos_th*u(iy1(k)); 
  }
  u = ur;
}

void TOperator_I::Scale_Rotate_U2_Procrustes( VectorXd& u )
{
  double tr, rt, cos_th, sin_th;
  VectorXd ur(2*fN2);
  double ugu;
  assert( u.rows() == 2*fN2 );

  ugu = u.transpose()*fG2*u;
  tr = u.transpose()*fG2*fW2; 
  rt = u.transpose()*fG2*fW2_d;
  cos_th = tr/ugu; // sqrt(tr*tr + rt*rt);
  sin_th = rt/ugu; // sqrt(tr*tr + rt*rt);
  for( int k = 0; k < fN2; ++k ){
    ur(ix2(k)) = cos_th*u(ix2(k)) - sin_th*u(iy2(k)); 
    ur(iy2(k)) = sin_th*u(ix2(k)) + cos_th*u(iy2(k)); 
  }
  u = ur;
}


void TOperator_I::Trans_UU1_center( VectorXd& uu )
{
  VectorXd centerU(2);
  VectorXd centerW(2);

  centerU = VectorXd::Zero(2); 
  centerW = VectorXd::Zero(2);

  for( int i = 0; i < fN1; ++i ){  
    centerW += fW1.segment(2*i,2);
    centerU += uu.segment(2*i,2);
  }
  centerW /= (double)fN1;
  centerU /= (double)fN1;

  for( int i = 0; i < fNN1; ++i )  
    uu.segment(2*i,2) += (-centerU + centerW);
}

void TOperator_I::Trans_UU2_center( VectorXd& uu )
{
  VectorXd centerU(2);
  VectorXd centerW(2);

  centerU = VectorXd::Zero(2); 
  centerW = VectorXd::Zero(2);

  for( int i = 0; i < fN2; ++i ){  
    centerW += fW2.segment(2*i,2);
    centerU += uu.segment(2*i,2);
  }
  centerW /= (double)fN2;
  centerU /= (double)fN2;

  for( int i = 0; i < fNN2; ++i )  
    uu.segment(2*i,2) += (-centerU + centerW);
}


bool TOperator_I::Check_Intersection( VectorXd& u, VectorXd& ua )
{
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;
  for( int i = 0; i < fN; ++ i ){
    for( int j = 0; j < i; ++ j ){
      x1 = u(ix(i));
      y1 = u(iy(i));
      x2 = u(ix(i+1));
      y2 = u(iy(i+1));
      x3 = u(ix(j));
      y3 = u(iy(j));
      x4 = u(ix(j+1));
      y4 = u(iy(j+1));

      tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
      value1 = tc*td;

      tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
      value2 = tc*td;

      if( value1 < 0 && value2 < 0 )
	return true;
    }
  }

  for( int i = 0; i < fNa+1; ++ i ){
    for( int j = 0; j < i; ++ j ){
      x1 = ua(2*i);
      y1 = ua(2*i+1);
      x2 = ua(2*(i+1));
      y2 = ua(2*(i+1)+1);
      x3 = ua(2*j);
      y3 = ua(2*j+1);
      x4 = ua(2*(j+1));
      y4 = ua(2*(j+1)+1);

      tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
      value1 = tc*td;

      tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
      value2 = tc*td;

      if( value1 < 0 && value2 < 0 )
	return true;
    }
  }

  for( int i = 0; i < fN; ++ i ){
    for( int j = 0; j < fNa+1; ++ j ){
      x1 = u(ix(i));
      y1 = u(iy(i));
      x2 = u(ix(i+1));
      y2 = u(iy(i+1));
      x3 = ua(2*j);
      y3 = ua(2*j+1);
      x4 = ua(2*(j+1));
      y4 = ua(2*(j+1)+1);

      tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
      value1 = tc*td;

      tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
      value2 = tc*td;

      if( value1 < 0 && value2 < 0 )
	return true;
    }
  }


  return false;
}

