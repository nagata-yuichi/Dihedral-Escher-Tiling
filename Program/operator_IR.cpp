/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_IR__
#include "operator_IR.h"
#endif


TOperator_IR::TOperator_IR( int N1, int N1_in, int N2, int N2_in ) : TOperator_I( N1, N1_in, N2, N2_in )
{
  fR1_cos = VectorXd::Zero(fNN1); 
  fR1_sin = VectorXd::Zero(fNN1); 
  fRot1 = new MatrixXd [fNN1];
  for( int i = 0; i < fNN1; ++i )
    fRot1[i] = MatrixXd(2,2);
  fRotD1 = new MatrixXd [fNN1];
  for( int i = 0; i < fNN1; ++i )
    fRotD1[i] = MatrixXd(2,2);

  fGG1 = MatrixXd::Zero(2*fNN1, 2*fNN1);
  fAA1 = SparseMatrix<double>(2*fNN1,2*fNN1);
  fAA1.reserve(8*fNN1);
  fTTWW1 = VectorXd::Zero(2*fNN1); 
  fTTWW1_d = VectorXd::Zero(2*fNN1);
  fWW1 = VectorXd::Zero(2*fNN1);
  fWW1_d = VectorXd::Zero(2*fNN1);

  fR2_cos = VectorXd::Zero(fNN2); 
  fR2_sin = VectorXd::Zero(fNN2); 
  fRot2 = new MatrixXd [fNN2];
  for( int i = 0; i < fNN2; ++i )
    fRot2[i] = MatrixXd(2,2);
  fRotD2 = new MatrixXd [fNN2];
  for( int i = 0; i < fNN2; ++i )
    fRotD2[i] = MatrixXd(2,2);

  fGG2 = MatrixXd::Zero(2*fNN2, 2*fNN2);
  fAA2 = SparseMatrix<double>(2*fNN2,2*fNN2);
  fAA2.reserve(8*fNN2);
  fTTWW2 = VectorXd::Zero(2*fNN2); 
  fTTWW2_d = VectorXd::Zero(2*fNN2);
  fWW2 = VectorXd::Zero(2*fNN2);
  fWW2_d = VectorXd::Zero(2*fNN2);
} 

TOperator_IR::~TOperator_IR() 
{

} 

void TOperator_IR::SetParameter() 
{
  TOperator_I::SetParameter();
}

void TOperator_IR::SetInit( VectorXd w1,  VectorXd w1_in, MatrixXi kk1i, VectorXd w2,  VectorXd w2_in, MatrixXi kk2i )
{
  TOperator_I::SetInit( w1, w1_in, kk1i, w2, w2_in, kk2i );

  fWW1.head(2*fN1) = fW1;
  fWW1.tail(2*fN1_in) = fW1_in;
  fWW1_d.head(2*fN1) = fW1_d;
  for( int i = 0; i < fN1_in; ++i ){
    fWW1_d(iix1(fN1+i)) = fW1_in(iiy1(i));
    fWW1_d(iiy1(fN1+i)) = -fW1_in(iix1(i));
  }

  fWW2.head(2*fN2) = fW2;
  fWW2.tail(2*fN2_in) = fW2_in;
  fWW2_d.head(2*fN2) = fW2_d;
  for( int i = 0; i < fN2_in; ++i ){
    fWW2_d(iix2(fN2+i)) = fW2_in(iiy2(i));
    fWW2_d(iiy2(fN2+i)) = -fW2_in(iix2(i));
  }

  this->SetK1(); // operator_Iとの重複がある（なおしたい）
  this->SetK2(); 
}

void TOperator_IR::Cal_u_SGD()
{
  double eval; 
  VectorXd u1(2*fN1);
  VectorXd u2(2*fN2);
  VectorXd u(2*fN);     // points of the template
  VectorXd ua(2*fNa+4); // points of the ua (including both ends)

  VectorXd u1_in(2*fN1_in); 
  VectorXd uu1(2*fNN1);
  VectorXd uu1_d(2*fNN1); 
  VectorXd u2_in(2*fN2_in);
  VectorXd uu2(2*fNN2);
  VectorXd uu2_d(2*fNN2);


  ///////////////////////////////////////////////////////

  for( int st0 = 0; st0 < fN; ++st0 ){
    this->Construct_B1_B2( st0 ); 
    
    for( int st1 = 0; st1 < fN1; ++st1 ){ 
      this->Set_B1j( st1 );
      for( int st2 = 0; st2 < fN2 ; ++st2 ){  
	this->Set_B2j( st2 );

	eval = this->SGD_main( u1, u2, u, ua, uu1, uu2 );

	fEval = 999999999.9;
	if( eval < fEval_best ){ 
	  if( Check_Intersection( u, ua ) == false ){ 
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



void TOperator_IR::Cal_u_SGD( int st0, int st1, int st2 ) 
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

  eval = this->SGD_main( u1, u2, u, ua, uu1, uu2 );

  fEval = 999999999.9;
  if( eval < fEval_best ){ 
    if( Check_Intersection( u, ua ) == false ){
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


double TOperator_IR::SGD_main( VectorXd& u1, VectorXd& u2, VectorXd& u, VectorXd& ua, VectorXd& uu1, VectorXd& uu2 )
{
  VectorXd gzai(fMMd+fMMs), gzai_prev(fMMd+fMMs);
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
  double eval1_const, eval2_const;
  double alpha, alpha_best, alpha_init, alpha_min, alpha_max;
  
  VectorXd hw1(2*fN1);
  VectorXd hw1_d(2*fN1);  
  VectorXd hw2(2*fN2);
  VectorXd hw2_d(2*fN2);  

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
  double d1, d2, alpha1, beta1, alpha2, beta2;
  double ratio = (double)fN2/(double)fN1;

  static VectorXd gzai_previous = VectorXd::Random(2*(fN1+fN2));
  static int MMd_previous = 2*(fN1+fN2);
  static int MMs_previous = 2*(fN1+fN2);

  if( fMMd == MMd_previous && fMMs == MMs_previous ) 
    gzai_curr = gzai_previous.head(fMMd+fMMs); // option (useful in LS)
  else
    gzai_curr = VectorXd::Random(fMMd+fMMs);
  gzai_curr = gzai_curr.normalized();


  fTTWW1 = fGG1*fWW1;
  fTTWW1_d = fGG1*fWW1_d;
  fTTWW2 = fGG2*fWW2;
  fTTWW2_d = fGG2*fWW2_d;
 
  int iter_IR;
  iter_IR = 0;

 START:
  eval1_const = fEval1_wGGw - fTTWW1.bottomRows(2*fN1_in).transpose()*fG2_inv1*fTTWW1.bottomRows(2*fN1_in);
  eval2_const = fEval2_wGGw - fTTWW2.bottomRows(2*fN2_in).transpose()*fG2_inv2*fTTWW2.bottomRows(2*fN2_in);

  // p1, p1_d 
  hw1 = fTTWW1.topRows(2*fN1) + fH1.transpose()*fTTWW1.bottomRows(2*fN1_in); 
  p1.head(fMMd) = fB1dj.transpose()*hw1;
  p1.tail(fMMs) = fB1sj.transpose()*hw1;
  hw1_d = fTTWW1_d.topRows(2*fN1) + fH1.transpose()*fTTWW1_d.bottomRows(2*fN1_in); 
  p1_d.head(fMMd) = fB1dj.transpose()*hw1_d;
  p1_d.tail(fMMs) = fB1sj.transpose()*hw1_d;

  // p2, p2_d 
  hw2 = fTTWW2.topRows(2*fN2) + fH2.transpose()*fTTWW2.bottomRows(2*fN2_in); 
  p2.head(fMMd) = fB2dj.transpose()*hw2;
  p2.tail(fMMs) = fB2sj.transpose()*hw2;
  hw2_d = fTTWW2_d.topRows(2*fN2) + fH2.transpose()*fTTWW2_d.bottomRows(2*fN2_in); 
  p2_d.head(fMMd) = fB2dj.transpose()*hw2_d;
  p2_d.tail(fMMs) = fB2sj.transpose()*hw2_d;

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

  eval1 = eval1_const - gag1_curr/gcg1_curr; 
  eval2 = eval2_const - gag2_curr/gcg2_curr;
  eval1 /= fEval1_wGw;  // 保留 （とりあえず）
  eval2 /= fEval2_wGw;  // 保留 （とりあえず）
  eval_curr = eval1 + ratio*eval2; 
  

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

	eval1 = eval1_const - gag1/gcg1;	
	eval2 = eval2_const - gag2/gcg2;
	eval1 /= fEval1_wGw;  // 保留 （とりあえず）
	eval2 /= fEval2_wGw;  // 保留 （とりあえず）
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

    eval1 = eval1_const - gag1/gcg1;	
    eval2 = eval2_const - gag2/gcg2;
    eval1 /= fEval1_wGw;  // 保留 （とりあえず）
    eval2 /= fEval2_wGw;  // 保留 （とりあえず）
    eval = eval1 + ratio*eval2;
    assert( eval1 > 0.0 );
    assert( eval2 > 0.0 );
    
    // printf( "iter_IR = %d i=%d (%lf): %lf -> %lf \n", iter_IR, iter, alpha, eval_curr, eval );
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

  VectorXd uur1(2*fNN1), uur2(2*fNN2);
  VectorXd ur1(2*fN1), ur2(2*fN2);
  VectorXd ur1_in(2*fN1_in), ur2_in(2*fN2_in);
  double ugw1, ugw1_d, ugw2, ugw2_d;
  double ugu1, ugu2;
  double cos1_th, sin1_th, cos2_th, sin2_th;

  ugw1 = u1.dot(hw1);
  ugw1_d = u1.dot(hw1_d);
  ugw2 = u2.dot(hw2);
  ugw2_d = u2.dot(hw2_d);

  tmp1 = fG1*u1;
  ugu1 = u1.dot(tmp1);
  tmp2 = fG2*u2;
  ugu2 = u2.dot(tmp2);

  alpha1 = ugw1/ugu1;
  beta1 = ugw1_d/ugu1;
  alpha2 = ugw2/ugu2;
  beta2 = ugw2_d/ugu2;
  d1 = sqrt(alpha1*alpha1 + beta1*beta1);
  d2 = sqrt(alpha2*alpha2 + beta2*beta2);
  cos1_th = alpha1/d1;
  sin1_th = beta1/d1;
  cos2_th = alpha2/d2;
  sin2_th = beta2/d2;

  for( int k = 0; k < fN1; ++k ){
    ur1(ix1(k)) = cos1_th*u1(ix1(k)) - sin1_th*u1(iy1(k)); 
    ur1(iy1(k)) = sin1_th*u1(ix1(k)) + cos1_th*u1(iy1(k)); 
  }
  for( int k = 0; k < fN2; ++k ){
    ur2(ix2(k)) = cos2_th*u2(ix2(k)) - sin2_th*u2(iy2(k)); 
    ur2(iy2(k)) = sin2_th*u2(ix2(k)) + cos2_th*u2(iy2(k)); 
  }
  ur1 *= d1;
  ur2 *= d2;

  ur1_in = fG2_inv1*fTTWW1.bottomRows(2*fN1_in) + fH1*ur1;  
  uur1.head(2*fN1) = ur1;        
  uur1.tail(2*fN1_in) = ur1_in;  
  ur2_in = fG2_inv2*fTTWW2.bottomRows(2*fN2_in) + fH2*ur2;  
  uur2.head(2*fN2) = ur2;        
  uur2.tail(2*fN2_in) = ur2_in;
  
  this->Set_RotationMatrices1( uur1 );
  this->Set_RotationMatrices2( uur2 );

  if( iter_IR == 0 ){
    gzai_previous.head(fMMd+fMMs) = gzai;
    MMd_previous = fMMd;
    MMs_previous = fMMs;
  }

  ++iter_IR;
  if( iter_IR < 3 ) 
    goto START;

  eval = this->Cal_eval( u1, u2 );
  // this->Check_eval( u1, u2, eval ); // for check

  u = fBd*gzai.head(fMd) + fBs*gzai.segment(fMd,fMs);

  VectorXd tmp(2*fN1);
  tmp = fB1d*gzai.head(fMMd) + fB1s*gzai.segment(fMMd,fMMs);
  ua.head(2*fNa+2) = tmp.tail(2*fNa+2);
  ua(2*(fNa+1)) = tmp(0);
  ua(2*(fNa+1)+1) = tmp(1);

  for( int k = 0; k < fNN1; ++k ){
    uu1(iix1(k)) = cos1_th*uur1(iix1(k)) + sin1_th*uur1(iiy1(k)); 
    uu1(iiy1(k)) = -sin1_th*uur1(iix1(k)) + cos1_th*uur1(iiy1(k)); 
  }
  uu1 *= 1.0/d1;
  for( int k = 0; k < fNN2; ++k ){
    uu2(iix2(k)) = cos2_th*uur2(iix2(k)) + sin2_th*uur2(iiy2(k)); 
    uu2(iiy2(k)) = -sin2_th*uur2(iix2(k)) + cos2_th*uur2(iiy2(k)); 
  }
  uu2 *= 1.0/d2;
  

  if( fReverse2 == 1 )
    for( int k = 0; k < fNN2; ++k )
      uu2(iix2(k)) *= -1.0;
  
  return eval;
}


double TOperator_IR::Cal_eval( VectorXd u1, VectorXd u2 )
{
  VectorXd hw1(2*fN1);
  VectorXd hw1_d(2*fN1);  
  VectorXd hw2(2*fN2);
  VectorXd hw2_d(2*fN2);  
  double p1, p1_d, p2, p2_d;
  double eval1_const, eval2_const;
  double gcg1, gag1, gcg2, gag2;
  double eval1, eval2, eval;
  double ratio = (double)fN2/(double)fN1;
  
  eval1_const = fEval1_wGGw - fTTWW1.bottomRows(2*fN1_in).transpose()*fG2_inv1*fTTWW1.bottomRows(2*fN1_in);
  eval2_const = fEval2_wGGw - fTTWW2.bottomRows(2*fN2_in).transpose()*fG2_inv2*fTTWW2.bottomRows(2*fN2_in);

  hw1 = fTTWW1.topRows(2*fN1) + fH1.transpose()*fTTWW1.bottomRows(2*fN1_in);
  hw1_d = fTTWW1_d.topRows(2*fN1) + fH1.transpose()*fTTWW1_d.bottomRows(2*fN1_in); 
  p1 = u1.transpose()*hw1;
  p1_d = u1.transpose()*hw1_d;

  hw2 = fTTWW2.topRows(2*fN2) + fH2.transpose()*fTTWW2.bottomRows(2*fN2_in);
  hw2_d = fTTWW2_d.topRows(2*fN2) + fH2.transpose()*fTTWW2_d.bottomRows(2*fN2_in); 
  p2 = u2.transpose()*hw2;
  p2_d = u2.transpose()*hw2_d;

  gcg1 = u1.transpose()*fG1*u1; 
  gag1 = p1*p1 + p1_d*p1_d;
  gcg2 = u2.transpose()*fG2*u2; 
  gag2 = p2*p2 + p2_d*p2_d;

  eval1 = eval1_const - gag1/gcg1; 
  eval2 = eval2_const - gag2/gcg2;
  eval1 /= fEval1_wGw;  // 保留
  eval2 /= fEval2_wGw;  // 保留
  eval = eval1 + ratio*eval2; 

  return eval;
}


void TOperator_IR::Check_eval( VectorXd u1,  VectorXd u2, double eval )
{
  VectorXd hw1(2*fN1);
  VectorXd hw1_d(2*fN1);  
  VectorXd hw2(2*fN2);
  VectorXd hw2_d(2*fN2);
  VectorXd uu1(2*fNN1), uu2(2*fNN2);
  VectorXd ur1_in(2*fN1_in), ur2_in(2*fN2_in);
  VectorXd ur1(2*fN1), ur2(2*fN2);
  double ugw1, ugw1_d, ugw2, ugw2_d;
  double ugu1, ugu2;
  double cos1_th, sin1_th, cos2_th, sin2_th;
  VectorXd v(2); 
  double sum1, sum2, sum;
  VectorXd ttmp1(2*fN1), ttmp2(2*fN2);
  double ratio = (double)fN2/(double)fN1;

  hw1 = fTTWW1.topRows(2*fN1) + fH1.transpose()*fTTWW1.bottomRows(2*fN1_in); 
  hw1_d = fTTWW1_d.topRows(2*fN1) + fH1.transpose()*fTTWW1_d.bottomRows(2*fN1_in); 
  hw2 = fTTWW2.topRows(2*fN2) + fH2.transpose()*fTTWW2.bottomRows(2*fN2_in); 
  hw2_d = fTTWW2_d.topRows(2*fN2) + fH2.transpose()*fTTWW2_d.bottomRows(2*fN2_in); 

  ugw1 = u1.dot(hw1);
  ugw1_d = u1.dot(hw1_d);
  ugw2 = u2.dot(hw2);
  ugw2_d = u2.dot(hw2_d);

  ugu1 = u1.transpose()*fG1*u1;
  ugu2 = u2.transpose()*fG2*u2;

  cos1_th = ugw1/ugu1;
  sin1_th = ugw1_d/ugu1;
  cos2_th = ugw2/ugu2;
  sin2_th = ugw2_d/ugu2; 

  for( int k = 0; k < fN1; ++k ){
    ur1(ix1(k)) = cos1_th*u1(ix1(k)) - sin1_th*u1(iy1(k)); 
    ur1(iy1(k)) = sin1_th*u1(ix1(k)) + cos1_th*u1(iy1(k)); 
  }
  for( int k = 0; k < fN2; ++k ){
    ur2(ix2(k)) = cos2_th*u2(ix2(k)) - sin2_th*u2(iy2(k)); 
    ur2(iy2(k)) = sin2_th*u2(ix2(k)) + cos2_th*u2(iy2(k)); 
  }
  
  ur1_in = fG2_inv1*fTTWW1.bottomRows(2*fN1_in) + fH1*ur1;  
  uu1.head(2*fN1) = ur1;        
  uu1.tail(2*fN1_in) = ur1_in;  
  ur2_in = fG2_inv2*fTTWW2.bottomRows(2*fN2_in) + fH2*ur2;  
  uu2.head(2*fN2) = ur2;        
  uu2.tail(2*fN2_in) = ur2_in;

  // fR1_cos(i)などが計算済でなくてはならない (
  sum1 = 0.0;
  for( int i = 0; i < fNN1; ++i ){
    int size = fNei1_Size(i);

    for( int s = 0; s < size; ++s ){
      int j = fNei1(i,s); 
      double aij = fA1(i,j);
      assert( i != j ); 
      v(0) = uu1(iix1(j)) - uu1(iix1(i)) - fR1_cos(i)*(fWW1(iix1(j)) - fWW1(iix1(i))) + fR1_sin(i)*(fWW1(iiy1(j)) - fWW1(iiy1(i)));
      v(1) = uu1(iiy1(j)) - uu1(iiy1(i)) - fR1_sin(i)*(fWW1(iix1(j)) - fWW1(iix1(i))) - fR1_cos(i)*(fWW1(iiy1(j)) - fWW1(iiy1(i)));
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
      v(0) = uu2(iix2(j)) - uu2(iix2(i)) - fR2_cos(i)*(fWW2(iix2(j)) - fWW2(iix2(i))) + fR2_sin(i)*(fWW2(iiy2(j)) - fWW2(iiy2(i)));
      v(1) = uu2(iiy2(j)) - uu2(iiy2(i)) - fR2_sin(i)*(fWW2(iix2(j)) - fWW2(iix2(i))) - fR2_cos(i)*(fWW2(iiy2(j)) - fWW2(iiy2(i)));
      sum2 += aij*(v(0)*v(0)+v(1)*v(1));
    }
  }
	
  // v(0) = uu1(iix1(0)) - fWW1(iix1(0));
  // v(1) = uu1(iiy1(0)) - fWW1(iiy1(0));
  // sum1 += v(0)*v(0)+v(1)*v(1);
  // v(0) = uu2(iix2(0)) - fWW2(iix2(0));
  // v(1) = uu2(iiy2(0)) - fWW2(iiy2(0));
  // sum2 += v(0)*v(0)+v(1)*v(1);

  sum = sum1/fEval1_wGw + ratio*sum2/fEval2_wGw; // 保留
  
  printf( "%lf %lf -> %lf\n", sum , eval, eval-sum ); 
  assert( fabs(eval-sum) < 0.000001 ); 
}

void TOperator_IR::Scale_Rotate_UU1_Procrustes( VectorXd u1, VectorXd& uu1 )
{
  VectorXd hw1(2*fN1);
  VectorXd hw1_d(2*fN1);  
  VectorXd ur1_in(2*fN1_in);
  VectorXd ur1(2*fN1);
  double ugw1, ugw1_d, ugu1; 
  double cos1_th, sin1_th;

  // fTTWW1は計算済でなくてはならない
  
  hw1 = fTTWW1.topRows(2*fN1) + fH1.transpose()*fTTWW1.bottomRows(2*fN1_in); 
  hw1_d = fTTWW1_d.topRows(2*fN1) + fH1.transpose()*fTTWW1_d.bottomRows(2*fN1_in); 

  ugw1 = u1.dot(hw1);
  ugw1_d = u1.dot(hw1_d);
  ugu1 = u1.transpose()*fG1*u1;

  cos1_th = ugw1/ugu1;
  sin1_th = ugw1_d/ugu1;

  for( int k = 0; k < fN1; ++k ){
    ur1(ix1(k)) = cos1_th*u1(ix1(k)) - sin1_th*u1(iy1(k)); 
    ur1(iy1(k)) = sin1_th*u1(ix1(k)) + cos1_th*u1(iy1(k)); 
  }
  
  ur1_in = fG2_inv1*fTTWW1.bottomRows(2*fN1_in) + fH1*ur1;  
  uu1.head(2*fN1) = ur1;        
  uu1.tail(2*fN1_in) = ur1_in;  
}

void TOperator_IR::Scale_Rotate_UU2_Procrustes( VectorXd u2, VectorXd& uu2 )
{
  VectorXd hw2(2*fN2);
  VectorXd hw2_d(2*fN2);  
  VectorXd ur2_in(2*fN2_in);
  VectorXd ur2(2*fN2);
  double ugw2, ugw2_d, ugu2;
  double cos2_th, sin2_th;

  // fTTWW2は計算済でなくてはならない
  
  hw2 = fTTWW2.topRows(2*fN2) + fH2.transpose()*fTTWW2.bottomRows(2*fN2_in); 
  hw2_d = fTTWW2_d.topRows(2*fN2) + fH2.transpose()*fTTWW2_d.bottomRows(2*fN2_in); 

  ugw2 = u2.dot(hw2);
  ugw2_d = u2.dot(hw2_d);
  ugu2 = u2.transpose()*fG2*u2;

  cos2_th = ugw2/ugu2;
  sin2_th = ugw2_d/ugu2; 

  for( int k = 0; k < fN2; ++k ){
    ur2(ix2(k)) = cos2_th*u2(ix2(k)) - sin2_th*u2(iy2(k)); 
    ur2(iy2(k)) = sin2_th*u2(ix2(k)) + cos2_th*u2(iy2(k)); 
  }
  
  ur2_in = fG2_inv2*fTTWW2.bottomRows(2*fN2_in) + fH2*ur2;  
  uu2.head(2*fN2) = ur2;        
  uu2.tail(2*fN2_in) = ur2_in;
}


void TOperator_IR::SetK1() 
{
  TOperator_I::SetK1();

  for( int i = 0; i < fNN1; ++i ){
    for( int j = 0; j < fNN1; ++j ){
      fGG1(iix1(i),iix1(j)) = fKK1(i,j);
      fGG1(iiy1(i),iiy1(j)) = fKK1(i,j);
    }
  }
  
  fEval1_wGGw = fWW1.transpose()*fGG1*fWW1;

  MatrixXd AA =  MatrixXd::Zero(2*fNN1,2*fNN1); 
  for( int i = 0; i < fNN1; ++i ){
    for( int j = 0; j < fNN1; ++j ){
      AA(iix1(i),iix1(j)) = fA1(i,j);
      AA(iiy1(i),iiy1(j)) = fA1(i,j);
    }
  }

  double val;
  fTripletList.clear();
  for( int j = 0; j < 2*fNN1; ++j ){
    for( int i = 0; i < 2*fNN1; ++i ){
      val = AA(i,j); 
      if( fabs(val) > 0.000000001 ) 
	fTripletList.push_back(T(i,j,val)); 
    }
  }
  fAA1.setFromTriplets(fTripletList.begin(), fTripletList.end());
}

void TOperator_IR::SetK2() 
{
  TOperator_I::SetK2();
 
  for( int i = 0; i < fNN2; ++i ){
    for( int j = 0; j < fNN2; ++j ){
      fGG2(iix2(i),iix2(j)) = fKK2(i,j);
      fGG2(iiy2(i),iiy2(j)) = fKK2(i,j);
    }
  }
  
  fEval2_wGGw = fWW2.transpose()*fGG2*fWW2;

  MatrixXd AA =  MatrixXd::Zero(2*fNN2,2*fNN2); 
  for( int i = 0; i < fNN2; ++i ){
    for( int j = 0; j < fNN2; ++j ){
      AA(iix2(i),iix2(j)) = fA2(i,j);
      AA(iiy2(i),iiy2(j)) = fA2(i,j);
    }
  }

  double val;
  fTripletList.clear();
  for( int j = 0; j < 2*fNN2; ++j ){
    for( int i = 0; i < 2*fNN2; ++i ){
      val = AA(i,j); 
      if( fabs(val) > 0.000000001 ) 
	fTripletList.push_back(T(i,j,val)); 
    }
  }
  fAA2.setFromTriplets(fTripletList.begin(), fTripletList.end());
}


void TOperator_IR::Set_RotationMatrices1( VectorXd uu )
{
  int j,size;
  VectorXd ud(2*fNN1); 
  VectorXd wd(2*fNN1);
  VectorXd wd_c(2*fNN1);
  double tr, rt, sqr;
  double di, aij;
  VectorXd RRWW(2*fNN1);

  for( int i = 0; i < fNN1; ++i ){
    size = fNei1_Size(i);
    for( int s = 0; s < size; ++s ){
      j = fNei1(i,s); 
      aij = sqrt(fA1(i,j)); // 注意
      ud(2*s) = aij*(uu(iix1(j)) - uu(iix1(i)));  
      ud(2*s+1) = aij*(uu(iiy1(j)) - uu(iiy1(i))); 
      wd(2*s) = aij*(fWW1(iix1(j)) - fWW1(iix1(i)));
      wd(2*s+1) = aij*(fWW1(iiy1(j)) - fWW1(iiy1(i)));
      wd_c(2*s) = wd(2*s+1); 
      wd_c(2*s+1) = -wd(2*s);
    }

    tr =ud.head(2*size).dot(wd.head(2*size));
    rt =ud.head(2*size).dot(wd_c.head(2*size));
    sqr = sqrt( tr*tr + rt*rt );
    fR1_cos(i) = tr/sqr;  // thetaはUを回転させてWに一致させた場合
    fR1_sin(i) = -rt/sqr; // Wを回転させるので、-thetaの回転
    // fR1_cos(i) = 1.0; // test
    // fR1_sin(i) = 0.0;
    fRot1[i](0,0) = fR1_cos(i);  
    fRot1[i](0,1) = -fR1_sin(i);
    fRot1[i](1,0) = fR1_sin(i);
    fRot1[i](1,1) = fR1_cos(i);
  }

  for( int i = 0; i < fNN1; ++i ){
    size = fNei1_Size(i);
    di = fDout1(i,i);
    fRotD1[i] = di*fRot1[i]; 

    for( int s = 0; s < size; ++s ){
      j = fNei1(i,s); 
      fRotD1[i] += fA1(j,i)*fRot1[j];  // ここ注意： fA(i,j)でない
    }
  }
  // fRotD1[0](0,0) += 1.0; // 正則化の項
  // fRotD1[0](1,1) += 1.0;


  fTTWW1 = VectorXd::Zero(2*fNN1);
  RRWW = VectorXd::Zero(2*fNN1);
  for( int i = 0; i < fNN1; ++i ){
    fTTWW1.segment(2*i,2) += fRotD1[i]*fWW1.segment(2*i,2); 
    RRWW.segment(2*i,2) += fRot1[i]*fWW1.segment(2*i,2); 
  }
  fTTWW1 += -fAA1.transpose()*RRWW;

  VectorXd tmp(2*fNN1);
  tmp = fAA1*fWW1;
  for( int i = 0; i < fNN1; ++i )
    fTTWW1.segment(2*i,2) -= fRot1[i]*tmp.segment(2*i,2);

  for( int i = 0; i < fNN1; ++i ){
    fTTWW1_d(iix1(i)) = fTTWW1(iiy1(i));
    fTTWW1_d(iiy1(i)) = -fTTWW1(iix1(i));
  }

}

void TOperator_IR::Set_RotationMatrices2( VectorXd uu )
{
  int j,size;
  VectorXd ud(2*fNN2); 
  VectorXd wd(2*fNN2);
  VectorXd wd_c(2*fNN2);
  double tr, rt, sqr;
  double di, aij;
  VectorXd RRWW(2*fNN2);

  for( int i = 0; i < fNN2; ++i ){
    size = fNei2_Size(i);
    for( int s = 0; s < size; ++s ){
      j = fNei2(i,s); 
      aij = sqrt(fA2(i,j)); // 注意
      ud(2*s) = aij*(uu(iix2(j)) - uu(iix2(i)));  
      ud(2*s+1) = aij*(uu(iiy2(j)) - uu(iiy2(i))); 
      wd(2*s) = aij*(fWW2(iix2(j)) - fWW2(iix2(i)));
      wd(2*s+1) = aij*(fWW2(iiy2(j)) - fWW2(iiy2(i)));
      wd_c(2*s) = wd(2*s+1); 
      wd_c(2*s+1) = -wd(2*s);
    }

    tr =ud.head(2*size).dot(wd.head(2*size));
    rt =ud.head(2*size).dot(wd_c.head(2*size));
    sqr = sqrt( tr*tr + rt*rt );
    fR2_cos(i) = tr/sqr;  // thetaはUを回転させてWに一致させた場合
    fR2_sin(i) = -rt/sqr; // Wを回転させるので、-thetaの回転
    // fR2_cos(i) = 1.0; // test
    //  fR2_sin(i) = 0.0;
    fRot2[i](0,0) = fR2_cos(i);  
    fRot2[i](0,1) = -fR2_sin(i);
    fRot2[i](1,0) = fR2_sin(i);
    fRot2[i](1,1) = fR2_cos(i);
  }

  for( int i = 0; i < fNN2; ++i ){
    size = fNei2_Size(i);
    di = fDout2(i,i);
    fRotD2[i] = di*fRot2[i]; 

    for( int s = 0; s < size; ++s ){
      j = fNei2(i,s); 
      fRotD2[i] += fA2(j,i)*fRot2[j];  // ここ注意： fA(i,j)でない
    }
  }
  // fRotD2(0,0) += 1.0; // 正則化の項
  // fRotD2(1,1) += 1.0;


  fTTWW2 = VectorXd::Zero(2*fNN2);
  RRWW = VectorXd::Zero(2*fNN2);
  for( int i = 0; i < fNN2; ++i ){
    fTTWW2.segment(2*i,2) += fRotD2[i]*fWW2.segment(2*i,2); 
    RRWW.segment(2*i,2) += fRot2[i]*fWW2.segment(2*i,2); 
  }
  fTTWW2 += -fAA2.transpose()*RRWW; 

  VectorXd tmp(2*fNN2);
  tmp = fAA2*fWW2;
  for( int i = 0; i < fNN2; ++i )
    fTTWW2.segment(2*i,2) -= fRot2[i]*tmp.segment(2*i,2); 

  for( int i = 0; i < fNN2; ++i ){
    fTTWW2_d(iix2(i)) = fTTWW2(iiy2(i));
    fTTWW2_d(iiy2(i)) = -fTTWW2(iix2(i));
  }
}

// I distanceに基いて、W1にmatchするようスケールと回転変換（重心も）
void TOperator_IR::Scale_Rotate_UU1( VectorXd& uu1 )
{
  VectorXd uu1_t(2*fNN1);
  double ugu, ugw, ugw_d;
  ugu = uu1.transpose()*fGG1*uu1;
  ugw = uu1.transpose()*fGG1*fWW1;
  ugw_d = uu1.transpose()*fGG1*fWW1_d;
  double cos_th = ugw / ugu;
  double sin_th = ugw_d / ugu;

  for( int k = 0; k < fNN1; ++k ){
    uu1_t(iix1(k)) = cos_th*uu1(iix1(k)) - sin_th*uu1(iiy1(k)); 
    uu1_t(iiy1(k)) = sin_th*uu1(iix1(k)) + cos_th*uu1(iiy1(k));
  }
  uu1 = uu1_t;

  this->Trans_UU1_center( uu1 );
}

void TOperator_IR::Scale_Rotate_UU2( VectorXd& uu2 )
{
  VectorXd uu2_t(2*fNN2);
  double ugu, ugw, ugw_d;
  ugu = uu2.transpose()*fGG2*uu2;
  ugw = uu2.transpose()*fGG2*fWW2;
  ugw_d = uu2.transpose()*fGG2*fWW2_d;
  double cos_th = ugw / ugu;
  double sin_th = ugw_d / ugu;

  for( int k = 0; k < fNN2; ++k ){
    uu2_t(iix2(k)) = cos_th*uu2(iix2(k)) - sin_th*uu2(iiy2(k)); 
    uu2_t(iiy2(k)) = sin_th*uu2(iix2(k)) + cos_th*uu2(iiy2(k));
  }
  uu2 = uu2_t;

  this->Trans_UU2_center( uu2 );
}
