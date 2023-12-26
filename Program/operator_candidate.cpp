/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_CANDIDATE__
#include "operator_candidate.h"
#endif

TOperator_Candidate::TOperator_Candidate( int N1, int N2 )
{
  fN1 = N1;
  fN2 = N2;

   int minN;
  if( N1 < N2 ) minN = N1;
  else minN = N2;
  
  fEval_J12_edge = new double** [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    fEval_J12_edge[ i ] = new double* [ fN2 ];
  for( int i = 0; i < fN1; ++i )
    for( int j = 0; j < fN2; ++j )
      fEval_J12_edge[ i ][ j ] = new double [ minN+1 ];

  fEval_J1_edge = new double** [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    fEval_J1_edge[ i ] = new double* [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    for( int j = 0; j < fN1; ++j )
      fEval_J1_edge[ i ][ j ] = new double [ fN1+1 ];

  fEval_J2_edge = new double** [ fN2 ];
  for( int i = 0; i < fN2; ++i )
    fEval_J2_edge[ i ] = new double* [ fN2 ];
  for( int i = 0; i < fN2; ++i )
    for( int j = 0; j < fN2; ++j )
      fEval_J2_edge[ i ][ j ] = new double [ fN2+1 ];
  

  fEval_S1_edge = new double* [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    fEval_S1_edge[ i ] = new double [ fN1+1 ]; 
  fEval_S2_edge = new double* [ fN2 ];
  for( int i = 0; i < fN2; ++i )
    fEval_S2_edge[ i ] = new double [ fN2+1 ];

  fEval_J12_GR_edge = new double** [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    fEval_J12_GR_edge[ i ] = new double* [ fN2 ];
  for( int i = 0; i < fN1; ++i )
    for( int j = 0; j < fN2; ++j )
      fEval_J12_GR_edge[ i ][ j ] = new double [ minN+1 ];

  fEval_J1_GR_edge = new double** [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    fEval_J1_GR_edge[ i ] = new double* [ fN1 ];
  for( int i = 0; i < fN1; ++i )
    for( int j = 0; j < fN1; ++j )
      fEval_J1_GR_edge[ i ][ j ] = new double [ fN1+1 ];

  fEval_J2_GR_edge = new double** [ fN2 ];
  for( int i = 0; i < fN2; ++i )
    fEval_J2_GR_edge[ i ] = new double* [ fN2 ];
  for( int i = 0; i < fN2; ++i )
    for( int j = 0; j < fN2; ++j )
      fEval_J2_GR_edge[ i ][ j ] = new double [ fN2+1 ];

  //  fn = new int [6];

  fW1 = VectorXd::Zero(2*fN1);
  fW2 = VectorXd::Zero(2*fN2);
}

TOperator_Candidate::~TOperator_Candidate() 
{
} 

void TOperator_Candidate::SetParameter() 
{

}

void TOperator_Candidate::SetInit( VectorXd w1, VectorXd w2 )
{
  double length1 = 0.0;
  for( int i = 0; i < fN1; ++i ){
    double x = w1(ix1(i+1))-w1(ix1(i));
    double y = w1(iy1(i+1))-w1(iy1(i));
    length1 += sqrt(x*x+y*y);
  }
  length1 /= (double)fN1;
  
  double length2 = 0.0;
  for( int i = 0; i < fN2; ++i ){
    double x = w2(ix2(i+1))-w2(ix2(i));
    double y = w2(iy2(i+1))-w2(iy2(i));
    length2 += sqrt(x*x+y*y);
  }
  length2 /= (double)fN2;

  fW1 = w1;
  fW2 = w2;
  fW1 /= (double)length1;
  fW2 /= (double)length2;
}


void TOperator_Candidate::Set_eval_JS()
{
  int minN;
  if( fN1 < fN2 ) minN = fN1;
  else minN = fN2;
  double eval;
  VectorXd w1(4*fN1+2);
  VectorXd w1r(2*fN1);
  VectorXd w2(4*fN2+2);
  VectorXd w2r(2*fN2);
  

  if( fReverse1 == 0 ){
    w1.head(2*fN1) = fW1;
    w1.segment(2*fN1,2*fN1) = fW1;
    w1.tail(2) = fW1.head(2);
  }
  else if( fReverse1 == 1 ){
    for( int i = 0; i < fN1; ++i ){
      w1r.row(ix1(i)) = -fW1.row(ix1(-i));
      w1r.row(iy1(i)) = fW1.row(iy1(-i));
    }
    w1.head(2*fN1) = w1r;
    w1.segment(2*fN1,2*fN1) = w1r;
    w1.tail(2) = w1r.head(2);
  }

  
  if( fReverse2 == 0 ){
    w2.head(2*fN2) = fW2;
    w2.segment(2*fN2,2*fN2) = fW2;
    w2.tail(2) = fW2.head(2);
  }
  else if( fReverse2 == 1 ){
    for( int i = 0; i < fN2; ++i ){
      w2r.row(ix2(i)) = -fW2.row(ix2(-i));
      w2r.row(iy2(i)) = fW2.row(iy2(-i));
    }
    w2.head(2*fN2) = w2r;
    w2.segment(2*fN2,2*fN2) = w2r;
    w2.tail(2) = w2r.head(2);
  }

  // J12 edge
  for( int len = 0; len <= minN; ++len ){
    for( int a1 = 0; a1 < fN1; ++a1 ){
      for( int a2 = 0; a2 < fN2; ++a2 ){
	eval = this->Diff_segment_rotation( w1.segment(2*a1,2*len+2), w2.segment(2*a2,2*len+2), len );
	// printf( "%d %d %d : %lf\n", a1, a2, len, eval );
	fEval_J12_edge[a1][a2][len] = eval;
      }
    }
  }

  // J1 edge (重なりを許しているのがムダ）
  for( int len = 0; len <= fN1; ++len ){
    for( int a1 = 0; a1 < fN1; ++a1 ){
      for( int a2 = 0; a2 < fN1; ++a2 ){
	eval = this->Diff_segment_parallel( w1.segment(2*a1,2*len+2), w1.segment(2*a2,2*len+2), len );
	fEval_J1_edge[a1][a2][len] = eval;
      }
    }
  }

  // J2 edge (重なりを許しているのがムダ）
  for( int len = 0; len <= fN2; ++len ){
    for( int a1 = 0; a1 < fN2; ++a1 ){
      for( int a2 = 0; a2 < fN2; ++a2 ){
	eval = this->Diff_segment_parallel( w2.segment(2*a1,2*len+2), w2.segment(2*a2,2*len+2), len );
	fEval_J2_edge[a1][a2][len] = eval;
      }
    }
  }

  // S1 edge 
  for( int len = 0; len <= fN1; ++len ){
    for( int a = 0; a < fN1; ++a ){
      eval = this->Diff_segment_parallel( w1.segment(2*a,2*len+2), -w1.segment(2*a,2*len+2), len );
      fEval_S1_edge[a][len] = eval;
    }
  }

  // S2 edge 
  for( int len = 0; len <= fN2; ++len ){
    for( int a = 0; a < fN2; ++a ){
      eval = this->Diff_segment_parallel( w2.segment(2*a,2*len+2), -w2.segment(2*a,2*len+2), len );
      fEval_S2_edge[a][len] = eval;
    }
  }

  // J12 edge (GR)
  for( int len = 0; len <= minN; ++len ){
    for( int a1 = 0; a1 < fN1; ++a1 ){
      for( int a2 = 0; a2 < fN2; ++a2 ){
	eval = this->Diff_segment_rotation_GR( w1.segment(2*a1,2*len+2), w2.segment(2*a2,2*len+2), len );
	// printf( "%d %d %d : %lf\n", a1, a2, len, eval );
	fEval_J12_GR_edge[a1][a2][len] = eval;
      }
    }
  }

  // J1 edge (GR)
  for( int len = 0; len <= fN1; ++len ){
    for( int a1 = 0; a1 < fN1; ++a1 ){
      for( int a2 = 0; a2 < fN1; ++a2 ){
	eval = this->Diff_segment_rotation_GR( w1.segment(2*a1,2*len+2), w1.segment(2*a2,2*len+2), len );
	fEval_J1_GR_edge[a1][a2][len] = eval;
      }
    }
  }

  // J2 edge (GR)
  for( int len = 0; len <= fN2; ++len ){
    for( int a1 = 0; a1 < fN2; ++a1 ){
      for( int a2 = 0; a2 < fN2; ++a2 ){
	eval = this->Diff_segment_rotation_GR( w2.segment(2*a1,2*len+2), w2.segment(2*a2,2*len+2), len );
	fEval_J2_GR_edge[a1][a2][len] = eval;
      }
    }
  }
  
}

// rotationあり 
double TOperator_Candidate::Diff_segment_rotation( VectorXd segment1, VectorXd segment2, int len )
{
  if( len == 0 )
    return 0.0;

  VectorXd w1(2*len);
  VectorXd w2(2*len);
  VectorXd w1_d(2*len);
  VectorXd w2_d(2*len);
  VectorXd center_w1(2);
  VectorXd center_w2(2);

  for( int i = 0; i < len; ++i ){
    w1(2*i) = segment1[2*(i+1)] - segment1[2*i]; 
    w1(2*i+1) = segment1[2*(i+1)+1] - segment1[2*i+1];
    w2(2*i) = segment2[2*(len-i-1)] - segment2[2*(len-i)]; 
    w2(2*i+1) = segment2[2*(len-i-1)+1] - segment2[2*(len-i)+1]; 
  }

  center_w1 = VectorXd::Zero(2); 
  center_w2 = VectorXd::Zero(2);
  for( int i = 0; i < len; ++i ){  
    center_w1 += w1.segment(2*i,2);
    center_w2 += w2.segment(2*i,2);
  }
  center_w1 /= (double)len;
  center_w2 /= (double)len;
  /*
  for( int i = 0; i < len; ++i ){
    w1.segment(2*i,2) -= center_w1;
    w2.segment(2*i,2) -= center_w2;
  }
  */

  // test
  /*
  printf( "ws1 = \n" );
  cout << w1 << endl;
  printf( "ws2 = \n" );
  cout << w2 << endl;
  */

  for( int i = 0; i < len; ++i ){
    w1_d(2*i) = w1(2*i+1);
    w1_d(2*i+1) = -w1(2*i);
    w2_d(2*i) = w2(2*i+1);
    w2_d(2*i+1) = -w2(2*i);
  }

  double p, p_d;
  p = w2.dot(w1);
  p_d = w2.dot(w1_d);
  double ww1 = w1.dot(w1);
  double ww2 = w2.dot(w2);

  double eval = ww1 + ww2 - 2.0*sqrt(p*p+p_d*p_d);
  assert( eval > -0.000001 );
  
  return eval;
}


double TOperator_Candidate::Diff_segment_parallel( VectorXd segment1, VectorXd segment2, int len )
{
  if( len == 0 )
    return 0.0;

  VectorXd w1(2*len);
  VectorXd w2(2*len);
  VectorXd center_w1(2);
  VectorXd center_w2(2);

  for( int i = 0; i < len; ++i ){
    w1(2*i) = segment1[2*(i+1)] - segment1[2*i]; 
    w1(2*i+1) = segment1[2*(i+1)+1] - segment1[2*i+1];
    w2(2*i) = segment2[2*(len-i-1)] - segment2[2*(len-i)]; 
    w2(2*i+1) = segment2[2*(len-i-1)+1] - segment2[2*(len-i)+1]; 
  }

  center_w1 = VectorXd::Zero(2); 
  center_w2 = VectorXd::Zero(2);
  for( int i = 0; i < len; ++i ){  
    center_w1 += w1.segment(2*i,2);
    center_w2 += w2.segment(2*i,2);
  }
  center_w1 /= (double)len;
  center_w2 /= (double)len;
  /*
  for( int i = 0; i < len; ++i ){  
    w1.segment(2*i,2) -= center_w1;
    w2.segment(2*i,2) -= center_w2;
  }
  */

  // test
  /*
  printf( "ws1 = \n" );
  cout << w1 << endl;
  printf( "ws2 = \n" );
  cout << w2 << endl;
  */

  
  VectorXd w12 = w1 - w2;
  double eval = w12.dot(w12);
  assert( eval > -0.000001 );
  
  return eval;
}

// IH6などのglide reflectionタイプのJ edge
double TOperator_Candidate::Diff_segment_rotation_GR( VectorXd segment1, VectorXd segment2, int len )
{
  if( len == 0 )
    return 0.0;

  VectorXd w1(2*len);
  VectorXd w2(2*len);
  VectorXd w1_d(2*len);
  VectorXd w2_d(2*len);
  VectorXd center_w1(2);
  VectorXd center_w2(2);

  for( int i = 0; i < len; ++i ){
    w1(2*i) = segment1[2*(i+1)] - segment1[2*i]; 
    w1(2*i+1) = segment1[2*(i+1)+1] - segment1[2*i+1];
    w2(2*i) = segment2[2*(i+1)] - segment2[2*i]; 
    w2(2*i+1) = -(segment2[2*(i+1)+1] - segment2[2*i+1]);
  }

  center_w1 = VectorXd::Zero(2); 
  center_w2 = VectorXd::Zero(2);
  for( int i = 0; i < len; ++i ){  
    center_w1 += w1.segment(2*i,2);
    center_w2 += w2.segment(2*i,2);
  }
  center_w1 /= (double)len;
  center_w2 /= (double)len;
  /*
  for( int i = 0; i < len; ++i ){
    w1.segment(2*i,2) -= center_w1;
    w2.segment(2*i,2) -= center_w2;
  }
  */

  // test
  /*
  printf( "ws1 = \n" );
  cout << w1 << endl;
  printf( "ws2 = \n" );
  cout << w2 << endl;
  */

  for( int i = 0; i < len; ++i ){
    w1_d(2*i) = w1(2*i+1);
    w1_d(2*i+1) = -w1(2*i);
    w2_d(2*i) = w2(2*i+1);
    w2_d(2*i+1) = -w2(2*i);
  }

  double p, p_d;
  p = w2.dot(w1);
  p_d = w2.dot(w1_d);
  double ww1 = w1.dot(w1);
  double ww2 = w2.dot(w2);

  double eval = ww1 + ww2 - 2.0*sqrt(p*p+p_d*p_d);
  assert( eval > -0.000001 );
  
  return eval;
}


