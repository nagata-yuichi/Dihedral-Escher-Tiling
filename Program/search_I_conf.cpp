/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_CONF__
#include "search_I_conf.h"
#endif


TSearch_I_conf::TSearch_I_conf( int N1, int N1_in, int N2, int N2_in ) 
{
  fN1 = N1;
  fN1_in = N1_in;
  fN2 = N2;
  fN2_in = N2_in;
  fNN1 = N1+N1_in;
  fNN2 = N2+N2_in;
}

				
TSearch_I_conf::~TSearch_I_conf()
{

}

void TSearch_I_conf::SetParameter()
{
  fNum_top = 100;
  fAlpha = 0.5;
  fNum_threads = 20;
}

void TSearch_I_conf::Define()
{
  fEval_top = new double [ fNum_top ];
  fIndex_top = new int [ fNum_top ];

  fNv_top = new int [ fNum_top ];
  fIH_top = new int [ fNum_top ];
  fN_top = new int [ fNum_top ];
  fNa_top = new int [ fNum_top ];
  fSt0_top = new int [ fNum_top ];
  fSt1_top = new int [ fNum_top ];
  fSt2_top = new int [ fNum_top ];
  ffn_top = new int* [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    ffn_top[i] = new int [ 6 ];
  fReverse_top = new int [ fNum_top ];
  
  fU_top = new VectorXd [ fNum_top ];
  fUa_top = new VectorXd [ fNum_top ];
  fUU1_top = new VectorXd [ fNum_top ];
  fUU2_top = new VectorXd [ fNum_top ];

  for( int i = 0; i < fNum_top; ++i ){
    fU_top[i] = VectorXd::Zero(2*(fN1+fN2));
    fUa_top[i] = VectorXd::Zero(2*(fN1+fN2));
    fUU1_top[i] = VectorXd::Zero(2*(fN1+fN1_in));
    fUU2_top[i] = VectorXd::Zero(2*(fN2+fN2_in));
  }
}

void TSearch_I_conf::SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2, VectorXd w2_in, MatrixXi kk2i )
{
  tOp = new CLASS_NAME* [fNum_threads];
  for( int i = 0; i < fNum_threads; ++i )
    tOp[i] = new CLASS_NAME( fN1, fN1_in, fN2, fN2_in );
  
  for( int i = 0; i < fNum_threads; ++i ){
    tOp[i]->SetParameter();
    tOp[i]->fAlpha = fAlpha;
    tOp[i]->SetInit( w1, w1_in, kk1i, w2, w2_in, kk2i );
  }

  omp_set_num_threads(fNum_threads);
  printf( "#threads = %d \n", fNum_threads );
}

void TSearch_I_conf::DoIt()
{
  fEval_best = 999999999.9;
  for( int i = 0; i < fNum_top; ++i ){
    fEval_top[i] = 999999999.9;
    fIndex_top[i] = i;
}

  printf( "IH1\n");
  this->IH1();
  printf( "IH2\n");
  this->IH2();
  printf( "IH3\n");
  this->IH3();
  printf( "IH4\n");
  this->IH4();
  printf( "IH5\n");
  this->IH5();
  printf( "IH6\n");
  this->IH6();
}

void TSearch_I_conf::IH4() 
{
  int na, k1, k2, k3, k4, st0, st1, st2, reverse;
  int num_candi;

  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH4();
  }

  num_candi = fCandi_end_IH4-fCandi_start_IH4;
#pragma omp parallel for private(na, k1, k2, k3, k4, st0, st1, st2, reverse)
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH4 + (s*num_candi)/fNum_threads; i < fCandi_start_IH4 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      k3 = fki_candi[i][2]; 
      k4 = fki_candi[i][3];
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];

      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH4( k1, k2, k3, k4 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}

void TSearch_I_conf::IH5() 
{
  int na, k1, k2, k3, st0, st1, st2, reverse;
  int num_candi;

  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH5();
  }

  num_candi = fCandi_end_IH5-fCandi_start_IH5;
#pragma omp parallel for private( na, k1, k2, k3, st0, st1, st2, reverse )
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH5 + (s*num_candi)/fNum_threads; i < fCandi_start_IH5 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      k3 = fki_candi[i][2]; 
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];

      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH5( k1, k2, k3 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}


void TSearch_I_conf::IH6() 
{
  int na, k1, k2, k3, st0, st1, st2, reverse;
  int num_candi;

  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH6();
  }

  num_candi = fCandi_end_IH6-fCandi_start_IH6;
#pragma omp parallel for private( na, k1, k2, k3, st0, st1, st2, reverse )
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH6 + (s*num_candi)/fNum_threads; i < fCandi_start_IH6 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      k3 = fki_candi[i][2]; 
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];

      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH6( k1, k2, k3 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}


void TSearch_I_conf::IH1() 
{
  int na, k1, k2, st0, st1, st2, reverse;
  int num_candi;
  
  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH1();
  }

  num_candi = fCandi_end_IH1-fCandi_start_IH1;
#pragma omp parallel for private( na, k1, k2, st0, st1, st2, reverse )
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH1 + (s*num_candi)/fNum_threads; i < fCandi_start_IH1 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];

      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH1( k1, k2 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}

void TSearch_I_conf::IH2() 
{
  int na, k1, k2, st0, st1, st2, reverse;
  int num_candi;
  
  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH2();
  }

  num_candi = fCandi_end_IH2-fCandi_start_IH2;
#pragma omp parallel for private( na, k1, k2, st0, st1, st2, reverse )
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH2 + (s*num_candi)/fNum_threads; i < fCandi_start_IH2 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];
    
      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH2( k1, k2 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}

void TSearch_I_conf::IH3() 
{
  int na, k1, k2, st0, st1, st2, reverse;
  int num_candi;
  
  for( int s = 0; s < fNum_threads; ++s ){
    tOp[s]->fReverse1 = 0;
    tOp[s]->Make_Bv_Be_IH3();
  }

  num_candi = fCandi_end_IH3-fCandi_start_IH3;
#pragma omp parallel for private( na, k1, k2, st0, st1, st2, reverse )
  for( int s = 0; s < fNum_threads; ++s ){
    for( int i = fCandi_start_IH3 + (s*num_candi)/fNum_threads; i < fCandi_start_IH3 + ((s+1)*num_candi)/fNum_threads; ++i ){
      na = fNa_candi[i];
      k1 = fki_candi[i][0]; 
      k2 = fki_candi[i][1];
      st0 = fSt0_candi[i];
      st1 = fSt1_candi[i];
      st2 = fSt2_candi[i];
      reverse = fReverse_candi[i];

      tOp[s]->Set_values( na, reverse, fEval_best );
      tOp[s]->Make_B_IH3( k1, k2 );
      tOp[s]->Cal_u_SGD( st0, st1, st2 );
      if( tOp[s]->fEval < fEval_best ){
#pragma omp critical
	this->Update_UU_top( tOp[s] );
      }
    }
  }
}


void TSearch_I_conf::Update_UU_top( CLASS_NAME* op )
{
  if( op->fEval >= fEval_best )
    return;

  int N, Na, Nv;
  N = op->fN;
  Na = op->fNa;
  Nv = op->fNv;
  
  static int count_imp = 0;

  for( int s = 0; s < fNum_top; ++s ){
    if( fabs( op->fEval - fEval_top[fIndex_top[s]] ) < 0.0000001 ) 
      return;
  }

  int orderIns;
  for( int s = 0; s < fNum_top; ++s ){
    if( op->fEval < fEval_top[fIndex_top[s]] ){ 
      orderIns = s;
      break;
    }
  }
 
  int index_stock = fIndex_top[fNum_top-1];
  for( int s = fNum_top-1; s >= orderIns+1; --s )
    fIndex_top[s] = fIndex_top[s-1];

  fIndex_top[orderIns] = index_stock;
  fEval_top[fIndex_top[orderIns]] = op->fEval;
  fEval_best = fEval_top[fIndex_top[fNum_top-1]];

  fSt0_top[fIndex_top[orderIns]] = op->fSt0;
  fSt1_top[fIndex_top[orderIns]] = op->fSt1;
  fSt2_top[fIndex_top[orderIns]] = op->fSt2;
  fReverse_top[fIndex_top[orderIns]] = op->fReverse2;

  fUU1_top[fIndex_top[orderIns]] = op->fUU1;
  fUU2_top[fIndex_top[orderIns]] = op->fUU2;

  fIH_top[fIndex_top[orderIns]] = op->fIH;
  fNv_top[fIndex_top[orderIns]] = op->fNv;
  for( int h = 0; h < Nv; ++h )
    ffn_top[fIndex_top[orderIns]][h] = op->fn[h];

  fN_top[fIndex_top[orderIns]] = op->fN;
  fNa_top[fIndex_top[orderIns]] = op->fNa;
  fU_top[fIndex_top[orderIns]].head(2*N) = op->fU.head(2*N);
  fUa_top[fIndex_top[orderIns]].head(2*Na+4) = op->fUa.head(2*Na+4);

  printf( "IH%d(%d): eval=%lf", op->fIH, count_imp, op->fEval );
  printf(": a = %d", op->fNa ); 
  printf(": k=(");
  for( int i = 0; i < Nv; ++i )
    printf( "%d ", op->fn[i] );
  printf(")");
  printf( ": st = %d %d %d : r = %d\n", op->fSt0, op->fSt1, op->fSt2, op->fReverse2 );
  
  ++count_imp;
}


void TSearch_I_conf::Set_candi( int num_candi, double* eval_candi, int* IH_candi, int* Na_candi, int** fn_candi, int* st0_candi, int* st1_candi, int* st2_candi, int* reverse_candi )
{
  int checked;
  int index;
  // num_candi = 1000000; // 読込configuration数の指定
  fNum_candi = num_candi;
  fIH_candi = new int [ fNum_candi ];
  fEval_candi = new double [ fNum_candi ];
  fNa_candi = new int [ fNum_candi ];
  fSt0_candi = new int [ fNum_candi ];
  fSt1_candi = new int [ fNum_candi ];
  fSt2_candi = new int [ fNum_candi ];
  fki_candi = new int* [ fNum_candi ];
  for( int i = 0; i < fNum_candi; ++i )
    fki_candi[i] = new int [ 4 ];
  fReverse_candi = new int [ fNum_candi ];

  checked = 0;

  fCandi_start_IH4 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 4 ){
      fIH_candi[checked] = 4;
      fEval_candi[checked] = eval_candi[s];
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][2];
      fki_candi[checked][3] = fn_candi[s][4];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH4 = checked;

  fCandi_start_IH5 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 5 ){
      fIH_candi[checked] = 5;
      fEval_candi[checked] = eval_candi[s];
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][4];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH5 = checked;

  fCandi_start_IH6 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 6 ){
      fIH_candi[checked] = 6;
      fEval_candi[checked] = eval_candi[s];
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][3];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH6 = checked;

  fCandi_start_IH1 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 1 ){
      fIH_candi[checked] = 1;
      fEval_candi[checked] = eval_candi[s];
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH1 = checked;

  fCandi_start_IH2 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 2 ){
      fIH_candi[checked] = 2;
      fEval_candi[checked] = eval_candi[s];
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH2 = checked;

  fCandi_start_IH3 = checked;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 3 ){
      fIH_candi[checked] = 3;
      fEval_candi[checked] = eval_candi[s];      
      fNa_candi[checked] = Na_candi[s];
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt0_candi[checked] = st0_candi[s];
      fSt1_candi[checked] = st1_candi[s];
      fSt2_candi[checked] = st2_candi[s];
      fReverse_candi[checked] = reverse_candi[s];
      ++checked;
    }
  }
  fCandi_end_IH3 = checked;
  
  assert( checked == num_candi );
}
