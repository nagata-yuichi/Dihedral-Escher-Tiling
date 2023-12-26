/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_GET__
#include "search_I_get_conf.h"
#endif


TSearch_I_get::TSearch_I_get( int N1, int N2 ) : TOperator_Candidate( N1, N2 )    
{
  fN1 = N1;
  fN2 = N2;


  fComp_S1_S1 = new double* [fN1];  
  for( int i = 0; i < fN1; ++i )
    fComp_S1_S1[i] = new double [fN1+1];

  fComp_S2_S2 = new double* [fN2];  
  for( int i = 0; i < fN2; ++i )
    fComp_S2_S2[i] = new double [fN2+1];

  fComp_S1_S12 = new double*** [fN1];
  for( int i1 = 0; i1 < fN1; ++i1 ){
    fComp_S1_S12[i1] = new double** [fN1+1];
    for( int i2 = 0; i2 <= fN1; ++i2 ){
      fComp_S1_S12[i1][i2] = new double* [fN2];
      for( int i3 = 0; i3 < fN2; ++i3 ){
	fComp_S1_S12[i1][i2][i3] = new double [fN2+1];
      }
    }
  }
  
  fComp_S2_S21 = new double*** [fN2];
  for( int i1 = 0; i1 < fN2; ++i1 ){
    fComp_S2_S21[i1] = new double** [fN2+1];
    for( int i2 = 0; i2 <= fN2; ++i2 ){
      fComp_S2_S21[i1][i2] = new double* [fN1];
      for( int i3 = 0; i3 < fN1; ++i3 ){
	fComp_S2_S21[i1][i2][i3] = new double [fN1+1];
      }
    }
  }

  fComp_S21_S1 = new double*** [fN1]; // 引数の与え方が不自然（要検討）
  for( int i1 = 0; i1 < fN1; ++i1 ){
    fComp_S21_S1[i1] = new double** [fN1+1];
    for( int i2 = 0; i2 <= fN1; ++i2 ){
      fComp_S21_S1[i1][i2] = new double* [fN2];
      for( int i3 = 0; i3 < fN2; ++i3 ){
	fComp_S21_S1[i1][i2][i3] = new double [fN2+1];
      }
    }
  }

  fComp_S12_S2 = new double*** [fN2];
  for( int i1 = 0; i1 < fN2; ++i1 ){
    fComp_S12_S2[i1] = new double** [fN2+1];
    for( int i2 = 0; i2 <= fN2; ++i2 ){
      fComp_S12_S2[i1][i2] = new double* [fN1];
      for( int i3 = 0; i3 < fN1; ++i3 ){
	fComp_S12_S2[i1][i2][i3] = new double [fN1+1];
      }
    }
  }

  fComp_S12_JJ = new double*** [fN1];
  for( int i1 = 0; i1 < fN1; ++i1 ){
    fComp_S12_JJ[i1] = new double** [fN1+1];
    for( int i2 = 0; i2 <= fN1; ++i2 ){
      fComp_S12_JJ[i1][i2] = new double* [fN2];
      for( int i3 = 0; i3 < fN2; ++i3 ){
	fComp_S12_JJ[i1][i2][i3] = new double [fN2+1];
      }
    }
  }

  fComp_S12_C = new double*** [fN1];
  for( int i1 = 0; i1 < fN1; ++i1 ){
    fComp_S12_C[i1] = new double** [fN1+1];
    for( int i2 = 0; i2 <= fN1; ++i2 ){
      fComp_S12_C[i1][i2] = new double* [fN2];
      for( int i3 = 0; i3 < fN2; ++i3 ){
	fComp_S12_C[i1][i2][i3] = new double [fN2+1];
      }
    }
  }

  fComp_S21_C = new double*** [fN2];
  for( int i1 = 0; i1 < fN2; ++i1 ){
    fComp_S21_C[i1] = new double** [fN2+1];
    for( int i2 = 0; i2 <= fN2; ++i2 ){
      fComp_S21_C[i1][i2] = new double* [fN1];
      for( int i3 = 0; i3 < fN1; ++i3 ){
	fComp_S21_C[i1][i2][i3] = new double [fN1+1];
      }
    }
  }

  int N = this->Max(fN1,fN2);
  
  fComp2_1 = new double* [N+1];  
  for( int i = 0; i <= N; ++i )
    fComp2_1[i] = new double [N+1];
  fComp2_5 = new double* [N+1];  
  for( int i = 0; i <= N; ++i )
    fComp2_5[i] = new double [N+1];
}
				 
TSearch_I_get::~TSearch_I_get()
{

}

void TSearch_I_get::SetParameter()
{
  fNum_top = 10000000;  // Get Conf
  TOperator_Candidate::SetParameter();
}

void TSearch_I_get::Define()
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

  // Update_Conf_top
  eval_tmp = new double [ fNum_top ];
  IH_tmp = new int [ fNum_top ];
  nv_tmp = new int [ fNum_top ];
  na_tmp = new int [ fNum_top ];
  st0_tmp = new int [ fNum_top ];
  st1_tmp = new int [ fNum_top ];
  st2_tmp = new int [ fNum_top ];
  fn_tmp = new int* [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    fn_tmp[i] = new int [ 6 ];
  reverse_tmp = new int [ fNum_top ];


  fEval_chunk = new double [ fNum_top ];
  fIH_chunk = new int [ fNum_top ];
  fNv_chunk = new int [ fNum_top ];
  fNa_chunk = new int [ fNum_top ];
  fSt0_chunk = new int [ fNum_top ];
  fSt1_chunk = new int [ fNum_top ];
  fSt2_chunk = new int [ fNum_top ];
  ffn_chunk = new int* [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    ffn_chunk[i] = new int [ 6 ];
  fReverse_chunk = new int [ fNum_top ];
  fIndex_chunk = new int [ fNum_top ];
  fIndex_new = new int [ fNum_top ];

  fn = new int [6];
}

void TSearch_I_get::SetInit( VectorXd w1, VectorXd w2 )
{
  TOperator_Candidate::SetInit( w1, w2 );
  fNum_threads = 1;
}


// All types (fReflection is specified)
void TSearch_I_get::DoIt_IH456()
{
  fEval_best = 999999999.9;
  for( int i = 0; i < fNum_top; ++i ){
    fEval_top[i] = 999999999.9;
    fIndex_top[i] = i;
  }
  fCount_update = 0;
  
  for ( int iter = 0; iter <= 1; ++iter ){ 
      
    if( fReflection == 0 ){
      if( iter == 0 ){
	fReverse = 0;
	fReverse1 = 0; 
	fReverse2 = 0;
      }
      else{
	fReverse = 1;
	fReverse1 = 0;
	fReverse2 = 1;
      }
    }

    if( fReflection == 1 ){
      if( iter == 0 ){
	fReverse = 0;
	fReverse2 = 0;
	fReverse1 = 0;
      }
      else{
	fReverse = 1;
	fReverse2 = 0;
	fReverse1 = 1;
      }
    }

    Set_eval_JS();
    Set_Comp();

#pragma omp parallel  num_threads(fNum_threads)
    {

#pragma omp sections nowait
      {
#pragma omp section
	this->IH4_type1();
#pragma omp section
	this->IH4_type2();
#pragma omp section
	this->IH4_type5();
#pragma omp section
	this->IH4_type9();
#pragma omp section
	this->IH5_type15();
#pragma omp section    
	this->IH6_type5();
#pragma omp section
	this->IH6_type7();
#pragma omp section
	this->IH6_type10();
#pragma omp section
	this->IH6_type1();
#pragma omp section      
	this->IH6_type9();
      }

      
#pragma omp sections 
      {
#pragma omp section
	this->IH5_type2();
#pragma omp section
	this->IH5_type1();
#pragma omp section
	this->IH6_type8();
#pragma omp section
	this->IH6_type12();
#pragma omp section
	this->IH6_type15();
      }

    }
  }
  
  this->Update_Conf_top( 99999999.9, -1, -1, -1 ); // Get Conf
}


void TSearch_I_get::DoIt_IH123()
{
  assert( fN1 == fN2 );

  fEval_best = 999999999.9;
  for( int i = 0; i < fNum_top; ++i ){
    fEval_top[i] = 999999999.9;
    fIndex_top[i] = i;
  }
  fCount_update = 0;
  

  for ( int iter = 0; iter <= 1; ++iter ){ 

    if( iter == 0 ){
      fReverse = 0;
      fReverse1 = 0;
      fReverse2 = 0;
    }
    else{
      fReverse = 1;
      fReverse1 = 0;
      fReverse2 = 1;
    }

    this->Set_eval_JS(); 
    this->Set_Comp();    

    this->IH1();
    this->IH2();
    this->IH3();
  }

  this->Update_Conf_top( 99999999.9, -1, -1, -1 ); // Get Conf
}


void TSearch_I_get::Set_Comp()
{
#pragma omp parallel sections
  {
#pragma omp section
    {
      this->Set_Comp1_1();
      this->Set_Comp1_2();
    }
#pragma omp section
    this->Set_Comp1_3();
#pragma omp section
    this->Set_Comp1_4();
#pragma omp section
    this->Set_Comp1_5();
#pragma omp section
    this->Set_Comp1_6();
#pragma omp section
    this->Set_Comp1_7();
#pragma omp section
    this->Set_Comp1_8();
#pragma omp section
    this->Set_Comp1_9();
  }
}

void TSearch_I_get::Set_Comp1_1()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;
  
  for( int j1 = 0; j1 < fN1; ++j1 ){
    // this->Set_fj( j1, 0 );
    for( int k23 = 0; k23 <= fN1; ++k23 ){
      double eval_min = 999999999.9;
      for( int k2 = 0; k2 <= k23; ++k2 ){
	a1 = j1;
	a2 = a1+k2;
	int k3 = k23-k2;;

	eval = this->Eval_S1_edge(a1,k2) + this->Eval_S1_edge(a2,k3);
	if( eval < eval_min )
	  eval_min = eval;
      }
      fComp_S1_S1[index1(a1)][k23] = eval_min;
    }
  }
}

void TSearch_I_get::Set_Comp1_2()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j2 = 0; j2 < fN2; ++j2 ){
    // this->Set_fj( 0, j2 );
    for( int k45 = 0; k45 <= fN2; ++k45 ){
      double eval_min = 999999999.9;
      for( int k4 = 0; k4 <= k45; ++k4 ){
	b4 = j2;
	b5 = b4+k4;
	int k5 = k45-k4;;

	eval = this->Eval_S2_edge(b4,k4) + this->Eval_S2_edge(b5,k5);
	if( eval < eval_min )
	  eval_min = eval;
      }
      fComp_S2_S2[index2(b4)][k45] = eval_min;
    }
  }
}


void TSearch_I_get::Set_Comp1_3()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      for( int s2 = 0; s2 <= fN2; ++s2 ){
	for( int len = 0; len <= fN1; ++len ){
	  eval_min = 99999999.9;
	  for( int k4 = 0; k4 <= len; ++k4 ){

	    a4 = j1;
	    a5 = a4+k4;
	    b_s = j2;
	    int s1_d = len-k4; 

	    eval = 0.0;
	    eval += this->Eval_S1_edge(a4,k4);
	    eval += this->Eval_S12_sprit(a5,s1_d,b_s,s2);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S1_S12[index1(a4)][len][index2(b_s)][s2] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_4()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //this->Set_fj( j1, j2 );
      for( int s1 = 0; s1 <= fN1; ++s1 ){
	for( int len = 0; len <= fN2; ++len ){
	  eval_min = 99999999.9;
	  for( int k2 = 0; k2 <= len; ++k2 ){

	    b1 = j2;
	    b2 = b1+k2;
	    a_s = j1;
	    int s2_d = len-k2; 

	    eval = 0.0;
	    eval += this->Eval_S2_edge(b1,k2);
	    eval += this->Eval_S21_sprit(b2,s2_d,a_s,s1);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S2_S21[index2(b1)][len][index1(a_s)][s1] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_5()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      for( int s2_d = 0; s2_d <= fN2; ++s2_d ){
	for( int len = 0; len <= fN1; ++len ){
	  eval_min = 99999999.9;
	  for( int k3 = 0; k3 <= len; ++k3 ){

	    int s1 = len-k3; 
	    a_s = j1;
	    a2 = a_s+s1;
	    a3 = a2+k3;
	    b_e = j2;
	    b1 = fN2-s2_d+j2;

	    eval = 0.0;
	    eval += this->Eval_S21_sprit(b1,s2_d,a_s,s1);
	    eval += this->Eval_S1_edge(a2,k3);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S21_S1[index1(a3)][len][index2(b_e)][s2_d] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_6()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      for( int s1_d = 0; s1_d <= fN1; ++s1_d ){
	for( int len = 0; len <= fN2; ++len ){
	  eval_min = 99999999.9;
	  for( int k5 = 0; k5 <= len; ++k5 ){

	    int s2 = len-k5; 
	    b_s = j2;
	    b5 = b_s+s2;
	    b0 = b5+k5;
	    a_e = j1;
	    a4 = fN1-s1_d+j1;

	    eval = 0.0;
	    eval += this->Eval_S12_sprit(a4,s1_d,b_s,s2);
	    eval += this->Eval_S2_edge(b5,k5);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S12_S2[index2(b0)][len][index1(a_e)][s1_d] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_7()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      for( int len1 = 0; len1 <= fN1; ++len1 ){
	for( int len2 = 0; len2 <= fN2; ++len2 ){
      
	  eval_min = 99999999.9;
	  int k0_max = this->Min(len1,len2);
	  for( int k0 = 0; k0 <= k0_max; ++k0 ){
	    int s1_d = len1-k0;
	    int s2 = len2-k0;

	    a2 = j1;
	    a_e = a2+s1_d;
	    b_e = j2;
	    b_s = b_e+k0;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a_e,b_e,k0);
	    eval += this->Eval_S12_sprit(a2,s1_d,b_s,s2);

	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S12_JJ[index1(a2)][len1][index2(b_e)][len2] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_8()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){  
      //  this->Set_fj( j1, j2 );
      for( int len1 = 0; len1 <= fN1; ++len1 ){
	for( int len2 = 0; len2 <= fN2; ++len2 ){
      
	  eval_min = 99999999.9;
	  int k0_max = this->Min(len1,len2);
	  for( int k0 = 0; k0 <= k0_max; ++k0 ){
	    int s1_d = len1-k0;
	    int s2 = len2-k0;

	    a2 = j1;
	    a_e = a2+s1_d;

	    b_e = j2;
	    b_s = b_e+k0;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a_e,b_e,k0);
	    eval += this->Eval_S12_sprit(a2,s1_d,b_s,s2);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S12_C[index1(a2)][len1][index2(b_e)][len2] = eval_min;
	}
      }
    }
  }
}

void TSearch_I_get::Set_Comp1_9()
{
  double eval, evalA, eval_min;
  int a0, a1, a2, a3, a4, a5;
  int b0, b1, b2, b3, b4, b5;
  int a_s, a_e, b_s, b_e;
  int len_m;

  for( int j2 = 0; j2 < fN2; ++j2 ){  
    for( int j1 = 0; j1 < fN1; ++j1 ){
      // this->Set_fj( j1, j2 );
      for( int len2 = 0; len2 <= fN2; ++len2 ){
	for( int len1 = 0; len1 <= fN1; ++len1 ){
      
	  eval_min = 99999999.9;
	  int k0_max = this->Min(len1,len2);
	  for( int k0 = 0; k0 <= k0_max; ++k0 ){
	    int s2_d = len2-k0;
	    int s1 = len1-k0;

	    a_e = j1;
	    a_s = a_e+k0;

	    b4 = j2;
	    b_e = b4+s2_d;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a_e,b_e,k0);
	    eval += this->Eval_S21_sprit(b4,s2_d,a_s,s1);
	    if( eval < eval_min )
	      eval_min = eval;
	  }
	  fComp_S21_C[index2(b4)][len2][index1(a_e)][len1] = eval_min;
	}
      }
    }
  }
}



// IH_6_type1
void TSearch_I_get::Set_Comp2_1( int j1, int j2 )
{
  int a1, a2;
  int b1, b2, b3, b4;
  double eval, eval_min;
  
  int k2_max = this->Min(fN1,fN2);
  for( int k2 = 0; k2 <= k2_max; ++k2 ){
    for( int len = k2; len <= fN2; ++len ){

      eval_min = 99999999.9;
      for( int k3 = 0; k3 <= len-k2; ++k3 ){
	int k4 = len-(k2+k3);
	a1 = j1;
	a2 = a1+k2;
	b1 = j2;
	b2 = b1+k3;
	b3 = b2+k2;
	b4 = b3+k4;

	eval = 0.0;
	eval += this->Eval_J12_GR_edge(a1,b2,k2);
	eval += this->Eval_S2_edge(b1,k3);
	eval += this->Eval_S2_edge(b3,k4);
	if( eval < eval_min )
	  eval_min = eval;
      }
      fComp2_1[k2][len] = eval_min;
    }
  }
}

// IH_6_type5
void TSearch_I_get::Set_Comp2_5( int j1, int j2 )
{
  int a1, a2, a3;
  int b1, b2, b3;
  double eval, eval_min;
  
  int N = this->Min(fN1,fN2);
  for( int len1 = 0; len1 <= N; ++len1 ){
    for( int len2 = 0; len2 <= N; ++len2 ){

      eval_min = 99999999.9;
      int k1_max = this->Min(len1,len2);
      for( int k1 = 0; k1 <= k1_max; ++k1 ){
	int k3 = len1-k1;
	int k4 = len2-k1;
	
	a1 = j1;
	a2 = a1+k1;
	a3 = a1+k3;

	b1 = j2;
	b2 = b1+k4;
	b3 = b2+k1;

	eval = 0.0;
	eval += this->Eval_J12_GR_edge(a1,b2,k1);
	eval += this->Eval_S1_edge(a2,k3);
	eval += this->Eval_S2_edge(b1,k4);
	if( eval < eval_min )
	  eval_min = eval;
      }
      fComp2_5[len1][len2] = eval_min;
    }
  }
}


double TSearch_I_get::Eval_J1_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J1_edge[index1(pos1)][index1(pos2)][len];
  return eval;
}

double TSearch_I_get::Eval_J2_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J2_edge[index2(pos1)][index2(pos2)][len];
  return eval;
}

double TSearch_I_get::Eval_S1_edge( int pos, int len )
{
  double eval;
  eval = fEval_S1_edge[index1(pos)][len];
  return eval;
}

double TSearch_I_get::Eval_S2_edge( int pos, int len )
{
  double eval;
  eval = fEval_S2_edge[index2(pos)][len];
  return eval;
}

double TSearch_I_get::Eval_J12_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J12_edge[index1(pos1)][index2(pos2)][len];
  return eval;
}

double TSearch_I_get::Eval_J21_sprit1( int b2, int s2_d, int a_s, int s1, int b1 )
{
  int b_e_comp;
  double eval;

  b_e_comp = b1+s1;
  
  eval = 0.0;
  eval += this->Eval_J12_edge(a_s, b1, s1);
  eval += this->Eval_J2_edge(b2, b_e_comp, s2_d);
  return eval;
}

double TSearch_I_get::Eval_J12_sprit1( int a1, int s1_d, int b_s, int s2, int b1 )
{
  int b_s_comp;
  double eval;

  b_s_comp = b1+s2;
  
  eval = 0.0;
  eval += this->Eval_J12_edge(a1, b_s_comp, s1_d);
  eval += this->Eval_J2_edge(b_s, b1, s2);
  return eval;
}

double TSearch_I_get::Eval_J12_sprit2( int a_s, int s1, int b3, int s2_d, int a3, int s1_d, int b_s, int s2 )
{
  int a_s_comp, a_e_comp, b_s_comp, b_e_comp;
  int len_m;
  double eval;

  eval = 0.0;
  if( s1 >= s1_d ){
    len_m = s1-s1_d;
    a_e_comp = a_s+len_m;
    b_e_comp = b_s+len_m;
    eval += this->Eval_J2_edge(b3, b_e_comp, s2_d);
    eval += this->Eval_J12_edge(a_s, b_s, len_m);
    eval += this->Eval_J1_edge(a_e_comp, a3, s1_d);
  }
  else{
    len_m = s1_d-s1;
    a_s_comp = a3+s1;
    b_s_comp = b3+s2;
    eval += this->Eval_J2_edge(b3, b_s, s2);
    eval += this->Eval_J12_edge(a_s_comp, b_s_comp, len_m);
    eval += this->Eval_J1_edge(a_s, a3, s1);
  }
  return eval;
}

double TSearch_I_get::Eval_S21_sprit( int b3, int s2_d, int a_s, int s1 )
{
  int len_m;
  double eval;

  eval = 0.0;
  if( s2_d > s1 ){
    len_m = s2_d-s1;
    eval += this->Eval_J12_edge(a_s, b3, s1);
    eval += this->Eval_S2_edge(b3+s1, len_m);
  }
  else{
    len_m = s1-s2_d;
    eval += this->Eval_J12_edge(a_s+len_m, b3, s2_d);
    eval += this->Eval_S1_edge(a_s, len_m);
  }
  return eval;
}

double TSearch_I_get::Eval_S12_sprit( int a3, int s1_d, int b_s, int s2 )
{
  int len_m;
  double eval;

  eval = 0.0;
  if( s1_d > s2 ){
    len_m = s1_d-s2;
    eval += this->Eval_J12_edge(a3, b_s, s2);
    eval += this->Eval_S1_edge(a3+s2, len_m);
  }
  else{
    len_m = s2-s1_d;
    eval += this->Eval_J12_edge(a3, b_s+len_m, s1_d);
    eval += this->Eval_S2_edge(b_s, len_m);
  }
  return eval;
}

// GR
double TSearch_I_get::Eval_J1_GR_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J1_GR_edge[index1(pos1)][index1(pos2)][len];
  return eval;
}

double TSearch_I_get::Eval_J2_GR_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J2_GR_edge[index2(pos1)][index2(pos2)][len];
  return eval;
}

double TSearch_I_get::Eval_J12_GR_edge( int pos1, int pos2, int len )
{
  double eval;
  eval = fEval_J12_GR_edge[index1(pos1)][index2(pos2)][len];
  return eval;
}

// 引数の順番は変更を検討
double TSearch_I_get::Eval_J12_GR_sprit2( int a_s, int s1, int b4, int s2_d, int a2, int s1_d, int b_s, int s2, int a_e, int b_e )
{
  int a_s_comp, a_e_comp, b_s_comp, b_e_comp;
  int len_m;
  double eval;

  eval = 0.0;
  if( s2 <= s1 ){
    len_m = s1-s2;
    a_e_comp = a_s+len_m;
    a_s_comp = a_e-len_m;
    eval += this->Eval_J12_GR_edge(a2, b4, s2_d);
    eval += this->Eval_J1_GR_edge(a_s, a_s_comp, len_m);
    eval += this->Eval_J12_GR_edge(a_e_comp, b_s, s2);
  }
  else{
    len_m = s2-s1;
    b_s_comp = b_e-len_m;
    b_e_comp = b_s+len_m;
    eval += this->Eval_J12_GR_edge(a2, b4, s1_d);
    eval += this->Eval_J2_GR_edge(b_s, b_s_comp, len_m);
    eval += this->Eval_J12_GR_edge(a_s, b_e_comp, s1);
  }
  
  return eval;
}

double TSearch_I_get::Eval_J21_2_GR_sprit1( int b2, int s2_d, int a_s, int s1, int b1 )
{
  int b_e_comp;
  double eval;

  b_e_comp = b1+s2_d;
  
  eval = 0.0;
  eval += this->Eval_J2_GR_edge(b2, b1, s2_d);
  eval += this->Eval_J12_GR_edge(a_s, b_e_comp, s1);
  return eval;
}

double TSearch_I_get::Eval_J12_2_GR_sprit1( int a1, int s1_d, int b_s, int s2, int b1 )
{
  int b_s_comp;
  double eval;

  b_s_comp = b1+s1_d;
  
  eval = 0.0;
  eval += this->Eval_J12_GR_edge(a1, b1, s1_d);
  eval += this->Eval_J2_GR_edge(b_s, b_s_comp, s2);
  return eval;
}

double TSearch_I_get::Eval_J21_1_GR_sprit1( int b2, int s2_d, int a_s, int s1, int a1 )
{
  int a_s_comp;
  double eval;

  a_s_comp = a1+s2_d;
  
  eval = 0.0;
  eval += this->Eval_J12_GR_edge(a1, b2, s2_d);
  eval += this->Eval_J1_GR_edge(a_s, a_s_comp, s1);
  return eval;
}

double TSearch_I_get::Eval_J12_1_GR_sprit1( int a1, int s1_d, int b_s, int s2, int a2 )
{
  int a_e_comp;
  double eval;

  a_e_comp = a2+s1_d;
  
  eval = 0.0;
  eval += this->Eval_J12_GR_edge(a_e_comp, b_s, s2);
  eval += this->Eval_J1_GR_edge(a1, a2, s1_d);
  return eval;
}


//////// IH4 ////////

void TSearch_I_get::IH4_type1()
{
  printf( "IH4_1 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 4;
  fNv = 6;
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max ; ++k0 ){
	int k1_max = (fN1+fN2)/2-k0;
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_max = this->Min(fN1-k0, k1);
	  for( int s1 = 0; s1 <= s1_max; ++s1 ){
	    int s2_max = this->Min(fN2-k0, k1);
	    for( int s2 = 0; s2 <= s2_max; ++s2 ){

	      int s2_d = k1-s1;
	      int s1_d = k1-s2;

	      int k23 = fN1-(k0+s1+s1_d);
	      int k45 = fN2-(k0+s2+s2_d);
	      if( k23 < 0 || k45 < 0 )
		continue;

	      a1 = j1;
	      a_s = fN1-s1+j1;
	      a_e = a_s-k0;
	      a3 = a_e-s1_d;

	      b4 = j2;
	      b_s = fN2-s2+j2;
	      b_e = b_s-k0;
	      b0 = b_e-s2_d;
	      
	      eval = 0.0;
	      eval += this->Eval_J12_edge( a_e, b_e, k0 );
	      eval += this->Eval_J12_sprit2( a_s, s1, b0, s2_d, a3, s1_d, b_s, s2 );
	      evalA = eval;

	      // Prune
	      if( evalA + fComp_S1_S1[index1(a1)][k23] + fComp_S2_S2[index2(b4)][k45] >= fEval_best )
		continue;

	      for( int k2 = 0; k2 <= k23; ++k2 ){
		int k3 = k23-k2;
		a2 = a1+k2;

		eval = evalA;
		eval += this->Eval_S1_edge(a1, k2);
		eval += this->Eval_S1_edge(a2, k3);
		evalB = eval;

		// assert( eval-evalA > fComp_S1_S1[index1(a1)][k23]-0.000001 );

		if( evalB + fComp_S2_S2[index2(b4)][k45] >= fEval_best )
		  continue;
		
		for( int k4 = 0; k4 <= k45; ++k4 ){
		  int k5 = k45-k4;
		  b5 = b4+k4;

		  eval = evalB;
		  eval += this->Eval_S2_edge(b4, k4);
		  eval += this->Eval_S2_edge(b5, k5);

		  // assert( eval-evalB > fComp_S2_S2[index2(b4)][k45]-0.000001 );
		  
		  assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 && 0 <= k5 );
		  assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		  assert( s1+s2_d == k1 );
		  assert( s2+s1_d == k1 );
		  assert( k0+k2+k3+s1+s1_d == fN1 );
		  assert( k0+k4+k5+s2+s2_d == fN2 );
		  if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 || k5 == 0 )
		    continue;
		  
		  ++count;

		  if( eval < fEval_best ){
#pragma omp critical (update)
		    {
		      fIH = 4;
		      fNv = 6;
		      fNa=k0-1;
		      if( fReflection == 0 ){
			fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k4-1; fn[5]=k5-1;
			s0 = s2_d;

			this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		      }
		      else{
			fn[0]=k1-1; fn[1]=k4-1; fn[2]=k5-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
			s0 = s1_d;
			
			this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH4_1: count = %ld time = %lf\n", count, end_t-start_t );
}

void TSearch_I_get::IH4_type2()
{
  printf( "IH4_2 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 4;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );

      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max ; ++k0 ){
	int k1_max = this->Min(fN1-k0,fN2-k0);
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_d_max = fN1-(k0+k1);
	  for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){
	    int s2_d_max = fN2-(k0+k1);
	    for( int s2_d = 0; s2_d <= s2_d_max; ++s2_d ){

	      a_s = j1;
	      a_e = fN1-k0+j1;
	      a4 = a_e-s1_d;
	      a3 = a4-k1;

	      b_s = j2;
	      b_e = fN2-k0+j2;
	      b1 = b_e-s2_d;
	      b0 = b1-k1;

	      eval = 0.0;
	      eval += this->Eval_J12_edge(a_e, b_e, k0);
	      eval += this->Eval_J12_edge(a3, b0, k1);
	      evalA = eval;

	      // prune
	      if( evalA + fComp_S12_S2[index2(b0)][b0-b_s][index1(a_e)][s1_d] + fComp_S21_S1[index1(a3)][a3-a_s][index2(b_e)][s2_d] >= fEval_best )
		continue;

	      int k5_max = fN2-(k0+k1+s2_d);
	      for( int k5 = 0; k5 <= k5_max; ++k5 ){
		int s2 = fN2-(k0+k1+k5+s2_d);
		int k4 = s2+s1_d;
		b5 = b_s+s2;

		eval = evalA;
		eval += this->Eval_S2_edge( b5, k5 );
		eval += this->Eval_S12_sprit( a4, s1_d, b_s, s2 );
		evalB = eval;

		// assert( eval - evalA > fComp_S12_S2[index2(b0)][b0-b_s][index1(a_e)][s1_d] - 0.000001 );
		
		// prune
		if( evalB + fComp_S21_S1[index1(a3)][a3-a_s][index2(b_e)][s2_d] >= fEval_best )
		  continue;

		int k3_max = fN1-(k0+k1+s1_d);
		for( int k3 = 0; k3 <= k3_max; ++k3 ){
		  int s1 = fN1-(k0+k1+k3+s1_d);
		  int k2 = s1+s2_d;
		  a2 = a_s+s1;

		  eval = evalB;
		  eval += this->Eval_S1_edge(a2, k3);
		  eval += this->Eval_S21_sprit( b1, s2_d, a_s, s1 );

		  // assert( eval - evalB > fComp_S21_S1[index1(a3)][a3-a_s][index2(b_e)][s2_d] - 0.000001 );
		  
		  assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 && 0 <= k5 );
		  assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		  assert( s1+s2_d == k2 );
		  assert( s2+s1_d == k4 );
		  assert( k0+k1+k3+s1+s1_d == fN1 );
		  assert( k0+k1+k5+s2+s2_d == fN2 );
		  if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 || k5 == 0 )
		    continue;

		  ++count;

		  if( eval < fEval_best ){
#pragma omp critical (update)
		    {
		      fIH = 4;
		      fNv = 6;
		      fNa=k0-1;
		      if( fReflection == 0 ){
			fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k4-1; fn[5]=k5-1;
			s0 = k1+s2_d;
			this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		      }
		      else{
			fn[0]=k1-1; fn[1]=k4-1; fn[2]=k5-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
			s0 = k1+s1_d;
			this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH4_2: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH4_type5()
{
  printf( "IH4_5 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 4;
  fNv = 6;
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      
      int k1_max = fN2/2;
      for( int k1 = 0; k1 <= k1_max; ++k1 ){
	int s1_max = fN1;
	for( int s1 = 0; s1 <= s1_max; ++s1 ){
	  int s2_d_max = fN2-2*k1;
	  for( int s2_d = 0; s2_d <= s2_d_max; ++s2_d ){
	    int len_max = fN2-(2*k1+s2_d);
	    for( int len = 0; len <= len_max; ++len ){
      
	      int k2 = s1+s2_d;
	      int k45 = fN2-(2*k1+s2_d+len);

	      a_s = j1;
	      a2 = a_s+s1;

	      b_e = j2;
	      b1 = fN2-s2_d+j2;
	      b0 = b1-k1;
	      b3 = b_e+len;
	      b4 = b3+k1; 
	      
	      eval = 0.0;
	      eval += this->Eval_J2_edge(b0, b3, k1);
	      eval += this->Eval_S21_sprit(b1, s2_d, a_s, s1);
	      evalA = eval;
	    
	      if( eval + fComp_S12_JJ[index1(a2)][fN1-s1][index2(b_e)][len] + fComp_S2_S2[index2(b4)][k45] > fEval_best )
		continue;

	      int k0_max = this->Min(fN1-s1, len);
	      for( int k0 = 0; k0 <= k0_max ; ++k0 ){

		int s1_d = fN1-(k0+s1);
		int s2 = len-k0;
		int k3 = s1_d+s2;

		a_e = a2+s1_d;
		b_s = b_e+k0;

		eval = evalA;
		eval += this->Eval_J12_edge(a_e, b_e, k0);
		eval += this->Eval_S12_sprit(a2, s1_d, b_s, s2);
		evalB = eval;

		// assert( eval - evalA > fComp_S12_JJ[index1(a2)][fN1-s1][index2(b_e)][len] - 0.000001 );

		if( eval + fComp_S2_S2[index2(b4)][k45] > fEval_best )
		  continue;

		for( int k4 = 0; k4 <= k45; ++k4 ){
		  int k5 = k45-k4;
		  
		  b1 = b4+k4;

		  eval = evalB;
		  eval += this->Eval_S2_edge(b4, k4);
		  eval += this->Eval_S2_edge(b1, k5);

		  // assert( eval-evalB > fComp_S2_S2[index2(b4)][k45] - 0.000001 );

		  assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 && 0 <= k5 );
		  assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		  assert( s1+s2_d == k2 );
		  assert( s2+s1_d == k3 );
		  assert( k0+s1+s1_d == fN1 );
		  assert( k0+2*k1+k4+k5+s2+s2_d == fN2 );
		  if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 || k5 == 0 )
		    continue;
		  
		  ++count;

		  if( eval < fEval_best ){
#pragma omp critical (update)
		    {
		      fIH = 4;
		      fNv = 6;
		      fNa=k0-1;
		      if( fReflection == 0 ){
			fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k4-1; fn[5]=k5-1;
			s0 = k1+s2_d;
			this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		      }
		      else{
			fn[0]=k1-1; fn[1]=k4-1; fn[2]=k5-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
			s0 = 2*k1+k4+k5+k2+s1_d;
			this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH4_5: count = %ld time = %lf\n", count, end_t-start_t );
}

void TSearch_I_get::IH4_type9()
{
  printf( "IH4_9 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 4;
  fNv = 6;
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max ; ++k0 ){
	int k1_max = this->Min(fN1-k0,fN2-k0);
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_max = fN1-(k0+k1);
	  for( int s1 = 0; s1 <= s1_max; ++s1 ){
	    int s1_d = fN1-(k0+k1+s1);
	    int len2_max = fN2-(k0+k1);
	    for( int len2 = 0; len2 <= len2_max; ++len2 ){
	      int len2_d = fN2-(k0+k1+len2);
	      
	      a_s = j1;
	      a_e = fN1-k0+j1;
	      a0 = a_s+s1;
	      a1 = a0+k1;
	      
	      b_s = j2;
	      b_e = fN2-k0+j2;
	      b3 = b_s+len2;
	      b4 = b3+k1;

	      eval = 0.0;
	      eval += this->Eval_J12_edge(a_e, b_e, k0);
	      eval += this->Eval_J12_edge(a0, b3, k1);
	      evalA = eval;

	      if( evalA + fComp_S12_S2[index2(b3)][len2][index1(a_e)][s1_d] + fComp_S2_S21[index2(b4)][len2_d][index1(a_s)][s1] > fEval_best )
		continue;

	      for( int k3 = 0; k3 <= len2; ++k3 ){
		int s2 = len2-k3;
		b2 = b_s+s2;
		int k2 = s1_d+s2;

		eval = evalA;
		eval += this->Eval_S2_edge(b2, k3);
		eval += this->Eval_S12_sprit( a1, s1_d, b_s, s2 );
		evalB = eval;

		// assert( eval - evalA > fComp_S12_S2[index2(b3)][len2][index1(a_e)][s1_d] - 0.000001 ) ;

		if ( evalB + fComp_S2_S21[index2(b4)][len2_d][index1(a_s)][s1] > fEval_best )
		  continue;

		for( int k4 = 0; k4 <= len2_d; ++k4 ){
		  int s2_d = len2_d-k4;
		  
		  // duplications are avoided
		  /*
		  if( s1_d > s2 && s1 < s2_d )
		    continue;
		  if( s1 < s2_d && s1_d < s2 )
		    if( k0 < k1 || k0 < s1 || k0 < s1_d )
		      continue;
		  */

		  b5 = b4+k4; 
		  int k5 = s1+s2_d;

		  eval = evalB;
		  eval += this->Eval_S2_edge(b4, k4);
		  eval += this->Eval_S21_sprit( b5, s2_d, a_s, s1 );

		  // assert( eval - evalB > fComp_S2_S21[index2(b4)][len2_d][index1(a_s)][s1] - 0.000001 ) ;

		  assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 && 0 <= k5 );
		  assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		  assert( s1+s2_d == k5 );
		  assert( s2+s1_d == k2 );
		  assert( k0+k1+s1+s1_d == fN1 );
		  assert( k0+k1+k3+k4+s2+s2_d == fN2 );
		  if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 || k5 == 0 )
		    continue;
		  
		  ++count;


		  if( eval < fEval_best ){
#pragma omp critical (update)
		    {
		      fIH = 4;
		      fNv = 6;
		      fNa=k0-1;
		      if( fReflection == 0 ){
			fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k4-1; fn[5]=k5-1;
			s0 = 2*k1+k2+k3+k4+s2_d;
			this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		      }
		      else{
			fn[0]=k1-1; fn[1]=k4-1; fn[2]=k5-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
			s0 = 2*k1+k4+k5+s1_d;
			this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH4_9: count = %ld time = %lf\n", count, end_t-start_t );
}


////// IH6 //////

void TSearch_I_get::IH6_type1()
{
  printf( "IH6_1 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      this->Set_Comp2_1( j1, j2 ); 
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max ; ++k0 ){
	int k2_max = this->Min(fN1-k0,fN2-k0);
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int k1_max = (fN1+fN2)/2-(k0+k2);
	  for( int k1 = 0; k1 <= k1_max; ++k1 ){
	    int s1_max = this->Min(fN1-k0-k2, k1);
	    for( int s1 = 0; s1 <= s1_max; ++s1 ){

	      int s2_d = k1-s1;
	      int s1_d = fN1-(k0+k2+s1);
	      int s2 = k1-s1_d;
	      if( s1 < 0 || s1_d < 0 || s2 < 0 || s2_d < 0 )
		continue;

	      a1 = j1;
	      a_s = fN1-s1+j1;
	      a_e = a_s-k0;
	      a2 = a_e-s1_d;

	      b1 = j2;
	      b_s = fN2-s2+j2;
	      b_e = b_s-k0;
	      b4 = b_e-s2_d;

	      eval = 0.0;
	      eval += this->Eval_J12_edge( a_e, b_e, k0 );
	      eval += this->Eval_J12_GR_sprit2( a_s, s1, b4, s2_d, a2, s1_d, b_s, s2, a_e, b_e );
	      evalA = eval;

	      if( evalA + fComp2_1[k2][b4-b1] > fEval_best )
		continue;

	      for( int k3 = 0; k3 <= fN2-(k0+k2+s2+s2_d); ++k3 ){

		b2 = b1+k3;
		b3 = b2+k2;

		int k4 = fN2-(k0+k2+k3+s2+s2_d);

		eval = evalA;
		eval += this->Eval_J12_GR_edge(a1,b2,k2);
		eval += this->Eval_S2_edge(b1,k3);
		eval += this->Eval_S2_edge(b3,k4);

		// assert( eval - evalA > fComp2_1[k2][b4-b1] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k1 );
		assert( s2+s1_d == k1 );
		assert( k0+k2+s1+s1_d == fN1 );
		assert( k0+k2+k3+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;

		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = k1+k2+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_1: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH6_type5()
{
  printf( "IH6_5 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      this->Set_Comp2_5( j1, j2 ); 
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max ; ++k0 ){
	int k2_max = this->Min((fN1+fN2)/2-k0, fN2-k0);
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int s1_max = this->Min(fN1-k0,k2);
	  for( int s1 = 0; s1 <= s1_max; ++s1 ){
	    int s2_max = this->Min(fN2-k0, k2);
	    for( int s2 = 0; s2 <= s2_max; ++s2 ){

	      int s1_d = k2-s2;
	      int s2_d = k2-s1;

	      int k1_max = this->Min(fN1-(k0+s1+s1_d),fN2-(k0+s2+s2_d));
	      if( k1_max < 0 )
		continue;

	      a1 = j1;
	      a_s = fN1-s1+j1;
	      a_e = a_s-k0;
	      a3 = a_e-s1_d;

	      b1 = j2;
	      b_s = fN2-s2+j2;
	      b_e = b_s-k0;
	      b3 = b_e-s2_d;

	      eval = 0.0;
	      eval += this->Eval_J12_edge(a_e, b_e, k0);
	      eval += this->Eval_J12_GR_sprit2( a_s, s1, b3, s2_d, a3, s1_d, b_s, s2, a_e, b_e );
	      evalA = eval;

	      if( evalA + fComp2_5[a3-a1][b3-b1] > fEval_best )
		continue;

	      for( int k1 = 0; k1 <= k1_max; ++k1 ){

		a2 = a1+k1;
		b2 = b3-k1;
		int k3 = a3-a2;
		int k4 = b2-b1;

		eval = evalA;
		eval += this->Eval_J12_GR_edge(a1, b2, k1);
		eval += this->Eval_S1_edge(a2, k3);
		eval += this->Eval_S2_edge(b1, k4);

		// assert( eval - evalA > fComp2_5[a3-a1][b3-b1] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k2 );
		assert( s2+s1_d == k2 );
		assert( k0+k1+k3+s1+s1_d == fN1 );
		assert( k0+k1+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = k1+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+k2+k3+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_5: count = %ld time = %lf\n", count, end_t-start_t );
}
	      
void TSearch_I_get::IH6_type7()
{
  printf( "IH6_7 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );
      int k1_max = fN2/2;
      for( int k1 = 0; k1 <= k1_max ; ++k1 ){
	int k2_max = this->Min(fN1, fN2-2*k1);
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int s2_max = fN2-(2*k1+k2);
	  for( int s2 = 0; s2 <= s2_max; ++s2 ){
	    int s1_d_max = fN1-k2;
	    for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){

	      int k0_max = this->Min(fN1-(k2+s1_d), fN2-(2*k1+k2+s2));
	      if( k0_max < 0 )
		continue;
	      
	      int k4 = s2+s1_d;

	      a_e = j1;
	      a2 = fN1-s1_d+j1;
	      a1 = a2-k2;

	      b4 = j2;
	      b3 = fN2-k1+j2;
	      b2 = b3-k2;
	      b1 = b2-k1;
	      b_s = b1-s2;

	      eval = 0.0;
	      eval += this->Eval_J2_GR_edge( b1, b3, k1 );
	      eval += this->Eval_J12_GR_edge( a1, b2, k2 );
	      eval += this->Eval_S12_sprit( a2, s1_d, b_s, s2 );
	      evalA = eval;

	      if( evalA + fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a1-a_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){

		int s1 = fN1-(k0+k2+s1_d);
		int s2_d = fN2-(k0+2*k1+k2+s2);
		int k3 = s1+s2_d;

		if( s1 < 0 || s2_d < 0 )
		  continue;

		a_s = a_e+k0;
		b_e = b4+s2_d;

		eval = evalA;
		eval += this->Eval_J12_edge(a_e, b_e, k0);
		eval += this->Eval_S21_sprit(b4, s2_d, a_s, s1);

		// assert( eval - evalA > fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a1-a_e] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k3 );
		assert( s2+s1_d == k4 );
		assert( k0+k2+s1+s1_d == fN1 );
		assert( k0+2*k1+k2+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;

		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = 2*k1+k2+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+2*k2+k3+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_7: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH6_type8()
{
  printf( "IH6_8 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k2_max = fN2/2;
      for( int k2 = 0; k2 <= k2_max ; ++k2 ){
	int k1_max = this->Min((fN1+fN2)/2-k2, fN2-2*k2);
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_d_max = this->Min(fN1, k1); 
	  for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){
	    int s2 = k1-s1_d;
	    int k3_max = fN2-(k1+2*k2+s2);
	    for( int k3 = 0; k3 <= k3_max; ++k3 ){

	      int k0_max = this->Min(fN1-s1_d, fN2-(k1+2*k2+k3+s2));
	      if( k0_max < 0 )
		continue;
	      
	      a_e = j1;
	      a1 = fN1-s1_d+j1;

	      b5 = j2;
	      b4 = fN2-k2+j2;
	      b3 = b4-k3;
	      b2 = b3-k1;
	      b1 = b2-k2;
	      b_s = b1-s2;

	      eval = 0.0;
	      eval += this->Eval_J12_2_GR_sprit1( a1, s1_d, b_s, s2, b2 );
	      eval += this->Eval_J2_GR_edge( b1, b4, k2 );
	      eval += this->Eval_S2_edge( b3, k3);
	      evalA = eval;

	      if( evalA + fComp_S21_C[index2(b5)][b_s-b5][index1(a_e)][a1-a_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1 = fN1-(k0+s1_d);
		int s2_d = fN2-(k0+k1+2*k2+k3+s2);
		int k4 = s1+s2_d;

		a_s = a_e+k0;
		b_e = b5+s2_d;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S21_sprit( b5, s2_d, a_s, s1 );
		
		// assert( eval - evalA > fComp_S21_C[index2(b5)][b_s-b5][index1(a_e)][a1-a_e] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k4 );
		assert( s2+s1_d == k1 );
		assert( k0+s1+s1_d == fN1 );
		assert( k0+k1+2*k2+k3+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = 2*k1+2*k2+k3+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();  
  // printf( "IH6_8: count = %ld time = %lf\n", count, end_t-start_t );
}

void TSearch_I_get::IH6_type9()
{
  printf( "IH6_9 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k1_max = this->Min(fN1,fN2);
      for( int k1 = 0; k1 <= k1_max ; ++k1 ){
	int k2_max = this->Min((fN1+fN2)/2-k1, fN2-k1);
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int s1_d_max = this->Min(fN1-k1, k2); 
	  for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){
	    int s2 = k2-s1_d;
	    int k3_max = fN2-(k1+k2+s2);
	    for( int k3 = 0; k3 <= k3_max; ++k3 ){

	      int k0_max = this->Min(fN1-(k1+s1_d), fN2-(k1+k2+k3+s2));
	      if( k0_max < 0 )
		continue;
	      
	      a_e = j1;
	      a2 = fN1-s1_d+j1;
	      a1 = a2-k1; 

	      b4 = j2;
	      b3 = fN2-k2+j2;
	      b2 = b3-k3;
	      b1 = b2-k1;
	      b_s = b1-s2;

	      eval = 0.0;
	      eval += this->Eval_J12_2_GR_sprit1( a2, s1_d, b_s, s2, b3 );
	      eval += this->Eval_J12_GR_edge( a1, b1, k1 );
	      eval += this->Eval_S2_edge( b2, k3 );
	      evalA = eval;

	      if( evalA + fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a1-a_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1 = fN1-(k0+k1+s1_d);
		int s2_d = fN2-(k0+k1+k2+k3+s2);
		int k4 = s1+s2_d;

		a_s = a_e+k0;
		b_e = b4+s2_d;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S21_sprit( b4, s2_d, a_s, s1 );
		
		// assert( eval - evalA > fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a1-a_e] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k4 );
		assert( s2+s1_d == k2 );
		assert( k0+k1+s1+s1_d == fN1 );
		assert( k0+k1+k2+k3+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = 2*k1+2*k2+k3+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = k1+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_9: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH6_type10()
{
  printf( "IH6_10 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k2_max = this->Min(fN1,fN2);
      for( int k2 = 0; k2 <= k2_max ; ++k2 ){
	int k1_max = this->Min((fN1+fN2)/2-k2, fN2-k2);
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_max = this->Min(fN1-k2, k1); 
	  for( int s1 = 0; s1 <= s1_max; ++s1 ){
	    int s2_d = k1-s1;
	    int k3_max = fN1-(k2+s1);
	    for( int k3 = 0; k3 <= k3_max; ++k3 ){

	      int k0_max = this->Min(fN1-(k2+k3+s1), fN2-(k1+k2+s2_d));
	      if( k0_max < 0 )
		continue;

	      a3 = j1;
	      a2 = fN1-k2+j1;
	      a1 = a2-k3;
	      a_s = a1-s1;

	      b_e = j2;
	      b3 = fN2-s2_d+j2; 
	      b2 = b3-k2;
	      b1 = b2-k1;

	      eval = 0.0;
	      eval += this->Eval_J21_2_GR_sprit1( b3, s2_d, a_s, s1, b1 );
	      eval += this->Eval_J12_GR_edge( a2, b2, k2 );
	      eval += this->Eval_S1_edge( a1, k3 );
	      evalA = eval;

	      if( evalA + fComp_S12_C[index1(a3)][a_s-a3][index2(b_e)][b1-b_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1_d = fN1-(k0+k2+k3+s1);
		int s2 = fN2-(k0+k1+k2+s2_d);
		int k4 = s2+s1_d;

		a_e = a3+s1_d;
		b_s = b1-s2;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S12_sprit( a3, s1_d, b_s, s2 );
		
		// assert( eval - evalA > fComp_S12_C[index1(a3)][a_s-a3][index2(b_e)][b1-b_e] - 0.000001 );
		
		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k1 );
		assert( s2+s1_d == k4 );
		assert( k0+k2+k3+s1+s1_d == fN1 );
		assert( k0+k1+k2+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = k1+k2+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+2*k2+k3+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_10: count = %ld time = %lf\n", count, end_t-start_t );
}




void TSearch_I_get::IH6_type12()
{
  printf( "IH6_12 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k2_max = fN2/2;
      for( int k2 = 0; k2 <= k2_max ; ++k2 ){
	int k1_max = this->Min((fN1+fN2)/2-k2, fN2-2*k2);
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int s1_max = this->Min(fN1, k1); 
	  for( int s1 = 0; s1 <= s1_max; ++s1 ){
	    int s2_d = k1-s1;
	    int k4_max = fN2-(k1+2*k2+s2_d);
	    for( int k4 = 0; k4 <= k4_max; ++k4 ){

	      int k0_max = this->Min(fN1-s1, fN2-(k1+2*k2+k4+s2_d));
	      if( k0_max < 0 )
		continue;

	      a1 = j1;
	      a_s = fN1-s1+j1;

	      b_e = j2;
	      b5 = fN2-s2_d+j2;
	      b4 = b5-k2;
	      b3 = b4-k1;
	      b2 = b3-k4;
	      b1 = b2-k2;

	      eval = 0.0;
	      eval += this->Eval_J21_2_GR_sprit1( b5, s2_d, a_s, s1, b3 );
	      eval += this->Eval_J2_GR_edge( b1, b4, k2 );
	      eval += this->Eval_S2_edge( b2, k4 );
	      evalA = eval;

	      if( evalA + fComp_S12_C[index1(a1)][a_s-a1][index2(b_e)][b1-b_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1_d = fN1-(k0+s1);
		int s2 = fN2-(k0+k1+2*k2+k4+s2_d);
		int k3 = s1_d+s2;

		a_e = a_s-k0;
		b_s = b_e+k0;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S12_sprit( a1, s1_d, b_s, s2 );
		
		// assert( eval - evalA > fComp_S12_C[index1(a1)][a_s-a1][index2(b_e)][b1-b_e] - 0.000001 );

		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k1 );
		assert( s2+s1_d == k3 );
		assert( k0+s1+s1_d == fN1 );
		assert( k0+k1+2*k2+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = k1+k2+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+k2+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_12: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH6_type15()
{
  printf( "IH6_15 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a1, a2, a3, a4, a5; 
  int b_s, b_e, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 6;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k1_max = this->fN2/2;
      for( int k1 = 0; k1 <= k1_max ; ++k1 ){
	int k2_max = this->Min((fN1+fN2)/2-k1, fN2-2*k1);
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int s1_d_max = this->Min(fN1, k2); 
	  for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){
	    int s2 = k2-s1_d;
	    int k4_max = fN2-(2*k1+k2+s2);
	    for( int k4 = 0; k4 <= k4_max; ++k4 ){

	      int k0_max = this->Min(fN1-s1_d, fN2-(2*k1+k2+k4+s2));
	      if( k0_max < 0 )
		continue;

	      a_e = j1;
	      a1 = fN1-s1_d+j1;

	      b5 = j2;
	      b4 = fN2-k1+j2;
	      b3 = b4-k2;
	      b2 = b3-k1;
	      b1 = b2-k4;
	      b_s = b1-s2;

	      eval = 0.0;
	      eval += this->Eval_J12_2_GR_sprit1( a1, s1_d, b_s, s2, b3 );
	      eval += this->Eval_J2_GR_edge( b2, b4, k1 );
	      eval += this->Eval_S2_edge( b1, k4 );
	      evalA = eval;

	      if( evalA + fComp_S21_C[index2(b5)][b_s-b5][index1(a_e)][a1-a_e] > fEval_best )
		continue;

	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1 = fN1-(k0+s1_d);
		int s2_d = fN2-(k0+2*k1+k2+k4+s2);
		int k3 = s1+s2_d;

		a_s = a_e+k0;
		b_e = b_s-k0;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S21_sprit( b5, s2_d, a_s, s1 );

		// assert( eval - evalA > fComp_S21_C[index2(b5)][b_s-b5][index1(a_e)][a1-a_e] - 0.000001 );
		
		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k3 );
		assert( s2+s1_d == k2 );
		assert( k0+s1+s1_d == fN1 );
		assert( k0+2*k1+k2+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 6;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k1-1; fn[3]=k3-1; fn[4]=k2-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = 2*k1+k2+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+k2+k3+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH6_15: count = %ld time = %lf\n", count, end_t-start_t );
}

//// IH5 ////
void TSearch_I_get::IH5_type1()
{
  printf( "IH5_1 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5;
  int b_s, b_e, b0, b1, b2, b3, b4, b5;
  double evalA, eval;

  fIH = 5;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max; ++k0 ){
	int k1_max = (fN2-k0)/2;
	for( int k1 = 0; k1 <= k1_max ; ++k1 ){
	  int k2_max = (fN1+fN2)/2-(k0+k1);
	  for( int k2 = 0; k2 <= k2_max; ++k2 ){
	    int s1_max = this->Min(fN1-k0, k2);
	    int s1_min = this->Max(k0+2*k1+k2-fN2, 0); 
	    for( int s1 = s1_min; s1 <= s1_max; ++s1 ){
	      int s2_d = k2-s1;
	      int s1_d = fN1-(k0+s1);
	      int s2 = k2-s1_d;
	      if( s2 < 0 )
		continue;

	      a_s = j1;
	      a_e = fN1-k0+j1;
	      a2 = a_s+s1;

	      b4 = j2;
	      b3 = fN2-k1+j2;
	      b_s = b3-s2;
	      b_e = b_s-k0;
	      b1 = b_e-s2_d;
	      b0 = b1-k1;

	      eval = 0.0;
	      eval += this->Eval_J12_edge( a_e, b_e, k0 );
	      eval += this->Eval_J12_GR_sprit2( a_s, s1, b1, s2_d, a2, s1_d, b_s, s2, a_e, b_e );
	      eval += this->Eval_J2_edge( b3, b0, k1 );
	      evalA = eval;

	      int len = b0-b4;
	      if( evalA + fComp_S2_S2[index2(b4)][len] > fEval_best )
		continue;

	      for( int k3 = 0; k3 <= len; ++k3 ){
		int k4 = len-k3;

		b5 = b4+k3;

		eval = evalA;
		eval += this->Eval_S2_edge( b4, k3 );
		eval += this->Eval_S2_edge( b5, k4 );

		// assert( eval - evalA > fComp_S2_S2[index2(b4)][len] - 0.000001 );
		
		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k2 );
		assert( s2+s1_d == k2 );
		assert( k0+s1+s1_d == fN1 );
		assert( k0+2*k1+k3+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 5;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k2-1; fn[3]=k1-1; fn[4]=k3-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = k1+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = k1+k2+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH5_1: count = %ld time = %lf\n", count, end_t-start_t );
}

void TSearch_I_get::IH5_type2()
{
  printf( "IH5_2 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5;
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, eval;

  fIH = 5;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k0_max = this->Min(fN1,fN2);
      for( int k0 = 0; k0 <= k0_max; ++k0 ){
	int k2_max = (fN2-k0)/2;
	for( int k2 = 0; k2 <= k2_max; ++k2 ){
	  int k1_max = (fN1+fN2)/2-(k0+k2);
	  for( int k1 = 0; k1 <= k1_max; ++k1 ){
	    int s1_max = this->Min(fN1-(k0+2*k2), k1);
	    int s1_min = this->Max(k0+k1-fN2, 0); 
	    for( int s1 = s1_min; s1 <= s1_max; ++s1 ){
	      int s2_d = k1-s1;
	      int s1_d = fN1-(k0+2*k2+s1);
	      int s2 = k1-s1_d;
	      assert( s1_d >= 0 );
	      if( s2 < 0 )
		continue;
	      
	      a_s = j1;
	      a_e = fN1-k0+j1;
	      a1 = a_s+s1;
	      a2 = a1+k2;
	      a3 = a2+k2;

	      b4 = j2;
	      b_s = fN2-s2+j2;
	      b_e = b_s-k0; 
	      b0 = b_e-s2_d;

	      eval = 0.0;
	      eval += this->Eval_J12_edge( a_e, b_e, k0 );
	      eval += this->Eval_J1_GR_edge( a1, a2, k2 );
	      eval += this->Eval_J12_sprit2( a_s, s1, b0, s2_d, a3, s1_d, b_s, s2 );
	      evalA = eval;

	      int len = b0-b4;
	      if( evalA + fComp_S2_S2[index2(b4)][len] > fEval_best )
		continue;

	      for( int k3 = 0; k3 <= len; ++k3 ){
		int k4 = len-k3;

		b5 = b4+k3;

		eval = evalA;
		eval += this->Eval_S2_edge( b4, k3 );
		eval += this->Eval_S2_edge( b5, k4 );

		// assert( eval - evalA > fComp_S2_S2[index2(b4)][len] - 0.000001 );
		
		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k1 );
		assert( s2+s1_d == k1 );
		assert( k0+2*k2+s1+s1_d == fN1 );
		assert( k0+k3+k4+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 )
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 5;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k2-1; fn[3]=k1-1; fn[4]=k3-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = k1+2*k2+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH5_2: count = %ld time = %lf\n", count, end_t-start_t );
}


void TSearch_I_get::IH5_type15()
{
  printf( "IH5_15 (%d)\n", fReflection );
  double start_t = omp_get_wtime();
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5;
  int b_s, b_e, b0, b1, b2, b3, b4, b5;
  double evalA, eval;

  fIH = 5;
  fNv = 6;

  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      //      this->Set_fj( j1, j2 );
      int k1_max = fN2/2;
      for( int k1 = 0; k1 <= k1_max; ++k1 ){
	int k2_max = (fN2-2*k1)/2;
	for( int k2 = 0; k2 <= k2_max ; ++k2 ){
	  int s2_max = fN2-(2*k1+2*k2);
	  for( int s2 = 0; s2 <= s2_max; ++s2 ){
	    int s1_d_max = fN1;
	    for( int s1_d = 0; s1_d <= s1_d_max; ++s1_d ){
	      int k4 = s2+s1_d;

	      a_e = j1;
	      a5 = fN1-s1_d+j1;

	      b4 = j2;
	      b3 = fN2-k1+j2;
	      b2 = b3-k2;
	      b1 = b2-k2;
	      b0 = b1-k1;
	      b_s = b0-s2;

	      eval = 0.0;
	      eval += this->Eval_J2_edge( b0, b3, k1 );
	      eval += this->Eval_J2_GR_edge( b1, b2, k2 );
	      eval += this->Eval_S12_sprit( a5, s1_d, b_s, s2 );
	      evalA = eval;

	      if( evalA + fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a5-a_e] > fEval_best )
		continue;

	      int k0_max = this->Min(fN1-s1_d, fN2-(2*k1+2*k2+s2));
	      for( int k0 = 0; k0 <= k0_max; ++k0 ){
		int s1 = fN1-(k0+s1_d);
		int s2_d = fN2-(k0+2*k1+2*k2+s2);
		int k3 = s1+s2_d;

		a_s = a_e+k0;

		b_e = b_s-k0;

		eval = evalA;
		eval += this->Eval_J12_edge( a_e, b_e, k0 );
		eval += this->Eval_S21_sprit( b4, s2_d, a_s, s1 );
		
		// assert( eval - evalA > fComp_S21_C[index2(b4)][b_s-b4][index1(a_e)][a5-a_e] - 0.000001 );
		
		assert( 0 <= k0 && 0 <= k1 && 0 <= k2 && 0 <= k3 && 0 <= k4 );
		assert( 0 <= s1 && 0 <= s2 && 0 <= s1_d && 0 <= s2_d  );
		assert( s1+s2_d == k3 );
		assert( s2+s1_d == k4 );
		assert( k0+s1+s1_d == fN1 );
		assert( k0+2*k1+2*k2+s2+s2_d == fN2 );
		if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 || k4 == 0 ) 
		  continue;
		
		++count;

		if( eval < fEval_best ){
#pragma omp critical (update)
		  {
		    fIH = 5;
		    fNv = 6;
		    fNa=k0-1;
		    fn[0]=k1-1; fn[1]=k2-1; fn[2]=k2-1; fn[3]=k1-1; fn[4]=k3-1; fn[5]=k4-1;
		    if( fReflection == 0 ){
		      s0 = 2*k1+2*k2+s2_d;
		      this->Update_Conf_top( eval, s0, index1(a_s), index2(b_s) );
		    }
		    else{
		      s0 = 2*k1+2*k2+k3+s1_d;
		      this->Update_Conf_top( eval, s0, index2(b_s), index1(a_s) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double end_t = omp_get_wtime();
  // printf( "IH5_15: count = %ld time = %lf\n", count, end_t-start_t );
}



void TSearch_I_get::IH1() 
{
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 1;
  fNv = 6;

  assert( fN1 == fN2 );
  assert( fReflection == 0 );
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );

      int k0_max = fN1;
      for( int k0 = 0; k0 <= k0_max; ++k0 ){
	int k1_max = fN1-k0;
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int k2_max = fN1-(k0+k1);
	  for( int k2 = 0; k2 <= k2_max; ++k2 ){
	    int k3 = fN1-(k0+k1+k2);
	    a0 = j1;
	    a1 = a0+k1;
	    a2 = a1+k2;
	    a3 = a2+k3;
	    b0 = j2;
	    b1 = b0+k1;
	    b2 = b1+k2;
	    b3 = b2+k3;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a3, b3, k0);
	    eval += this->Eval_J12_edge(a0, b0, k1);
	    eval += this->Eval_J12_edge(a1, b1, k2);
	    eval += this->Eval_J12_edge(a2, b2, k3);

	    assert( 0 <= k0 && 0<= k1 && 0 <= k2 && 0 <= k3 );
	    assert( k0+k1+k2+k3 == fN1 );
	    if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 )
	      continue;

	     ++count;

	     if( eval < fEval_best ){
	       fNa=k0-1;
	       fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
	       s0 = 0;
	       this->Update_Conf_top( eval, s0, index1(a0), index2(b0) );
	     }
	     
	  }
	}
      }
    }
  }
}

void TSearch_I_get::IH2() 
{
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;


  fIH = 2;
  fNv = 6;

  assert( fN1 == fN2 );
  assert( fReflection == 0 );
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );

      int k0_max = fN1;
      for( int k0 = 0; k0 <= k0_max; ++k0 ){
	int k1_max = fN1-k0;
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int k2_max = fN1-(k0+k1);
	  for( int k2 = 0; k2 <= k2_max; ++k2 ){
	    int k3 = fN1-(k0+k1+k2);
	    a0 = j1;
	    a1 = a0+k1;
	    a2 = a1+k2;
	    a3 = a2+k3;
	    b0 = j2;
	    b1 = b0+k3;
	    b2 = b1+k2;
	    b3 = b2+k1;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a3, b3, k0);
	    eval += this->Eval_J12_edge(a1, b1, k2);
	    eval += this->Eval_J12_GR_edge(a2,b0,k3);
	    eval += this->Eval_J12_GR_edge(a0,b2,k1);

	    assert( 0 <= k0 && 0<= k1 && 0 <= k2 && 0 <= k3 );
	    assert( k0+k1+k2+k3 == fN1 );
	    if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 )
	      continue;

	     ++count;

	     if( eval < fEval_best ){
	       fNa=k0-1;
	       fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
	       s0 = 0;
	       this->Update_Conf_top( eval, s0, index1(a0), index2(b0) );
	     }
	     
	  }
	}
      }
    }
  }
}

void TSearch_I_get::IH3() 
{
  int s0;
  long int count;
  int a_s, a_e, a0, a1, a2, a3, a4, a5; 
  int b_s, b_e, b0, b1, b2, b3, b4, b5; 
  double evalA, evalB, eval;

  fIH = 3;
  fNv = 6;

  assert( fN1 == fN2 );
  assert( fReflection == 0 );
  
  count = 0;
  for( int j1 = 0; j1 < fN1; ++j1 ){
    for( int j2 = 0; j2 < fN2; ++j2 ){
      // this->Set_fj( j1, j2 );

      int k0_max = fN1;
      for( int k0 = 0; k0 <= k0_max; ++k0 ){
	int k1_max = fN1-k0;
	for( int k1 = 0; k1 <= k1_max; ++k1 ){
	  int k2_max = fN1-(k0+k1);
	  for( int k2 = 0; k2 <= k2_max; ++k2 ){
	    int k3 = fN1-(k0+k1+k2);
	    a0 = j1;
	    a1 = a0+k1;
	    a2 = a1+k2;
	    a3 = a2+k3;
	    b0 = j2;
	    b1 = b0+k1;
	    b2 = b1+k3;
	    b3 = b2+k2;

	    eval = 0.0;
	    eval += this->Eval_J12_edge(a3, b3, k0);
	    eval += this->Eval_J12_edge(a0, b0, k1);
	    eval += this->Eval_J12_GR_edge(a1,b2,k2);
	    eval += this->Eval_J12_GR_edge(a2,b1,k3);

	    assert( 0 <= k0 && 0<= k1 && 0 <= k2 && 0 <= k3 );
	    assert( k0+k1+k2+k3 == fN1 );
	    if( k0 == 0 || k1 == 0 || k2 == 0 || k3 == 0 )
	      continue;

	     ++count;

	     if( eval < fEval_best ){
	       fNa=k0-1;
	       fn[0]=k1-1; fn[1]=k2-1; fn[2]=k3-1; fn[3]=k1-1; fn[4]=k2-1; fn[5]=k3-1;
	       s0 = 0;
	       this->Update_Conf_top( eval, s0, index1(a0), index2(b0) );
	     }
	     
	  }
	}
      }
    }
  }
}



void TSearch_I_get::Set_candi( int num_candi, int* IH_candi, int* Na_candi, int** fn_candi, int* st0_candi, int* st1_candi, int* st2_candi, int* reverse_candi )
{
  // これを入れとかないとエラー
}


// Get Conf
void TSearch_I_get::Update_Conf_top( double eval, int st0, int st1, int st2 )
{
  ++fCall_Update; // temp
  
  if( st0 != -1 ){ // st0 = -1 は最後の処理フラグ
    int n = fN1+fN2-2*(fNa+1);
    if ( st0 == n )
      st0 = 0;

    fEval_chunk[fCount_update] = eval;
    fIH_chunk[fCount_update] = fIH;
    fNv_chunk[fCount_update] = fNv;
    fNa_chunk[fCount_update] = fNa;
    for( int h = 0; h < fNv; ++h )
      ffn_chunk[fCount_update][h] = fn[h];
    fSt0_chunk[fCount_update] = st0;
    fSt1_chunk[fCount_update] = st1;
    fSt2_chunk[fCount_update] = st2;
    fReverse_chunk[fCount_update] = fReverse;

    ++fCount_update;
  }

  if( fCount_update == 0 ){
    assert( st0 == -1 );
    return;
  }

  if( fCount_update == fNum_top || st0 == -1 ){   // st0 = -1 は最後の処理フラグ
    for( int i = 0; i < fCount_update; ++i )  // required for QuickIndex
      fIndex_chunk[i] = i;    
    this->QuickIndex( fEval_chunk, fIndex_chunk, 0, fCount_update-1 );

    double e_top, e_chunk;    
    int i_top, i_chunk;
    int orderIns_chunk = 0;
    int orderIns_top = 0;
    int orderIns_new = 0;
    int num_used_chunk = 0;
    double eval_tail = -1.0;
    while( 1 ){
      assert( orderIns_top < fNum_top );

      i_top = fIndex_top[orderIns_top];
      e_top = fEval_top[i_top]; 

      if( orderIns_chunk < fCount_update ){
	i_chunk = fIndex_chunk[orderIns_chunk];
	e_chunk = fEval_chunk[i_chunk];
      }
      else{
	e_chunk = e_top + 1.0;
      }
      if( fabs( e_chunk - e_top ) < 0.0000001 ){// この場合、e_chunkは捨てる
	e_chunk = e_top + 1.0;
	++orderIns_chunk;
	continue;
      }

      // duplications are avoided
      if( fabs( e_chunk - eval_tail ) < 0.0000001 ){
	
	++orderIns_chunk;
	continue;
      }

      if( e_top <= e_chunk ){
	fIndex_new[orderIns_new] = i_top;
	++orderIns_top;
      }
      else {
	int i_new = fIndex_top[fNum_top-1-num_used_chunk];
	fIndex_new[orderIns_new] = i_new;
	fEval_top[i_new] = fEval_chunk[i_chunk];
	fIH_top[i_new] = fIH_chunk[i_chunk];
	fNv_top[i_new] = fNv_chunk[i_chunk];
	fNa_top[i_new] = fNa_chunk[i_chunk];
	for( int h = 0; h < fNv; ++h )
	  ffn_top[i_new][h] = ffn_chunk[i_chunk][h];
	fSt0_top[i_new] = fSt0_chunk[i_chunk];
	fSt1_top[i_new] = fSt1_chunk[i_chunk];
	fSt2_top[i_new] = fSt2_chunk[i_chunk];
	fReverse_top[i_new] = fReverse_chunk[i_chunk];
	++orderIns_chunk;
	++num_used_chunk;
      }

      eval_tail = fEval_top[fIndex_new[orderIns_new]];
      ++orderIns_new;
      if( orderIns_new == fNum_top )
	break;
    }

    for( int i = 0; i < fNum_top; ++i )
      fIndex_top[i] = fIndex_new[i];

    fCount_update = 0;
    fEval_best = fEval_top[fIndex_top[fNum_top-1]]; // とりえあえず

    for( int i = 0; i < fNum_top-1; ++i ){
      assert( fEval_top[fIndex_top[i]] <=fEval_top[fIndex_top[i+1]] );
    }
  }
}



// Get Conf
void TSearch_I_get::QuickIndex( double* Arg, int* indexOrderd, int begin, int end )
{
  int i, j, m;
  double pivot;
  int stock; 
  int flag; 

  flag = 0;

  for( m = (begin + end)/2; m < end; ++m ){
    if( Arg[ indexOrderd[m] ] != Arg[ indexOrderd[m+1] ] ){
      if( Arg[ indexOrderd[m] ] > Arg[ indexOrderd[m+1] ] )
	pivot = Arg[ indexOrderd[m] ];
      else
	pivot = Arg[ indexOrderd[m+1] ];
      flag = 1;
      break;
    }
  }
  if( flag == 0 ){
    for( m = (begin + end)/2; m > begin; --m ){
      if( Arg[ indexOrderd[m] ] != Arg[ indexOrderd[m-1] ] ){
	if( Arg[ indexOrderd[m] ] > Arg[ indexOrderd[m-1] ] )
	  pivot = Arg[ indexOrderd[m] ];
	else
	  pivot = Arg[ indexOrderd[m-1] ];
	flag = 1;
	break;
      }
    }
  }

  if( flag == 0 )
    return;


  i = begin;
  j = end;

  while(1)
  {

    while(1){
      if( Arg[ indexOrderd[i] ] >= pivot )
	break;
      ++i;
    }
    while(1){
      if( Arg[ indexOrderd[j] ] < pivot )
	break;
      --j;
    }

    if( i >= j ) 
      break;

    stock = indexOrderd[i];
    indexOrderd[i] = indexOrderd[j];
    indexOrderd[j] = stock;
  }

  if( begin < i-1 )
    this->QuickIndex( Arg, indexOrderd, begin, i-1 );
  if( j+1 < end )
    this->QuickIndex( Arg, indexOrderd, j+1, end );
}


int TSearch_I_get::Min( int a, int b )
{
  if( a < b )
    return a;
  else
    return b;
}

int TSearch_I_get::Max( int a, int b )
{
  if( a < b )
    return b;
  else
    return a;
}
