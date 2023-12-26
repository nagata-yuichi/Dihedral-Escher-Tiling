/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_GET__
#define __SEARCH_I_GET__

#ifndef  __OPERATOR_CANDIDATE__  
#include "operator_candidate.h"  
#endif

#include <omp.h>
#include <unistd.h>

class TSearch_I_get : TOperator_Candidate     
{  
 public:
  TSearch_I_get( int N1, int N2 );
  ~TSearch_I_get(); 

  void SetParameter();  // set parameter values
  void Define();
  void SetInit( VectorXd w1, VectorXd w2 );

  void DoIt_IH456();
  void DoIt_IH123();

  void Set_Comp();
  void Set_Comp1_1();
  void Set_Comp1_2();
  void Set_Comp1_3();
  void Set_Comp1_4();
  void Set_Comp1_5();
  void Set_Comp1_6();
  void Set_Comp1_7();
  void Set_Comp1_8();
  void Set_Comp1_9();

  void Set_Comp2_1( int j1, int j2 );
  void Set_Comp2_5( int j1, int j2 );

  double Eval_J1_edge( int pos1, int pos2, int len );
  double Eval_J2_edge( int pos1, int pos2, int len );
  double Eval_S1_edge( int pos, int len );
  double Eval_S2_edge( int pos, int len );
  double Eval_J12_edge( int pos1, int pos2, int len );
  double Eval_J21_sprit1( int b2, int s2_d, int a_s, int s1, int b1 );
  double Eval_J12_sprit1( int a1, int s1_d, int b_s, int s2, int b1 );
  double Eval_J12_sprit2( int a_s, int s1, int b3, int s2_d, int a3, int s1_d, int b_s, int s2 );
  
  double Eval_S21_sprit( int b3, int s2_d, int a_s, int s1 );
  double Eval_S12_sprit( int a3, int s1_d, int b_s, int s2 );

  double Eval_J1_GR_edge( int pos1, int pos2, int len );
  double Eval_J2_GR_edge( int pos1, int pos2, int len );
  double Eval_J12_GR_edge( int pos1, int pos2, int len );
  double Eval_J12_GR_sprit2( int a_s, int s1, int b4, int s2_d, int a2, int s1_d, int b_s, int s2, int a_e, int b_e );
  double Eval_J21_2_GR_sprit1( int b2, int s2_d, int a_s, int s1, int b1 );
  double Eval_J12_2_GR_sprit1( int a1, int s1_d, int b_s, int s2, int b1 );
  double Eval_J21_1_GR_sprit1( int b2, int s2_d, int a_s, int s1, int b1 );
  double Eval_J12_1_GR_sprit1( int a1, int s1_d, int b_s, int s2, int b1 );


  void IH4_type1();
  void IH4_type2();
  void IH4_type3();
  void IH4_type4();
  void IH4_type5();
  void IH4_type6();
  void IH4_type7();
  void IH4_type8();
  void IH4_type9();

  void IH6_type1();
  void IH6_type2();
  void IH6_type3();
  void IH6_type4();
  void IH6_type5();
  void IH6_type6();
  void IH6_type7();
  void IH6_type8();
  void IH6_type9();
  void IH6_type10();
  void IH6_type11();
  void IH6_type12();
  void IH6_type13();
  void IH6_type14();
  void IH6_type15();

  void IH5_type1();
  void IH5_type2();
  void IH5_type3();
  void IH5_type4();
  void IH5_type5();
  void IH5_type6();
  void IH5_type7();
  void IH5_type8();
  void IH5_type9();
  void IH5_type10();
  void IH5_type11();
  void IH5_type12();
  void IH5_type13();
  void IH5_type14();
  void IH5_type15();

  void IH1();
  void IH2();
  void IH3();


  void Set_candi( int num_candi, int* IH_candi, int* Na_candi, int** fn_candi, int* st0_candi, int* st1_candi, int* st2_candi, int* reverse_candi );
  void Update_Conf_top( double eval, int st0, int st1, int st2 );
  void Update_Sort();
  void QuickIndex( double* Arg, int* indexOrderd, int begin, int end );
  int Min( int a, int b );
  int Max( int a, int b );


  int fNum_top;      // the number of top solutions stored in the EST
  double fEval_best;
  double *fEval_top; // the evaluating values of the top solutions 
  VectorXd *fU_top;  // the top solutions (tile shape)
  VectorXd *fUa_top;  // the top solutions (tile shape)
  VectorXd *fUU1_top; // the top solutions (mesh representation)
  VectorXd *fUU2_top; // the top solutions (mesh representation) 
  int *fIndex_top;   // [i] -> index of the i-th best solution
  int **ffn_top;     // [s] -> fi[] of the tile shape fU_top[s]
  int *fNv_top;      // [s] -> fNv of the tile shape fU_top[s]
  int *fIH_top;      // [s] -> IH type of the tile shape fU_top[s]
  int *fN_top;       //
  int *fNa_top;       // 
  int *fSt0_top;     // [s] -> st0 of the tile shape fU_top[s]
  int *fSt1_top;     // [s] -> st1 of the tile shape fU_top[s]
  int *fSt2_top;     // [s] -> st2 of the tile shape fU_top[s]
  int *fReverse_top; //

  double **fComp_S1_S1; 
  double **fComp_S2_S2; 
  double ****fComp_S1_S12; 
  double ****fComp_S2_S21; 
  double ****fComp_S21_S1; 
  double ****fComp_S12_S2; 
  double ****fComp_S12_JJ; 
  double ****fComp_S12_C;  
  double ****fComp_S21_C;  

  double **fComp2_1;
  double **fComp2_5;

  int fType;         // template type  
  int fReflection;   // templateの反転 （W1とW2を入れ替えに対応）
  int fReverse;      // W1 or W2の反転

  int fCount_update;
  int fFlag_sort;
  int fThread_num_sort;
  double *eval_tmp;
  int *IH_tmp;
  int *nv_tmp;
  int *na_tmp;
  int *st0_tmp;
  int *st1_tmp;
  int *st2_tmp;
  int **fn_tmp;
  int *reverse_tmp;

  int fCall_Update; 
  double *fEval_chunk;
  int *fIH_chunk;
  int *fNv_chunk;
  int *fNa_chunk;
  int *fSt0_chunk;
  int *fSt1_chunk;
  int *fSt2_chunk;
  int **ffn_chunk;
  int *fReverse_chunk;
  int *fIndex_chunk;
  int *fIndex_new;
  int *fn;

  int fNum_threads;


 
  int inline index1( int i ){

    assert( 0 <= i && i <= 2*fN1 ); // 最終的には削除（計算コスト大）

    if ( i >= fN1 )
      i -= fN1;

    return i;
  }
  
  int inline index2( int i ){

    assert( 0 <= i && i <= 2*fN2 );  // 最終的には削除（計算コスト大）
    
    if ( i >= fN2 )
      i -= fN2;
    
    return i;
  }
  
};


#endif
