/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_CONF__
#define __SEARCH_I_CONF__

#define IR   // Set I or IR
// I: E_I distance, IR: E_IR distance

#ifdef I
  #ifndef  __OPERATOR_I__
  #include "operator_I.h"
  #define CLASS_NAME TOperator_I
  #endif
#endif
#ifdef IR
  #ifndef  __OPERATOR_IR__
  #include "operator_IR.h"
  #define CLASS_NAME TOperator_IR
  #endif
#endif


class TSearch_I_conf
{
 public:
  TSearch_I_conf( int N1, int N1_in, int N2, int N2_in );
  ~TSearch_I_conf(); 
  
  void SetParameter();  // set parameter values
  void Define();        // declare arryas
  void SetInit( VectorXd w1, VectorXd w1_in, MatrixXi kk1i, VectorXd w2, VectorXd w2_in, MatrixXi kk2i );

  void DoIt();          // main process of the EST
  void IH4();
  void IH5();
  void IH6();
  void IH1();
  void IH2();
  void IH3();

  void Update_UU_top( CLASS_NAME* op );
  void Set_candi( int num_candi, double* eval_candi, int* IH_candi, int* Na_candi, int** fn_candi, int* st0_candi, int* st1_candi, int* st2_candi, int* reverse_candi );


  int fN1, fN1_in, fNN1, fN2, fN2_in, fNN2;
  double fAlpha;
  CLASS_NAME **tOp;

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

  // Conf
  int fNum_candi;     // the number of template configurations read from the file
  double *fEval_candi;     // the number of template configurations read from the file
  int *fIH_candi;     // [s] -> IH type of the s-th configuration
  int *fNa_candi;     // [s] -> Na
  int **fki_candi;    // [s] -> (k_1, k_2, ...) of the s-th configuration
  int *fSt0_candi;    // [s] -> st1 (=j) of the the s-th configuration
  int *fSt1_candi;    // [s] -> st1 (=j) of the the s-th configuration
  int *fSt2_candi;    // [s] -> st1 (=j) of the the s-th configuration
  int *fReverse_candi; //  
  int fCurrent_candi; // the number of the template configuration that is currently being optimized

  int fCandi_start_IH1, fCandi_end_IH1;
  int fCandi_start_IH2, fCandi_end_IH2;
  int fCandi_start_IH3, fCandi_end_IH3;
  int fCandi_start_IH4, fCandi_end_IH4;
  int fCandi_start_IH5, fCandi_end_IH5;
  int fCandi_start_IH6, fCandi_end_IH6;

  int fNum_threads;
};


#endif
