/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#define __ENVIRONMENT_I__

#ifndef __SEARCH_I_GET__ 
#include "search_I_get_conf.h"
#endif

class TEnvironment_I {
 public:
  TEnvironment_I(); 
  ~TEnvironment_I(); 
  void SetTarget();          // set the goal figure
  virtual void SetInit();    // initial setting
  virtual void DoIt();       // main procedure
  void Write_Conf_data( int fw );
    
  TSearch_I_get *tSearch0;
  TSearch_I_get *tSearch1;
  
  int fN1, fN2;  // see operation_I.h
  int fN1_in, fN2_in; 
  int fNN1, fNN2;
  VectorXd fW1, fW2; 
  VectorXd fW1_in, fW2_in; 
  VectorXd fWW1, fWW2; 
  MatrixXi fKK1i, fKK2i;

  double fTimeStart, fTimeEnd; // for execution tme
  char *fInstanceName1;          // name of the target figure
  char *fInstanceName2;          // name of the target figure
  char *fOutputName = NULL;     // output file name
  char *fFileNameConfiguration; // file name of the configuration data

  int fIH;
  int fType;
  int fReflection;
  int fReverse;
  int fNum_top;

  int fNum_threads;
};

#endif
