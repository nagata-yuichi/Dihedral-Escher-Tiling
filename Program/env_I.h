/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#define __ENVIRONMENT_I__

#ifndef __SEARCH_I_CONF__ 
#include "search_I_conf.h"
#endif

class TEnvironment_I {
 public:
  TEnvironment_I(); 
  ~TEnvironment_I(); 
  void SetTarget();          // set the goal figure
  virtual void SetInit();    // initial setting
  virtual void DoIt();       // main procedure
  void Display_top();        // display top tile shapes 
  void Write_Tile();         // write top tilig shapes to a file
  void ReadConfigurations(); // read configuration data 

  TSearch_I_conf *tSearch;
  int fN1, fN2; 
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

  int fNum_top;
  double fAlpha;      // parameter alpha for the E_I distance
  int fNum_threads;
};

#endif
