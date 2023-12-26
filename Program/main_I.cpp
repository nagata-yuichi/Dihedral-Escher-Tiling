/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env_I.h"
#endif

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
  TEnvironment_I* tEnv = NULL;
  tEnv = new TEnvironment_I();

  tEnv->fInstanceName1 = argv[1]; // file name of the target shape
  tEnv->fInstanceName2 = argv[2]; // file name of the target shape 
  tEnv->fOutputName = argv[3];   // output file name 
  sscanf( argv[4], "%d", &(tEnv->fNum_top) ); 
  sscanf( argv[5], "%lf", &(tEnv->fAlpha) );
  tEnv->fFileNameConfiguration = argv[6]; // file name of the configuration data
  sscanf( argv[7], "%d", &(tEnv->fNum_threads) );
  tEnv->DoIt();
  
  return 0;
}
