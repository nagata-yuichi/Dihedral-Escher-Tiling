/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env.h"
#endif

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
  assert(  argc >= 4 );
  TEnvironment* tEnv = NULL;
  tEnv = new TEnvironment();

  tEnv->fInstanceName = argv[1]; // file name of the target shape 
  sscanf( argv[2], "%lf", &(tEnv->fAlpha) ); 
  sscanf( argv[3], "%lf", &(tEnv->fBeta) ); 
  tEnv->fOutputName = argv[4];   // output file name 
  
  tEnv->DoIt();
  
  return 0;
}
