/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env.h"
#endif

TEnvironment::TEnvironment()
{
  fOutputName = NULL;
}

TEnvironment::~TEnvironment()
{

}

void TEnvironment::SetTarget()
{
  FILE* fp = fopen( fInstanceName, "r" );
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }

  int dumy = fscanf( fp, "%d", &fN );
  double xd, yd;

  fW = VectorXd(2*fN);
 
  int count = 0;
  while( 1 )
  {
    dumy = fscanf( fp, "%lf%lf", &xd, &yd );
    if( feof( fp ) != 0 )
      break;
    fW(count++) = xd;
    fW(count++) = yd;
  }
  assert( count == 2*fN );
}

void TEnvironment::SetInit()
{
  this->SetTarget();
  tSearch = new TSearch( fN );
  tSearch->tDisplay = tDisplay;

  tSearch->fAlpha = fAlpha; 
  tSearch->fBeta = fBeta; 
  tSearch->SetInit( fW );
}


void TEnvironment::DoIt()
{

  char ch;
  
  this->SetInit();

  tDisplay.Disp_goal( fW, tSearch->fG );


  tSearch->DoIt();

  this->WriteData();

}


void TEnvironment::WriteData()
{
  string filename;

  if( fOutputName != NULL ){
    filename = fOutputName;
  }
  else {
    filename = fInstanceName;
    filename.pop_back();
    filename.pop_back();
    filename.pop_back();
    filename.pop_back();
  }

  filename += "_";
  filename += to_string(tSearch->fN_in); 
  filename += ".dat";

  FILE *fp;

  fp = fopen( filename.c_str(), "w" );
  fprintf( fp, "%d %d\n", tSearch->fN, tSearch->fN_in );
  for( int i = 0; i < tSearch->fNN; ++i ){
    fprintf( fp, "%.1f %.1f \n", tSearch->fWW(2*i), tSearch->fWW(2*i+1) ); 
  }

  for( int i = 0; i < tSearch->fNN; ++i ){
    for( int j = 0; j < tSearch->fNN; ++j ){
      int a = tSearch->fKKi(i,j);
      if( i == j )
	a = 0;
      if( a != 0 )
	a = 1;
      fprintf( fp, "%d ", a );
    }
    fprintf( fp, "\n" );
  }
  

  fclose( fp );


}
