/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#include "env_I_get_conf.h"     
#endif

TEnvironment_I::TEnvironment_I()
{
}


TEnvironment_I::~TEnvironment_I()
{
}


void TEnvironment_I::SetInit()
{
  fReflection = 0;
  this->SetTarget(); // fReflectionで読込ファイルを場合分け
  tSearch0 = new TSearch_I_get( fN1, fN2 );
  tSearch0->SetParameter();
  tSearch0->fNum_top = fNum_top;  
  tSearch0->Define();
  tSearch0->SetInit( fW1, fW2 );

  fReflection = 1;
  this->SetTarget(); // fReflectionで読込ファイルを場合分け
  tSearch1 = new TSearch_I_get( fN1, fN2 );
  tSearch1->SetParameter();
  tSearch1->fNum_top = fNum_top;  
  tSearch1->Define();
  tSearch1->SetInit( fW1, fW2 );
}



// All
void TEnvironment_I::DoIt()
{
  int num_threads_here; 

  fTimeStart = omp_get_wtime();
  
  this->SetInit();

  omp_set_nested(1);

  if( fNum_threads >= 2 ){
    num_threads_here = 2;
    tSearch0->fNum_threads = fNum_threads/2;
    tSearch1->fNum_threads = fNum_threads/2;    
  }
  else{
    num_threads_here = 1;
    tSearch0->fNum_threads = fNum_threads;
    tSearch1->fNum_threads = fNum_threads;    
  }
  

#pragma omp parallel sections  num_threads(num_threads_here)
  {

#pragma omp section
    {
      tSearch0->fReflection = 0;
        tSearch0->DoIt_IH456();
      // tSearch0->DoIt_IH123();
    }

#pragma omp section
    {
    tSearch1->fReflection = 1;
      tSearch1->DoIt_IH456();
      // tSearch1->DoIt_IH123();
    }
  }
  fTimeEnd = omp_get_wtime();
  
  this->Write_Conf_data( 0 );   
  printf( "time = %.2f \n", (fTimeEnd-fTimeStart) );
}



void TEnvironment_I::SetTarget()
{
  double xd, yd;
  int dumy, d;
  FILE* fp;


  //// goal 1 ////
  if( fReflection == 0 ) // ここ注意
    fp = fopen( fInstanceName1, "r" );
  else if ( fReflection == 1 )
    fp = fopen( fInstanceName2, "r" );
  else
    assert( 1 == 2 );
  
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }
  dumy = fscanf( fp, "%d%d", &fN1, &fN1_in );

  fNN1 = fN1+fN1_in;
  fW1 = VectorXd(2*fN1);
  fW1_in = VectorXd(2*fN1_in);
  fWW1 = VectorXd(2*fNN1);
  fKK1i = MatrixXi(fNN1,fNN1);
 
  for( int i = 0; i < fNN1; ++i ){
    dumy = fscanf( fp, "%lf%lf", &xd, &yd );
    fWW1(2*i) = xd;
    fWW1(2*i+1) = yd;
  }

  fW1 = fWW1.head(2*fN1);
  fW1_in = fWW1.tail(2*fN1_in);

  for( int i = 0; i < fNN1; ++i ){
    for( int j = 0; j < fNN1; ++j ){
      assert( feof( fp ) == 0 );
      dumy = fscanf( fp, "%d", &d );
      fKK1i(i,j) = d;
    }
  }
  dumy = fscanf( fp, "%d", &d );
  assert( feof( fp ) != 0 );

  fclose( fp );

  //// goal 2 ////
  if( fReflection == 0 ) // ここ注意
    fp = fopen( fInstanceName2, "r" );
  else if( fReflection == 1 )
    fp = fopen( fInstanceName1, "r" );
  else
    assert( 1 == 2 );
  
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }
  dumy = fscanf( fp, "%d%d", &fN2, &fN2_in );

  fNN2 = fN2+fN2_in;
  fW2 = VectorXd(2*fN2);
  fW2_in = VectorXd(2*fN2_in);
  fWW2 = VectorXd(2*fNN2);
  fKK2i = MatrixXi(fNN2,fNN2);
 
  for( int i = 0; i < fNN2; ++i ){
    dumy = fscanf( fp, "%lf%lf", &xd, &yd );
    fWW2(2*i) = xd;
    fWW2(2*i+1) = yd;
  }

  fW2 = fWW2.head(2*fN2);
  fW2_in = fWW2.tail(2*fN2_in);

  for( int i = 0; i < fNN2; ++i ){
    for( int j = 0; j < fNN2; ++j ){
      assert( feof( fp ) == 0 );
      dumy = fscanf( fp, "%d", &d );
      fKK2i(i,j) = d;
    }
  }
  dumy = fscanf( fp, "%d", &d );
  assert( feof( fp ) != 0 );

  fclose( fp );
}


// Get Conf
void TEnvironment_I::Write_Conf_data( int fw )
{
  FILE *fp;
  char filename[80];
  int index;
  static int flag_init = 0;

  if( fOutputName == NULL )
    return;

  string name1, name2;
  name1 = fInstanceName1;
  name2 = fInstanceName2;
  name1.pop_back();name1.pop_back();name1.pop_back();name1.pop_back();
  name2.pop_back();name2.pop_back();name2.pop_back();name2.pop_back();
  sprintf( filename, "%s_%s_%s.conf", fOutputName, name1.c_str(), name2.c_str());

  if( fw == 0 ) 
    fp = fopen( filename, "w" );
  else{ 
    if( flag_init == 0 ){
      fp = fopen( filename, "w" );
      fclose(fp);
      flag_init = 1;
    }
    fp = fopen( filename, "a" );
  }

  // time
  fprintf( fp, "%.2f \n", (fTimeEnd-fTimeStart) );
  fprintf( fp, "%d \n", tSearch0->fNum_top );

  // conf data


  double value, value0, value1;
  int num0, num1;
  int num_candi;
  int index0, index1;

  num0 = 0;
  num1 = 0;
  num_candi = 0;
  value = -1.0;

  while(1){
    index0 = tSearch0->fIndex_top[num0]; 
    index1 = tSearch1->fIndex_top[num1];

    value0 = tSearch0->fEval_top[index0];
    value1 = tSearch1->fEval_top[index1];
    
    if( value0 < value1 ){
      if ( fabs( value0 - value ) > 0.0000001 ){
	fprintf( fp, "%.9f %d %d\n", tSearch0->fEval_top[index0], tSearch0->fIH_top[index0], tSearch0->fNv_top[index0] );
	fprintf( fp, "%d : ", tSearch0->fNa_top[index0] );
	for( int i = 0; i < tSearch0->fNv_top[index0]; ++i )
	  fprintf( fp, "%d ", tSearch0->ffn_top[index0][i] );
	fprintf( fp, ": %d %d %d : %d", tSearch0->fSt0_top[index0], tSearch0->fSt1_top[index0], tSearch0->fSt2_top[index0], tSearch0->fReverse_top[index0] );
	fprintf( fp, "\n" );

	value = value0;
	++num_candi; 
      }

      ++num0;
    } 
    else{
      if ( fabs( value1 - value ) > 0.0000001 ){
	fprintf( fp, "%.9f %d %d\n", tSearch1->fEval_top[index1], tSearch1->fIH_top[index1], tSearch1->fNv_top[index1] );
	fprintf( fp, "%d : ", tSearch1->fNa_top[index1] );
	for( int i = 0; i < tSearch1->fNv_top[index1]; ++i )
	  fprintf( fp, "%d ", tSearch1->ffn_top[index1][i] );
	fprintf( fp, ": %d %d %d : %d", tSearch1->fSt0_top[index1], tSearch1->fSt1_top[index1], tSearch1->fSt2_top[index1], tSearch1->fReverse_top[index1] );
	fprintf( fp, "\n" );

	value = value1;
	++num_candi; 
      }
      
      ++num1;
    }

    if( num_candi == tSearch0->fNum_top )
      break;
    if( num0 == tSearch0->fNum_top || num1 == tSearch1->fNum_top )
      assert( 1 == 2 );
  }

  fclose(fp);

}
