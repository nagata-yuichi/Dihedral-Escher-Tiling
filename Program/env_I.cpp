/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#include "env_I.h"     
#endif

TEnvironment_I::TEnvironment_I()
{
}


TEnvironment_I::~TEnvironment_I()
{
}


void TEnvironment_I::SetInit()
{
  this->SetTarget();
  tSearch = new TSearch_I_conf( fN1, fN1_in, fN2, fN2_in );
  tSearch->SetParameter();
  tSearch->fNum_top = fNum_top; 
  tSearch->fAlpha = fAlpha;
  tSearch->fNum_threads = fNum_threads;

  tSearch->Define();
  tSearch->SetInit( fW1, fW1_in, fKK1i, fW2, fW2_in, fKK2i );
}


void TEnvironment_I::DoIt()
{
  this->SetInit();

  this->ReadConfigurations();  

  fTimeStart = omp_get_wtime();
  tSearch->DoIt();
  fTimeEnd = omp_get_wtime();

  printf( "time = %.2f \n", (fTimeEnd-fTimeStart) );

  this->Write_Tile();  
}


void TEnvironment_I::SetTarget()
{
  double xd, yd;
  int dumy, d;
  FILE* fp;


  //// goal 1 ////
  fp = fopen( fInstanceName1, "r" );
  if( fp == NULL ){
    printf( "file1 not found\n");
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
  fp = fopen( fInstanceName2, "r" );
  if( fp == NULL ){
    printf( "file2 not found\n");
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


void TEnvironment_I::Write_Tile()
{
  FILE *fp;
  char filename[80];
  int index;

  if( tSearch->fNum_top > 1000 )
    return;
  
  for( int order = 0; order < tSearch->fNum_top; ++order ){
    index = tSearch->fIndex_top[order];
    printf( "order = %2d eval=%lf: IH=%2d", order, tSearch->fEval_top[index], tSearch->fIH_top[index] );
    printf(": a = %d", tSearch->fNa_top[index] ); 
    printf(" (");
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      printf( "%d ", tSearch->ffn_top[index][i] );
    printf(")");
    printf(": st = %d %d %d : r = %d", tSearch->fSt0_top[index], tSearch->fSt1_top[index], tSearch->fSt2_top[index],  tSearch->fReverse_top[index] );
    printf("\n");
  }

  if( fOutputName == NULL )
    return;

  string name1, name2;
  name1 = fInstanceName1;
  name2 = fInstanceName2;
  name1.pop_back();name1.pop_back();name1.pop_back();name1.pop_back();
  name2.pop_back();name2.pop_back();name2.pop_back();name2.pop_back();
  sprintf( filename, "%s_%s_%s.tile", fOutputName, name1.c_str(), name2.c_str());

  //  sprintf( filename, "%s_%s_%s.tile", fOutputName, fInstanceName1, fInstanceName2 );

  fp = fopen( filename, "w" );

  // time
  fprintf( fp, "%.2f \n", (fTimeEnd-fTimeStart) );

  // instance
  fprintf( fp, "%d %d\n", fN1, fN1_in ); 
  for( int i = 0; i < fNN1; ++i )
    fprintf( fp, "%lf %lf\n", fWW1(2*i), fWW1(2*i+1) ); 

  // KKi matrix   
  for( int i = 0; i < fNN1; ++i ){
    for( int j = 0; j < fNN1; ++j ){
      fprintf( fp, "%d ", fKK1i(i,j) );
    }
    fprintf( fp, "\n" );
  }

  // instance
  fprintf( fp, "%d %d\n", fN2, fN2_in ); 
  for( int i = 0; i < fNN2; ++i )
    fprintf( fp, "%lf %lf\n", fWW2(2*i), fWW2(2*i+1) ); 

  // KKi matrix   
  for( int i = 0; i < fNN2; ++i ){
    for( int j = 0; j < fNN2; ++j ){
      fprintf( fp, "%d ", fKK2i(i,j) );
    }
    fprintf( fp, "\n" );
  }
  

  // tile shapes
  for( int s = 0; s < tSearch->fNum_top; ++s ){
    index = tSearch->fIndex_top[s];
    fprintf( fp, "%d %d %d %lf\n", tSearch->fIH_top[index], tSearch->fN_top[index], tSearch->fNv_top[index], tSearch->fEval_top[index] );
    fprintf( fp, "%d : ", tSearch->fNa_top[index] );
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      fprintf( fp, "%d ", tSearch->ffn_top[index][i] );
    fprintf( fp, ": " );
    fprintf( fp, "%d %d %d : %d", tSearch->fSt0_top[index], tSearch->fSt1_top[index], tSearch->fSt2_top[index], tSearch->fReverse_top[index] );
    fprintf( fp, "\n" );

    for( int i = 0; i < tSearch->fN_top[index]; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fU_top[index](2*i), (double)tSearch->fU_top[index](2*i+1) );
    }
    for( int i = 0; i < tSearch->fNa_top[index]+2; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fUa_top[index](2*i), (double)tSearch->fUa_top[index](2*i+1) );
    }
    
    for( int i = 0; i < fN1+fN1_in; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fUU1_top[index](2*i), (double)tSearch->fUU1_top[index](2*i+1) );
    }
    for( int i = 0; i < fN2+fN2_in; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fUU2_top[index](2*i), (double)tSearch->fUU2_top[index](2*i+1) );
    }
  } 
  fclose(fp);
}


void TEnvironment_I::ReadConfigurations() 
{
  int ih, nv, ni;
  double eval;
  double dd;
  int na;
  int st0, st1, st2;
  int reverse;
  FILE *fp;
  int ffi[7];
  double *eval_candi;
  int *IH_candi;
  int *Na_candi;
  int *Nv_candi;
  int **fn_candi;
  int *st0_candi;
  int *st1_candi;
  int *st2_candi;
  int *reverse_candi;
  int num_candi;
  double time;
  int dumy;
  char word[ 80 ];

  fp = fopen( fFileNameConfiguration, "r" );
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }

  dumy = fscanf( fp, "%lf", &time );
  dumy = fscanf( fp, "%d", &num_candi );


  eval_candi = new double [ num_candi ];
  IH_candi = new int [ num_candi ];
  Na_candi = new int [ num_candi ];
  Nv_candi = new int [ num_candi ];
  st0_candi = new int [ num_candi ];
  st1_candi = new int [ num_candi ];
  st2_candi = new int [ num_candi ];
  fn_candi = new int* [ num_candi ];
  for( int i = 0; i < num_candi; ++i )
    fn_candi[i] = new int [ 6 ];
  reverse_candi = new int [ num_candi ];

  int count = 0;
  while( 1 )
  {
    dumy = fscanf( fp, "%lf%d%d", &eval, &ih, &nv );
    if( feof( fp ) != 0 )
      break;
    dumy = fscanf( fp, "%d", &na );
    eval_candi[count] = eval;
    IH_candi[count] = ih;
    Nv_candi[count] = nv;
    Na_candi[count] = na;

    dumy = fscanf( fp, "%s", word ); 
    for( int i = 0; i < nv; ++i ) {
      dumy = fscanf( fp, "%d", &ni );
      fn_candi[ count ][i] = ni;
    }
    dumy = fscanf( fp, "%s", word ); 
    dumy = fscanf( fp, "%d%d%d", &st0, &st1, &st2 );
    st0_candi[ count ] = st0;
    st1_candi[ count ] = st1;
    st2_candi[ count ] = st2;
    dumy = fscanf( fp, "%s", word );
    dumy = fscanf( fp, "%d", &reverse );
    reverse_candi[ count ] = reverse;

    ++count;
  }
  printf( "num_candi = %d\n", count );
  fclose( fp );

  if( count != num_candi )
    assert( 1 == 2 );

  tSearch->Set_candi( num_candi, eval_candi, IH_candi, Na_candi, fn_candi, st0_candi, st1_candi, st2_candi, reverse_candi );
}

