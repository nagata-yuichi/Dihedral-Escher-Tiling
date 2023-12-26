/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env_display.h"
#endif


TEnvironment::TEnvironment()
{
  
}

TEnvironment::~TEnvironment()
{

}


void TEnvironment::ReadTile()
{
  FILE* fp = fopen( fResultFileName, "r" );
  double x, y;
  int n,n_in, na, ih, nv;
  double eval;
  int dumy; 
  int st0, st1, st2;
  int count;
  char word[ 80 ];
  
  dumy = fscanf( fp, "%lf", &fTime );

  dumy = fscanf( fp, "%d%d", &fN1, &fN1_in );
  fNN1 = fN1+fN1_in;

  fW1 = VectorXd(2*fN1);
  fWW1 = VectorXd(2*fNN1);  
  fKK1i = MatrixXi(fNN1,fNN1);
  
  count = 0;
  for( int i = 0; i < fNN1; ++i ){
    dumy = fscanf( fp, "%lf%lf", &x, &y );
    fWW1(count++) = x;
    fWW1(count++) = y;
  }
  fW1 = fWW1.head(2*fN1);

  if( fN1_in != 0 ){
    for( int i = 0; i < fNN1; ++i ){
      for( int j = 0; j < fNN1; ++j ){
	int a;
	dumy = fscanf( fp, "%d", &a );
	fKK1i(i,j) = a;
      }
    }
  }

  dumy = fscanf( fp, "%d%d", &fN2, &fN2_in );
  fNN2 = fN2+fN2_in;

  fW2 = VectorXd(2*fN2);
  fWW2 = VectorXd(2*fNN2);  
  fKK2i = MatrixXi(fNN2,fNN2);
  
  count = 0;
  for( int i = 0; i < fNN2; ++i ){
    dumy = fscanf( fp, "%lf%lf", &x, &y );
    fWW2(count++) = x;
    fWW2(count++) = y;
  }
  fW2 = fWW2.head(2*fN2);

  if( fN2_in != 0 ){
    for( int i = 0; i < fNN2; ++i ){
      for( int j = 0; j < fNN2; ++j ){
	int a;
	dumy = fscanf( fp, "%d", &a );
	fKK2i(i,j) = a;
      }
    }
  }

  int num_top = 1000;

  fEval_top = new double [ num_top ];
  fIndex_top = new int [ num_top ];

  fNv_top = new int [ num_top ];
  fIH_top = new int [ num_top ];
  fN_top = new int [ num_top ];
  fNa_top = new int [ num_top ];
  fSt0_top = new int [ num_top ];
  fSt1_top = new int [ num_top ];
  fSt2_top = new int [ num_top ];

  fU_top = new VectorXd [ num_top ];
  fUa_top = new VectorXd [ num_top ];
  fUU1_top = new VectorXd [ num_top ];
  fUU2_top = new VectorXd [ num_top ];
  ffn_top = new int* [ num_top ];
  for( int i = 0; i < num_top; ++i ){
    fU_top[i] = VectorXd::Zero(2*(fN1+fN2));
    fUa_top[i] = VectorXd::Zero(2*(fN1+fN2));
    fUU1_top[i] = VectorXd::Zero(2*(fN1+fN1_in));
    fUU2_top[i] = VectorXd::Zero(2*(fN2+fN2_in));
    ffn_top[i] = new int [ 6 ];
  }
  fReverse_top = new int [ num_top ];

  fNum_top = 0;
  while( 1 )
  {
    dumy = fscanf( fp, "%d%d%d%lf" ,&ih,&n,&nv,&eval );
    if( feof( fp ) != 0 )
      break;
    dumy = fscanf( fp, "%d" ,&na );

    dumy = fscanf( fp, "%s", word ); 
    fIndex_top[ fNum_top ] = fNum_top;
    fIH_top[ fNum_top ] = ih;
    fN_top[ fNum_top ] = n;
    fNa_top[ fNum_top ] = na;
    fEval_top[ fNum_top ] = eval;
    fNv_top[ fNum_top ] = nv;
    for( int i = 0; i < nv; ++i ){
      int a;
      dumy = fscanf( fp, "%d", &a );
      ffn_top[ fNum_top ][i]  = a;
    }

    dumy = fscanf( fp, "%s", word ); 
    dumy = fscanf( fp, "%d%d%d", &fSt0_top[ fNum_top ], &fSt1_top[ fNum_top ], &fSt2_top[ fNum_top ] );
    dumy = fscanf( fp, "%s", word );
    dumy = fscanf( fp, "%d", &fReverse_top[ fNum_top ] );

    for( int i = 0; i < n; ++i ) {
      dumy = fscanf( fp, "%lf%lf", &x, &y );
      fU_top[ fNum_top ](2*i) = x;
       fU_top[ fNum_top ](2*i+1) = y;
    }

    for( int i = 0; i < na+2; ++i ) { 
      dumy = fscanf( fp, "%lf%lf", &x, &y );
      fUa_top[ fNum_top ](2*i) = x;
      fUa_top[ fNum_top ](2*i+1) = y;
    }

    for( int i = 0; i < fNN1; ++i ) {
      dumy = fscanf( fp, "%lf%lf", &x, &y );
      fUU1_top[ fNum_top ](2*i) = x;
      fUU1_top[ fNum_top ](2*i+1) = y;
    }
    for( int i = 0; i < fNN2; ++i ) {
      dumy = fscanf( fp, "%lf%lf", &x, &y );
      fUU2_top[ fNum_top ](2*i) = x;
      fUU2_top[ fNum_top ](2*i+1) = y;
    }
    ++fNum_top;
  }
  printf( "num_top = %d\n", fNum_top );
  fclose( fp );

  this->SortIndex( fEval_top, fNum_top, fIndex_top, fNum_top );
}


void TEnvironment::SetInit()
{
  tDisplay = new TDisplay();
  tDisplay->SetInit( fN1, fN1_in, fWW1, fKK1i, fN2, fN2_in, fWW2, fKK2i );
  tDisplay->fScale_T = fScale_T;
}


void TEnvironment::DoIt()
{
  this->ReadTile();
  this->SetInit();

  tDisplay->Goal_mesh1();
  tDisplay->Goal_mesh2();

  Display_top();
  // tDisplay->All_tiles( fUU1_top, fUU2_top );

  while(1){
    char ch;
    fprintf(stderr, "Hit Return key to continue.\n");
    ch = getc(stdin);
    if( ch == 'q' )
      break;
  }
}


void TEnvironment::Display_top()
{
  int order = 0;
  int index; 
  char ch;
  int mode;
  int fIH;

  mode = 1;
  while( 1 ){
    index = fIndex_top[order];
    printf( "order = %2d eval=%lf: IH=%2d", order, fEval_top[index], fIH_top[index] );
    printf(": a = %d", fNa_top[index] ); 
    printf(" (");
    for( int i = 0; i < fNv_top[index]; ++i )
      printf( "%d ", ffn_top[index][i] );
    printf(")");
    printf(": st = %d %d %d : r = %d", fSt0_top[index], fSt1_top[index], fSt2_top[index], fReverse_top[index] ); 
    fIH = fIH_top[index];

    VectorXd uu1_t = fUU1_top[index];
    VectorXd uu2_t = fUU2_top[index];
    if( fReverse_top[index] == 1 )
      for( int i = 0; i < fNN2; ++i )
	uu2_t(2*i) *= -1.0;

    tDisplay->Tile_mesh1( uu1_t );
    tDisplay->Tile_mesh2( uu2_t );

    int n, na;
    n = fN_top[index];
    na = fNa_top[index];

    VectorXd u = fU_top[index].head(2*n);
    VectorXd ua = fUa_top[index].head(2*na+4);
    int st0 = fSt0_top[index];

    tDisplay->Dihedral( u, ua, st0, fNv_top[index], ffn_top[index], fUU1_top[index], fUU2_top[index]  );
    tDisplay->Dihedral2( u, ua, st0, fUU1_top[index], fUU2_top[index] );

	
    if( fIH == 4 )
      tDisplay->Tiling_IH4( u, ua, ffn_top[index], fNv_top[index], fUU1_top[index], fUU2_top[index] );
    if( fIH == 5 )
      tDisplay->Tiling_IH5( u, ua, ffn_top[index], fNv_top[index], fUU1_top[index], fUU2_top[index] );
    if( fIH == 6 )
      tDisplay->Tiling_IH6( u, ua, ffn_top[index], fNv_top[index], fUU1_top[index], fUU2_top[index] ); 
    if( fIH == 1 )
      tDisplay->Tiling_IH1( u, ua, ffn_top[index], fNv_top[index] );
    if( fIH == 2 )
      tDisplay->Tiling_IH2( u, ua, ffn_top[index], fNv_top[index] );
    if( fIH == 3 )
      tDisplay->Tiling_IH3( u, ua, ffn_top[index], fNv_top[index] );


    ch = getc(stdin);
    if( ch == 'f' )
      mode = 1;
    if( ch == 'b' )
      mode = -1;

    order += mode;

    if( order == fNum_top ) order = 0;
    if( order == -1 ) order = fNum_top-1;
  }

}



// Arg[]: input array
// numOfArg: the number of elements of Arg[]
// indexOrdered[i] <- the position of the i-th smallest element of Arg[]
// only the top numOfOrd elements are considered
void TEnvironment::SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest;
  double valueBest;
  int checked[numOfArg];

  assert( Arg[0] < 99999999999.9 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = 99999999999.9;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] < valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }
}


