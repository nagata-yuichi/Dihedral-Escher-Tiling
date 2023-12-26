/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH__
#include "search.h"
#endif

TSearch::TSearch( int N ) 
{
  fN = N;
  fNN = fN;
  fW = VectorXd::Zero(2*fN);
  fG = MatrixXd::Zero(fN,fN);
  fW_in = VectorXd::Zero(2*fN*9); 
  fWW = VectorXd::Zero(2*fN*10); 
  fKKi = MatrixXi::Zero(fNN,fNN); 
  fNei = MatrixXi::Zero(2*fN*10,20);
  fNei_Size = VectorXi::Zero(2*fN*10);
}

TSearch::~TSearch()
{

}

void TSearch::SetParameter() // call from env.cpp
{
  
}

void TSearch::SetInit( VectorXd w ) // call from env.cpp 
{
  fW = w;
  fG = MatrixXd::Identity(fN,fN);
  fN_in = 0;
  fNN = fN;
  fWW.head(2*fN) = fW;
  fKKi = MatrixXi::Zero(fNN,fNN); 
}

void TSearch::DoIt() 
{
  this->SetRegion(); 
  this->SetInitialInnerPoints();
  this->SetK();
  tDisplay.Disp_AdjacentRelation( fN, fWW, fKKi );
  this->ModiInnerPoints();
  tDisplay.Disp_AdjacentRelation2( fN, fWW, fKKi );
}



void TSearch::SetRegion() 
{
  double x_min, x_max, y_min, y_max;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < fN; ++i ){
    if( fW(ix(i)) < x_min )
      x_min = fW(ix(i));
    if( fW(ix(i)) > x_max )
      x_max = fW(ix(i));
    if( fW(iy(i)) < y_min )
      y_min = fW(iy(i));
    if( fW(iy(i)) > y_max )
      y_max = fW(iy(i));
  }
  
  fLengthX = x_max-x_min+1;
  fLengthY = y_max-y_min+1;
  fX_min = x_min;
  fY_min = y_min;

  for( int i = 0; i < fN; ++i ){ // 端を(0,0)にする
    fW(ix(i)) -= fX_min;
    fW(iy(i)) -= fY_min; 
  } 
  fWW.head(2*fN) = fW;

  fInnerRegion = MatrixXi::Zero((int)fLengthX,(int)fLengthY);

  for( int px = 0; px < (int)fLengthX; ++px ){
    for( int py = 0; py < (int)fLengthY; ++py ){
      if( this->Check_InnerResion( (double)px, (double)py ) )
	fInnerRegion(px,py) = 1;
      else {
	fInnerRegion(px,py) = 0;
      }
    }
  }

  // L
  fL = 0.0;
  for( int i = 0; i < fN; ++i )
    fL += sqrt(pow(fW(ix(i+1))-fW(ix(i)),2.0) + pow(fW(iy(i+1))-fW(iy(i)),2.0));
  fL /= (double)fN;
  fL_min = fL*fAlpha - fL*fBeta;
  fL_max = fL*fAlpha + fL*fBeta;
}

void TSearch::SetInitialInnerPoints() 
{
  double px, py;
  double value_best, value, d;
  int px_best, py_best;
  VectorXd ww(2*fN*10);

  fN_in = 0;
  fNN = fN;
  ww.head(2*fN) = fW;

  while(1){

    value_best = 0.0;
    for( int i = 0; i < (int)fLengthX; ++i ){
      for( int j = 0; j < (int)fLengthY; ++j ){
	px = (double)i; // 図形を囲うグリッドを網羅的にチェック
	py = (double)j;

	if( fInnerRegion(i,j) == 1 ){ // 図形内部の点 // outer (comment out) 
	  value = 0.0;
	  for( int k = 0; k < fNN; ++k ){
	    d = sqrt( pow(px-ww(iix(k)),2.0) + pow(py-ww(iiy(k)),2.0) );
	    // printf( "%d %d: %lf %lf \n", i,j, d, fL  );
	    // 現在割当済の点との距離が[fL_min,fL_max]の範囲となる点の個数が最大となる位置(px,py)を探索
	    // 同数の場合は距離がfLに近い場所を優先
	    if( d <= fL_min ){
	      value -= 9999999.9;
	    }
	    else if( d <= fL_max ) 
	      value += 1000.0 - fabs(d-fL);
	  }
	  value += 0.01 * drand48(); // tie break
	  if( value_best < value ){
	    value_best = value;
	    px_best = px;
	    py_best = py;
	  }
	} // outer (comment out) 
      }
    }

    if( value_best > 1.0 ){
      ww(2*fNN) = (double)px_best;
      ww(2*fNN+1) = (double)py_best;
      ++fNN;
      printf( "%d (%d,%d) %lf \n", fNN, px_best, py_best, value_best);
    }
    else 
      break;

  }

  fWW.resize(2*fNN);
  fWW = ww.head(2*fNN);
  fN_in = fNN-fN;
}

/*
void TSearch::ModiInnerPoints()
{
  int v1, v2;
  VectorXd diff = VectorXd(2);
  int iter;
  
  iter = 0;
  while(1){
    this->SetK();

    if( rand() % 2 == 0 ){
      v1 = fV_long(0);
      v2 = fV_long(1);

      diff(0) = fWW(iix(v2)) - fWW(iix(v1));
      diff(1) = fWW(iiy(v2)) - fWW(iiy(v1));
      if( v1 >= fN && v2 >= fN ){
	fWW(iix(v1)) += 0.005*diff(0);
	fWW(iiy(v1)) += 0.005*diff(1);
	fWW(iix(v2)) -= 0.005*diff(0);
	fWW(iiy(v2)) -= 0.005*diff(1);
      }
      else if( v1 >= fN && v2 < fN ){
	fWW(iix(v1)) += 0.01*diff(0);
	fWW(iiy(v1)) += 0.01*diff(1);
      }
      else if( v1 < fN && v2 >= fN ){
	fWW(iix(v2)) -= 0.01*diff(0);
	fWW(iiy(v2)) -= 0.01*diff(1);
      }
      else 
	assert( 1 == 2 );
    }
    else {
      v1 = fV_short(0);
      v2 = fV_short(1);

      diff(0) = fWW(iix(v2)) - fWW(iix(v1));
      diff(1) = fWW(iiy(v2)) - fWW(iiy(v1));
      if( v1 >= fN && v2 >= fN ){
	fWW(iix(v1)) -= 0.005*diff(0);
	fWW(iiy(v1)) -= 0.005*diff(1);
	fWW(iix(v2)) += 0.005*diff(0);
	fWW(iiy(v2)) += 0.005*diff(1);
      }
      else if( v1 >= fN && v2 < fN ){
	fWW(iix(v1)) -= 0.01*diff(0);
	fWW(iiy(v1)) -= 0.01*diff(1);
      }
      else if( v1 < fN && v2 >= fN ){
	fWW(iix(v2)) += 0.01*diff(0);
	fWW(iiy(v2)) += 0.01*diff(1);
      }
      else 
	assert( 1 == 2 );
      ++iter;
      printf( "iter = %d\n", iter );
      //      if( iter == 1000 )
      //	break;
    }

    tDisplay.Disp_AdjacentRelation2( fN, fWW, fKKi );
  }
}
*/

void TSearch::ModiInnerPoints()
{
  int v1, v2, v2_max, v2_min;
  VectorXd diff = VectorXd(2);
  VectorXd wa = VectorXd(2);
  int iter;
  double d, d_max, d_min;
  VectorXd tabu= VectorXd::Zero(fNN);

  
  initscr();
  timeout(0);

  printf("push 'q' to stop\n" );
  iter = 0;
  while(1){
    this->SetK();

    while (1){ 
      v1 = rand() % fN_in + fN; // 追加点をダンダムに選択
      if( tabu(v1) <= iter ){
	tabu(v1) = iter + (int)(fN_in*0.7); // tabu paremeter
	break;
      }
    }

    d_max = 0.0;
    d_min = 99999999.0;
    for( int k = 0; k < fNei_Size(v1); ++k ){
      v2 = fNei(v1,k);
      d = sqrt( pow(fWW(iix(v2))-fWW(iix(v1)),2.0) + pow(fWW(iiy(v2))-fWW(iiy(v1)),2.0) );
      if( d > d_max ){
	d_max = d;
	v2_max = v2;
      }
      if( d < d_min ){
	d_min = d;
	v2_min = v2;
      }
    }

    if( rand() % 2 == 0 ){
      v2 = v2_max;

      diff(0) = fWW(iix(v2)) - fWW(iix(v1));
      diff(1) = fWW(iiy(v2)) - fWW(iiy(v1));

      wa(0) = fWW(iix(v1)) + 0.01*diff(0);
      wa(1) = fWW(iiy(v1)) + 0.01*diff(1);
      if( this->Check_InnerResion( wa(0), wa(1) ) ){ 
	fWW(iix(v1)) = wa(0);
	fWW(iiy(v1)) = wa(1);
      }
    }
    else {
      v2 = v2_min;

      diff(0) = fWW(iix(v2)) - fWW(iix(v1));
      diff(1) = fWW(iiy(v2)) - fWW(iiy(v1));

      wa(0) = fWW(iix(v1)) - 0.01*diff(0);
      wa(1) = fWW(iiy(v1)) - 0.01*diff(1);
      if( this->Check_InnerResion( wa(0), wa(1) ) ){ 
	fWW(iix(v1)) = wa(0);
	fWW(iiy(v1)) = wa(1);
      }
    }

    tDisplay.Disp_AdjacentRelation2( fN, fWW, fKKi );

    ++iter;
    // printf( "iter = %d\n", iter );

    int ch = getch();
    if (ch == 'q') 
      break;

  }
  endwin();
  printf( "iter=%d: fN=%d fN_in=%d \n", iter, fN, fN_in );
}


bool TSearch::Check_InnerResion( double px, double py )
{
  int count;
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;
  
  x1 = px;
  y1 = py;
  x2 = x1+10000.0;
  y2 = y1+10.0;

  count = 0;
  for( int i = 0; i < fN; ++i ){
    x3 = fW(ix(i));
    y3 = fW(iy(i));
    x4 = fW(ix(i+1));
    y4 = fW(iy(i+1));

    tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
    td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
    value1 = tc*td;
    
    tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
    td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
    value2 = tc*td;

    if( value1 < 0 && value2 < 0 )
      ++count;
  }  

  if( count % 2 == 1 )  
    return true;
  else 
    return false;

}



void TSearch::SetK() 
{
  int num1, num2, w, dumy;
  
  // General: T1 T2
  double edgeL[fNN*100]; // selected edge length
  int edgeV[fNN*100][2]; // endpoints of the selected edges
  int sortedIndex[fNN*100]; // sorted indices
  int numEdgeV; // the number of the selected edges
  int selectedIndex[fNN*100];
  int numSelected;

  
  /////// edgeL and edgeV //////
  double interval_ave = 0.0;
  double d;
  for( int i = 0; i < fN-1; ++i ){
    d = pow(fW(ix(i+1))-fW(ix(i)),2.0) + pow(fW(iy(i+1))-fW(iy(i)),2.0);
    interval_ave += sqrt(d);
  }
  d = pow(fW(ix(0))-fW(ix(fN-1)),2.0) + pow(fW(iy(0))-fW(iy(fN-1)),2.0);
  interval_ave += sqrt(d);
  interval_ave /= (double)fN;
 
  numEdgeV = 0; 
  for( int i = 0; i < fN; ++i ){ // 輪郭部分
    edgeL[numEdgeV] = 0.0001;
    edgeV[numEdgeV][0] = i;
    edgeV[numEdgeV++][1] = (i+1) % fN;
  }

  for( int i = 0; i < fNN; ++i ){
    for( int j = i+1; j < fNN; ++j ){
      d = pow(fWW(2*i)-fWW(2*j),2.0) + pow(fWW(2*i+1)-fWW(2*j+1),2.0);
      d = sqrt(d);

      int flag = 0;
      if( i != j && d < 2.0*interval_ave &&
      	  ( i >= fN || j >= fN || (fabs(i-j) != 1 && fabs(i-j) != fN-1) )  // 輪郭は除外
      	  && this->Check_Internal_link( i, j) )
	flag = 1;
      //if( i < fN && i < j && this->Check_Internal_link(i, j) == 0 && d < 1.25*interval_ave ) // 外部linkの追加
      //flag = 1;

      if(flag == 1 ){
	edgeL[numEdgeV] = d;
	edgeV[numEdgeV][0] = i;
	edgeV[numEdgeV++][1] = j;
      }
    }
  }

  this->SortIndex( edgeL, numEdgeV, sortedIndex, numEdgeV );

  int flag_intersection;
  int a1, a2, b1, b2;
  numSelected = 0; 
  for( int s1 = 0; s1 < numEdgeV; ++s1 ){
    a1 = edgeV[sortedIndex[s1]][0];
    a2 = edgeV[sortedIndex[s1]][1];
    flag_intersection = 0;
    for( int s2 = 0; s2 < numSelected; ++s2 ){
      b1 = edgeV[selectedIndex[s2]][0];
      b2 = edgeV[selectedIndex[s2]][1];
      if( this->Intersection( a1, a2, b1, b2 ) ){
	flag_intersection = 1;
      }
    }
    if( flag_intersection == 0 ){ // normal // T1, T2
      selectedIndex[numSelected++] = sortedIndex[s1];
    }
  }

  /////// set fKKi, fKK, etc  //////
  fKKi.resize(fNN,fNN);
  fKKi = MatrixXi::Zero(fNN,fNN); 
  double alpha = 1.0;
  for( int s = 0; s < numSelected; ++s ){
    int i = edgeV[selectedIndex[s]][0];
    int j = edgeV[selectedIndex[s]][1];

    fKKi(i,i) += 1;
    fKKi(j,j) += 1;
    fKKi(i,j) -= 1;
    fKKi(j,i) -= 1;
  }

  for( int i = 0; i < fNN; ++i ){
    int num = 0;
    for( int j = 0; j < fNN; ++j ){
      if( fKKi(i,j) < 0 )
	fNei(i,num++) = j;
      assert( num < 20 );
    }
    fNei_Size(i) = num;
  }

  for( int s = 0; s < numSelected; ++s ){
    int i = edgeV[selectedIndex[s]][0];
    int j = edgeV[selectedIndex[s]][1];

    if( i >= fN || j >= fN ){ 
      fV_short(0) = i; // 輪郭を除く最小edge
      fV_short(1) = j;
      break;
    }
  }
  for( int s = numSelected-1; s >= 0; --s ){
    int i = edgeV[selectedIndex[s]][0];
    int j = edgeV[selectedIndex[s]][1];

    if( i >= fN || j >= fN ){ 
      fV_long(0) = i;
      fV_long(1) = j;
      break;
    }
  }

}

int TSearch::iix( int i ){ 
  if( i < 0 )
    i += fNN;
  else if ( i >= fNN )
    i -= fNN;
  return 2*i;
}

int TSearch::iiy( int i ){
  if( i < 0 )
    i += fNN;
  else if ( i >= fNN )
    i -= fNN;
  return 2*i+1;
}





bool TSearch::Check_Internal_link( int s1, int s2 ) 
{
  if( (fabs(s1-s2) == 1 || fabs(s1-s2) == fN-1) && s1 < fN && s2 < fN )
    return true;
  
  int count;

  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;

  // 輪郭と交差なし
  x1 = fWW(iix(s1)); // Inner
  y1 = fWW(iiy(s1));
  x2 = fWW(iix(s2));
  y2 = fWW(iiy(s2));
  for( int i = 0; i < fN; ++i ){
    x3 = fW(ix(i));
    y3 = fW(iy(i));
    x4 = fW(ix(i+1));
    y4 = fW(iy(i+1));

    tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
    td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
    value1 = tc*td;
    
    tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
    td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
    value2 = tc*td;

    if( value1 < 0 && value2 < 0 )
      return false;
  }  


  // return true; // outer (remove comment out) 

  // ポリゴン内部
  x1 = 0.5*(fWW(iix(s1))+fWW(iix(s2)));
  y1 = 0.5*(fWW(iiy(s1))+fWW(iiy(s2)));
  x2 = x1+10000.0;
  y2 = y1+10.0;

  count = 0;
  for( int i = 0; i < fN; ++i ){
    x3 = fW(ix(i));
    y3 = fW(iy(i));
    x4 = fW(ix(i+1));
    y4 = fW(iy(i+1));

    tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
    td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
    value1 = tc*td;
    
    tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
    td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
    value2 = tc*td;

    if( value1 < 0 && value2 < 0 )
      ++count;
  }  

  if( count % 2 == 1 )  
    return true;
  else 
    return false;
}


bool TSearch::Intersection( int a1, int a2, int b1, int b2 )
{
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;

  x1 = fWW(iix(a1)); // Inner
  y1 = fWW(iiy(a1));
  x2 = fWW(iix(a2));
  y2 = fWW(iiy(a2));
  x3 = fWW(iix(b1));
  y3 = fWW(iiy(b1));
  x4 = fWW(iix(b2));
  y4 = fWW(iiy(b2));

  tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
  td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
  value1 = tc*td;

  tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
  td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
  value2 = tc*td;

  if( value1 < -0.0 && value2 < -0.0 )
    return true;
  else
    return false;
}


void TSearch::SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest;
  double valueBest;
  int checked[50000];

  assert( Arg[0] < 99999999999.9 );
  assert( numOfArg < 50000 );

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





