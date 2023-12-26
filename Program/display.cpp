/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#include "display.h"
#endif

int open_flg_tiling = 0;
int d_num_tiling;

TDisplay::TDisplay()
{
  fNumDisplay = 1;
  fWidth = 600;
  fHeight = 600;
  fScale_T = 0.4;
}

TDisplay::~TDisplay()
{
}


void TDisplay::SetInit( int N1, int N1_in, VectorXd ww1, MatrixXi kk1i, int N2, int N2_in, VectorXd ww2, MatrixXi kk2i )
{
  fN1 = N1;
  fN1_in = N1_in;
  fNN1 = fN1+fN1_in;
  fWW1 = ww1;
  fW1 = fWW1.head(2*fN1);
  fKK1i = kk1i;

  fN2 = N2;
  fN2_in = N2_in;
  fNN2 = fN2+fN2_in;
  fWW2 = ww2;
  fW2 = fWW2.head(2*fN2);
  fKK2i = kk2i;
  
  fOpen_flag_tiling = 0;
  this->Set_range1();
  this->Set_range2();

  tOperator = new TOperator_ADD( fN1, fN1_in, fN2, fN2_in );
  tOperator->SetParameter();
  tOperator->SetInit( fW1, fWW1.tail(2*fN1_in), fKK1i, fW2, fWW2.tail(2*fN2_in), fKK2i );
}

void TDisplay::Set_range1()
{
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9;

  for( int i = 0; i < fN1; ++i ){
    if( fW1(2*i) < x_min )
      x_min = fW1(2*i);
    if( fW1(2*i) > x_max )
      x_max = fW1(2*i);
    if( fW1(2*i+1) < y_min )
      y_min = fW1(2*i+1);
    if( fW1(2*i+1) > y_max )
      y_max = fW1(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  fxs1 = (x_min+x_max)*0.5 - 0.6*length; 
  fxe1 = (x_min+x_max)*0.5 + 0.6*length; 
  fys1 = (y_min+y_max)*0.5 - 0.6*length; 
  fye1 = (y_min+y_max)*0.5 + 0.6*length;
}

void TDisplay::Set_range2()
{
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9;

  for( int i = 0; i < fN2; ++i ){
    if( fW2(2*i) < x_min )
      x_min = fW2(2*i);
    if( fW2(2*i) > x_max )
      x_max = fW2(2*i);
    if( fW2(2*i+1) < y_min )
      y_min = fW2(2*i+1);
    if( fW2(2*i+1) > y_max )
      y_max = fW2(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  fxs2 = (x_min+x_max)*0.5 - 0.6*length; 
  fxe2 = (x_min+x_max)*0.5 + 0.6*length; 
  fys2 = (y_min+y_max)*0.5 - 0.6*length; 
  fye2 = (y_min+y_max)*0.5 + 0.6*length;
}

/*
void TDisplay::Goal() 
{
  static int open_flg = 0;
  static int d_num = 0;


  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Goal" );
    g_window( d_num, fxs, fys, fxe, fye );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN; ++i )
    g_line( d_num, fW(ix(i)), fW(iy(i)), fW(ix(i+1)), fW(iy(i+1)), 3 );

  g_flush();
}
*/


void TDisplay::Goal_mesh1() 
{
  assert( fN1_in != 0 );
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Goal_mesh1" );
    g_window( d_num, fxs1, fys1, fxe1, fye1 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);    
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN1-1; ++i ){
    for( int j = i+1; j < fNN1; ++j ){
      if( fKK1i(i,j) != 0 ){
	g_line( d_num, fWW1(iix1(i)), fWW1(iiy1(i)), fWW1(iix1(j)), fWW1(iiy1(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN1; ++i )
    g_line( d_num, fWW1(ix1(i)), fWW1(iy1(i)), fWW1(ix1(i+1)), fWW1(iy1(i+1)), 3 );

  g_flush();
}

void TDisplay::Goal_mesh2() 
{
  assert( fN2_in != 0 );
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Goal_mesh2" );
    g_window( d_num, fxs2, fys2, fxe2, fye2 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);    
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN2-1; ++i ){
    for( int j = i+1; j < fNN2; ++j ){
      if( fKK2i(i,j) != 0 ){
	g_line( d_num, fWW2(iix2(i)), fWW2(iiy2(i)), fWW2(iix2(j)), fWW2(iiy2(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN2; ++i )
    g_line( d_num, fWW2(ix2(i)), fWW2(iy2(i)), fWW2(ix2(i+1)), fWW2(iy2(i+1)), 3 );

  g_flush();
}

void TDisplay::Tile_mesh1( VectorXd uu ) 
{
  assert( uu.rows() == 2*fNN1 );
  static int open_flg = 0;
  static int d_num = 0;

  // 目標図形にfitするようscale-rotation
  tOperator->Scale_Rotate_UU1( uu );

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Tile mesh1" );
    g_window( d_num, fxs1, fys1, fxe1, fye1 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN1-1; ++i ){
    for( int j = i+1; j < fNN1; ++j ){
      if( fKK1i(i,j) != 0 ){
	g_line( d_num, uu(iix1(i)), uu(iiy1(i)), uu(iix1(j)), uu(iiy1(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN1; ++i )
    g_line( d_num, uu(ix1(i)), uu(iy1(i)), uu(ix1(i+1)), uu(iy1(i+1)), 3 );

  g_flush();
}

void TDisplay::Tile_mesh2( VectorXd uu ) 
{
  assert( uu.rows() == 2*fNN2 );
  static int open_flg = 0;
  static int d_num = 0;

  // 目標図形にfitするようscale-rotation
  tOperator->Scale_Rotate_UU2( uu );

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Tile mesh2" );
    g_window( d_num, fxs2, fys2, fxe2, fye2 );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN2-1; ++i ){
    for( int j = i+1; j < fNN2; ++j ){
      if( fKK2i(i,j) != 0 ){
	g_line( d_num, uu(iix2(i)), uu(iiy2(i)), uu(iix2(j)), uu(iiy2(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN2; ++i )
    g_line( d_num, uu(ix2(i)), uu(iy2(i)), uu(ix2(i+1)), uu(iy2(i+1)), 3 );

  g_flush();
}


void TDisplay::Dihedral( VectorXd u, VectorXd ua, int st0, int nv, int ffn[], VectorXd uu1, VectorXd uu2 ) 
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends
  fN = n;

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 
  for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(2*i);
    if( uu1(2*i) > x_max )
      x_max = uu1(2*i);
    if( uu1(2*i+1) < y_min )
      y_min = uu1(2*i+1);
    if( uu1(2*i+1) > y_max )
      y_max = uu1(2*i+1);
  }
  for( int i = 0; i < fNN2; ++i ){
    if( uu2(2*i) < x_min )
      x_min = uu2(2*i);
    if( uu2(2*i) > x_max )
      x_max = uu2(2*i);
    if( uu2(2*i+1) < y_min )
      y_min = uu2(2*i+1);
    if( uu2(2*i+1) > y_max )
      y_max = uu2(2*i+1);
  }

  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);

  for( int i = 0; i < fNN1; ++i ){
    uu1(iix1(i)) -= x_mid;
    uu1(iiy1(i)) -= y_mid;
  }
  for( int i = 0; i < fNN2; ++i ){
    uu2(iix2(i)) -= x_mid;
    uu2(iiy2(i)) -= y_mid;
  }
  for( int i = 0; i < fN; ++i ){
    u(ix(i)) -= x_mid;
    u(iy(i)) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  double S = tOperator->Cal_Area1( uu1 ) + tOperator->Cal_Area2( uu2 );
  
  uu1 *= sqrt(fWidth*fHeight/S)*fScale_T;
  uu2 *= sqrt(fWidth*fHeight/S)*fScale_T;
  u *= sqrt(fWidth*fHeight/S)*fScale_T;
  ua *= sqrt(fWidth*fHeight/S)*fScale_T;
  
  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Dihedral" );
    g_window( d_num, -0.5*fWidth, -0.5*fHeight, 0.5*fWidth, 0.5*fHeight ); 
    g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-20-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);


  /*
  g_setColor(d_num,2);
  for( int i = 0; i < n-1; ++i )
    g_line( d_num, u(2*i), u(2*i+1), u(2*(i+1)), u(2*(i+1)+1), 4 );
  g_line( d_num, u(2*(n-1)), u(2*(n-1)+1), u(0), u(1), 4 );
  */

  g_setColor(d_num,2);
  for( int i = 0; i < fN1-na+1; ++i )
    g_line( d_num, u(ix(i+st0)), u(iy(i+st0)), u(ix(i+1+st0)), u(iy(i+1+st0)), 4 );
  g_setColor(d_num,3);
  for( int i = fN1-na+1; i < n; ++i )
    g_line( d_num, u(ix(i+st0)), u(iy(i+st0)), u(ix(i+1+st0)), u(iy(i+1+st0)), 4 );
  
  g_setColor(d_num,4);
  for( int i = 0; i < na-1; ++i )
    g_line( d_num, ua(2*i), ua(2*i+1), ua(2*(i+1)), ua(2*(i+1)+1), 4 );

  
  g_setColor(d_num,1);
  int h = 0;
  double boxL = 0.5;
  char str[80];
  for( int v = 0; v < nv; ++v ){
    g_box( d_num, u(2*h)-boxL, u(2*h+1)-boxL, u(2*h)+boxL, u(2*h+1)+boxL, 2 );
    sprintf( str, "%d", v );
    g_string( d_num, u(2*h), u(2*h+1), str );
    h += ffn[v]+1;
  }
  assert( h == n );
  
  g_flush();
}


void TDisplay::Dihedral2( VectorXd u, VectorXd ua, int st0, VectorXd uu1, VectorXd uu2 )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends
  fN = n;

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 
  for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(2*i);
    if( uu1(2*i) > x_max )
      x_max = uu1(2*i);
    if( uu1(2*i+1) < y_min )
      y_min = uu1(2*i+1);
    if( uu1(2*i+1) > y_max )
      y_max = uu1(2*i+1);
  }
  for( int i = 0; i < fNN2; ++i ){
    if( uu2(2*i) < x_min )
      x_min = uu2(2*i);
    if( uu2(2*i) > x_max )
      x_max = uu2(2*i);
    if( uu2(2*i+1) < y_min )
      y_min = uu2(2*i+1);
    if( uu2(2*i+1) > y_max )
      y_max = uu2(2*i+1);
  }

  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);

  for( int i = 0; i < fNN1; ++i ){
    uu1(iix1(i)) -= x_mid;
    uu1(iiy1(i)) -= y_mid;
  }
  for( int i = 0; i < fNN2; ++i ){
    uu2(iix2(i)) -= x_mid;
    uu2(iiy2(i)) -= y_mid;
  }
  for( int i = 0; i < fN; ++i ){
    u(ix(i)) -= x_mid;
    u(iy(i)) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }
  
  double S = tOperator->Cal_Area1( uu1 ) + tOperator->Cal_Area2( uu2 );
  
  uu1 *= sqrt(fWidth*fHeight/S)*fScale_T;
  uu2 *= sqrt(fWidth*fHeight/S)*fScale_T;
  u *= sqrt(fWidth*fHeight/S)*fScale_T;
  ua *= sqrt(fWidth*fHeight/S)*fScale_T;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Dihedral" );
    g_window( d_num, -0.5*fWidth, -0.5*fHeight, 0.5*fWidth, 0.5*fHeight ); 
    g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-20-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);

  /*
  fN = n; // ix(), iy()のために必要
  g_setColor(d_num,2);
  for( int i = 0; i < fN1-na+1; ++i )
    g_line( d_num, u(ix(i+st0)), u(iy(i+st0)), u(ix(i+1+st0)), u(iy(i+1+st0)), 4 );
  g_setColor(d_num,3);
  for( int i = fN1-na+1; i < n; ++i )
    g_line( d_num, u(ix(i+st0)), u(iy(i+st0)), u(ix(i+1+st0)), u(iy(i+1+st0)), 4 );

  
  g_setColor(d_num,4);
  for( int i = 0; i < na-1; ++i )
    g_line( d_num, ua(2*i), ua(2*i+1), ua(2*(i+1)), ua(2*(i+1)+1), 4 );
  */


  // mesh
  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN1-1; ++i ){
    for( int j = i+1; j < fNN1; ++j ){
      if( fKK1i(i,j) != 0 ){
	g_line( d_num, uu1(iix1(i)), uu1(iiy1(i)), uu1(iix1(j)), uu1(iiy1(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN1; ++i )
    g_line( d_num, uu1(ix1(i)), uu1(iy1(i)), uu1(ix1(i+1)), uu1(iy1(i+1)), 3 );

  /*
  double boxL = 2.0;
  g_setColor(d_num,1); // black
  for( int i = 0; i < fN1; ++i )
    g_box( d_num, uu1(2*i)-boxL, uu1(2*i+1)-boxL, uu1(2*i)+boxL, uu1(2*i+1)+boxL, 2 );
  */

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN2-1; ++i ){
    for( int j = i+1; j < fNN2; ++j ){
      if( fKK2i(i,j) != 0 ){
	g_line( d_num, uu2(iix2(i)), uu2(iiy2(i)), uu2(iix2(j)), uu2(iiy2(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN2; ++i )
    g_line( d_num, uu2(ix2(i)), uu2(iy2(i)), uu2(ix2(i+1)), uu2(iy2(i+1)), 3 );

  /*
  g_setColor(d_num,1); // black  
  for( int i = 0; i < fN2; ++i )
    g_box( d_num, uu2(2*i)-boxL, uu2(2*i+1)-boxL, uu2(2*i)+boxL, uu2(2*i+1)+boxL, 2 );
  */
  

  g_flush();
}


void TDisplay::Tiling_IH4( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd ua_ori_1(2*na);

  for( int i = 0; i < n; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }
  u_ori_1 = - u_ori_0;

  for( int i = 0; i < na; ++i ){ // new
    ua_ori_0(2*i) = ua(2*i);
    ua_ori_0(2*i+1) = ua(2*i+1);
  }
  ua_ori_1 = - ua_ori_0;

  diff(0) = u_ori_1(2*fi[3]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[3]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){  // new
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[n][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 %2 == 0 )
	g_setColor(d_num,10); 
      else
	g_setColor(d_num,15);
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_0(2*i)+dx, ua_ori_0(2*i+1)+dy, ua_ori_0(2*(i+1))+dx, ua_ori_0(2*(i+1)+1)+dy, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 %2 == 0 )
	g_setColor(d_num,19); 
      else
	g_setColor(d_num,29);
      
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new 
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_1(2*i)+dx, ua_ori_1(2*i+1)+dy, ua_ori_1(2*(i+1))+dx, ua_ori_1(2*(i+1)+1)+dy, 4 );
      }

    } 
  }

  g_flush();
}

void TDisplay::Tiling_IH5( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd u_ori_2(2*n);
  VectorXd u_ori_3(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd ua_ori_1(2*na);
  VectorXd ua_ori_2(2*na);
  VectorXd ua_ori_3(2*na);

  u_ori_0 = u;
  for( int i=0; i < n; ++i ){ 
    u_ori_1(2*i) = u_ori_0(2*i);
    u_ori_1(2*i+1) = -u_ori_0(2*i+1);
    u_ori_2(2*i) =  -u_ori_0(2*i);
    u_ori_2(2*i+1) = u_ori_0(2*i+1);
    u_ori_3(2*i) = -u_ori_0(2*i);
    u_ori_3(2*i+1) = -u_ori_0(2*i+1);
  }

  ua_ori_0 = ua;
  for( int i=0; i < na; ++i ){ 
    ua_ori_1(2*i) = ua_ori_0(2*i);
    ua_ori_1(2*i+1) = -ua_ori_0(2*i+1);
    ua_ori_2(2*i) =  -ua_ori_0(2*i);
    ua_ori_2(2*i+1) = ua_ori_0(2*i+1);
    ua_ori_3(2*i) = - -ua_ori_0(2*i);
    ua_ori_3(2*i+1) = - ua_ori_0(2*i+1);
  }

  diff(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[3]); 
  diff(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[3]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){    // new
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_2(2*fi[5]) - u_ori_1(2*fi[0]); 
  diff(1) = u_ori_2(2*fi[5]+1) - u_ori_1(2*fi[0]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){   // new
    ua_ori_2(2*i) -= diff(0);
    ua_ori_2(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_3(2*fi[2]) - u_ori_2(2*fi[1]); 
  diff(1) = u_ori_3(2*fi[2]+1) - u_ori_2(2*fi[1]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){   // new
    ua_ori_3(2*i) -= diff(0);
    ua_ori_3(2*i+1) -= diff(1);
   }
  
  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_3(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_3(2*fi[5]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[n][2];

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,10); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_0(2*i)+dx;
	double y1 = ua_ori_0(2*i+1)+dy;
	double x2 = ua_ori_0(2*(i+1))+dx;
	double y2 = ua_ori_0(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,15); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_1(2*i)+dx;
	double y1 = ua_ori_1(2*i+1)+dy;
	double x2 = ua_ori_1(2*(i+1))+dx;
	double y2 = ua_ori_1(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }
      

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,19); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_2(2*i)+dx;
	double y1 = ua_ori_2(2*i+1)+dy;
	double x2 = ua_ori_2(2*(i+1))+dx;
	double y2 = ua_ori_2(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }
      

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,29); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_3(2*i)+dx;
	double y1 = ua_ori_3(2*i+1)+dy;
	double x2 = ua_ori_3(2*(i+1))+dx;
	double y2 = ua_ori_3(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }      
    }
  }

  g_flush();
}


void TDisplay::Tiling_IH6( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd u_ori_2(2*n);
  VectorXd u_ori_3(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd ua_ori_1(2*na);
  VectorXd ua_ori_2(2*na);
  VectorXd ua_ori_3(2*na);

  u_ori_0 = u;
  for( int i=0; i < n; ++i ){ 
    u_ori_1(2*i) = -u_ori_0(2*i);
    u_ori_1(2*i+1) = u_ori_0(2*i+1);
    u_ori_2(2*i) =  u_ori_0(2*i);
    u_ori_2(2*i+1) = -u_ori_0(2*i+1);
    u_ori_3(2*i) = - u_ori_0(2*i);
    u_ori_3(2*i+1) = - u_ori_0(2*i+1);
  }

  ua_ori_0 = ua;
  for( int i=0; i < na; ++i ){ 
    ua_ori_1(2*i) = -ua_ori_0(2*i);
    ua_ori_1(2*i+1) = ua_ori_0(2*i+1);
    ua_ori_2(2*i) =  ua_ori_0(2*i);
    ua_ori_2(2*i+1) = -ua_ori_0(2*i+1);
    ua_ori_3(2*i) = - ua_ori_0(2*i);
    ua_ori_3(2*i+1) = - ua_ori_0(2*i+1);
  }


  diff(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){    // new
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_2(2*fi[0]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_2(2*fi[0]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){   // new
    ua_ori_2(2*i) -= diff(0);
    ua_ori_2(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_3(2*fi[3]) - u_ori_0(2*fi[4]); 
  diff(1) = u_ori_3(2*fi[3]+1) - u_ori_0(2*fi[4]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){   // new
    ua_ori_3(2*i) -= diff(0);
    ua_ori_3(2*i+1) -= diff(1);
   }
  

  diff_0(0) = u_ori_3(2*fi[0]) - u_ori_0(2*fi[5]);
  diff_0(1) = u_ori_3(2*fi[0]+1) - u_ori_0(2*fi[5]+1);
  diff_1(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[5]);
  diff_1(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[5]+1);

  double pointsD[n][2];

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,10); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_0(2*i)+dx;
	double y1 = ua_ori_0(2*i+1)+dy;
	double x2 = ua_ori_0(2*(i+1))+dx;
	double y2 = ua_ori_0(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,15); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_1(2*i)+dx;
	double y1 = ua_ori_1(2*i+1)+dy;
	double x2 = ua_ori_1(2*(i+1))+dx;
	double y2 = ua_ori_1(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }
      

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,19); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_2(2*i)+dx;
	double y1 = ua_ori_2(2*i+1)+dy;
	double x2 = ua_ori_2(2*(i+1))+dx;
	double y2 = ua_ori_2(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }
      

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD[i][1] = u_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num,29); 
      g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      // new
      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	double x1 = ua_ori_3(2*i)+dx;
	double y1 = ua_ori_3(2*i+1)+dy;
	double x2 = ua_ori_3(2*(i+1))+dx;
	double y2 = ua_ori_3(2*(i+1)+1)+dy;

	g_line( d_num, x1, y1, x2, y2, 4 );
      }      
    }
  }

  g_flush();
}



void TDisplay::Tiling_IH1( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd ua_ori_0(2*na);


  for( int i = 0; i < n; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }

  for( int i = 0; i < na; ++i ){
    ua_ori_0(2*i) = ua(2*i);
    ua_ori_0(2*i+1) = ua(2*i+1);
  }

  diff_0(0) = u_ori_0(2*fi[2]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[2]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[n][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      // g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_0(2*i)+dx, ua_ori_0(2*i+1)+dy, ua_ori_0(2*(i+1))+dx, ua_ori_0(2*(i+1)+1)+dy, 4 );
      }
    } 
  }

  g_flush();
}


void TDisplay::Tiling_IH2( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd u_ori_1(2*n);
  VectorXd ua_ori_1(2*na);


  for( int i = 0; i < n; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_0(2*i) = ua(2*i);
    ua_ori_0(2*i+1) = ua(2*i+1);
  }

  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) = u_ori_0(2*i);
    u_ori_1(2*i+1) = -u_ori_0(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) = ua_ori_0(2*i);
    ua_ori_1(2*i+1) = -ua_ori_0(2*i+1);
  }

  diff(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[3]); 
  diff(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[3]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[5]) - u_ori_0(2*fi[1]);
  diff_0(1) = u_ori_0(2*fi[5]+1) - u_ori_0(2*fi[1]+1);
  diff_1(0) = u_ori_1(2*fi[0]) - u_ori_0(2*fi[1]);
  diff_1(1) = u_ori_1(2*fi[0]+1) - u_ori_0(2*fi[1]+1);

  double pointsD[n][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      // g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_0(2*i)+dx, ua_ori_0(2*i+1)+dy, ua_ori_0(2*(i+1))+dx, ua_ori_0(2*(i+1)+1)+dy, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      // g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_1(2*i)+dx, ua_ori_1(2*i+1)+dy, ua_ori_1(2*(i+1))+dx, ua_ori_1(2*(i+1)+1)+dy, 4 );
      }
    } 
  }

  g_flush();
}


void TDisplay::Tiling_IH3( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd u_ori_1(2*n);
  VectorXd ua_ori_1(2*na);


  for( int i = 0; i < n; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_0(2*i) = ua(2*i);
    ua_ori_0(2*i+1) = ua(2*i+1);
  }

  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) = -u_ori_0(2*i);
    u_ori_1(2*i+1) = u_ori_0(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) = -ua_ori_0(2*i);
    ua_ori_1(2*i+1) = ua_ori_0(2*i+1);
  }

  diff(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[1]); 
  diff(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[1]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[n][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      // g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_0(2*i)+dx, ua_ori_0(2*i+1)+dy, ua_ori_0(2*(i+1))+dx, ua_ori_0(2*(i+1)+1)+dy, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      // g_polygonfill( d_num, pointsD, n );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_1(2*i)+dx, ua_ori_1(2*i+1)+dy, ua_ori_1(2*(i+1))+dx, ua_ori_1(2*(i+1)+1)+dy, 4 );
      }
    } 
  }

  g_flush();
}

void TDisplay::Tiling_IH27( VectorXd& u, VectorXd& ua, int* fn, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( u(2*i) < x_min )
      x_min = u(2*i);
    if( u(2*i) > x_max )
      x_max = u(2*i);
    if( u(2*i+1) < y_min )
      y_min = u(2*i+1);
    if( u(2*i+1) > y_max )
      y_max = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    if( ua(2*i) < x_min )
      x_min = ua(2*i);
    if( ua(2*i) > x_max )
      x_max = ua(2*i);
    if( ua(2*i+1) < y_min )
      y_min = ua(2*i+1);
    if( ua(2*i+1) > y_max )
      y_max = ua(2*i+1);
  }


  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);
  for( int i = 0; i < n; ++i ){
    u(2*i) -= x_mid;
    u(2*i+1) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;

  u *= 100.0/length;
  ua *= 100.0/length;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, -80, -80, 80, 80 );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd ua_ori_0(2*na);
  VectorXd u_ori_1(2*n);
  VectorXd ua_ori_1(2*na);
  VectorXd u_ori_2(2*n);
  VectorXd ua_ori_2(2*na);
  VectorXd u_ori_3(2*n);
  VectorXd ua_ori_3(2*na);


  for( int i = 0; i < n; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_0(2*i) = ua(2*i);
    ua_ori_0(2*i+1) = ua(2*i+1);
  }

  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) = -u_ori_0(2*i);
    u_ori_1(2*i+1) = u_ori_0(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) = -ua_ori_0(2*i);
    ua_ori_1(2*i+1) = ua_ori_0(2*i+1);
  }

  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) = u_ori_0(2*i);
    u_ori_2(2*i+1) = -u_ori_0(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_2(2*i) = ua_ori_0(2*i);
    ua_ori_2(2*i+1) = -ua_ori_0(2*i+1);
  }

  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) = -u_ori_0(2*i);
    u_ori_3(2*i+1) = -u_ori_0(2*i+1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_3(2*i) = -ua_ori_0(2*i);
    ua_ori_3(2*i+1) = -ua_ori_0(2*i+1);
  }
  

  diff(0) = u_ori_1(2*fi[4]) - u_ori_0(2*fi[1]); 
  diff(1) = u_ori_1(2*fi[4]+1) - u_ori_0(2*fi[1]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_1(2*i) -= diff(0);
    ua_ori_1(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_2(2*fi[1]) - u_ori_0(2*fi[4]); 
  diff(1) = u_ori_2(2*fi[1]+1) - u_ori_0(2*fi[4]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_2(2*i) -= diff(0);
    ua_ori_2(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_3(2*fi[3]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_3(2*fi[3]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < na; ++i ){
    ua_ori_3(2*i) -= diff(0);
    ua_ori_3(2*i+1) -= diff(1);
  }

  
  diff_0(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_2(2*fi[3]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_2(2*fi[3]+1) - u_ori_0(2*fi[0]+1);
  
  double pointsD[n][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
   
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_0(2*i)+dx, ua_ori_0(2*i+1)+dy, ua_ori_0(2*(i+1))+dx, ua_ori_0(2*(i+1)+1)+dy, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
     
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_1(2*i)+dx, ua_ori_1(2*i+1)+dy, ua_ori_1(2*(i+1))+dx, ua_ori_1(2*(i+1)+1)+dy, 4 );
      }

      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
     
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_2(2*i)+dx, ua_ori_2(2*i+1)+dy, ua_ori_2(2*(i+1))+dx, ua_ori_2(2*(i+1)+1)+dy, 4 );
      }
      
      // polygon
      for( int i = 0; i < n ; ++i ){
	pointsD[i][0] = u_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
     
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, n, 4 );

      g_setColor(d_num,4);
      for( int i = 0; i < na-1; ++i ){
	double dx = diff_0(0)*k1+diff_1(0)*k2;
	double dy = diff_0(1)*k1+diff_1(1)*k2;
	g_line( d_num, ua_ori_3(2*i)+dx, ua_ori_3(2*i+1)+dy, ua_ori_3(2*(i+1))+dx, ua_ori_3(2*(i+1)+1)+dy, 4 );
      }

      
    } 
  }

  g_flush();
}


void TDisplay::Tiling_IH4( VectorXd u, VectorXd ua, int* fn, int nv, VectorXd uu1, VectorXd uu2 )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends
  fN = n;

  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 
  for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(2*i);
    if( uu1(2*i) > x_max )
      x_max = uu1(2*i);
    if( uu1(2*i+1) < y_min )
      y_min = uu1(2*i+1);
    if( uu1(2*i+1) > y_max )
      y_max = uu1(2*i+1);
  }
  for( int i = 0; i < fNN2; ++i ){
    if( uu2(2*i) < x_min )
      x_min = uu2(2*i);
    if( uu2(2*i) > x_max )
      x_max = uu2(2*i);
    if( uu2(2*i+1) < y_min )
      y_min = uu2(2*i+1);
    if( uu2(2*i+1) > y_max )
      y_max = uu2(2*i+1);
  }

  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);

  for( int i = 0; i < fNN1; ++i ){
    uu1(iix1(i)) -= x_mid;
    uu1(iiy1(i)) -= y_mid;
  }
  for( int i = 0; i < fNN2; ++i ){
    uu2(iix2(i)) -= x_mid;
    uu2(iiy2(i)) -= y_mid;
  }
  for( int i = 0; i < fN; ++i ){
    u(ix(i)) -= x_mid;
    u(iy(i)) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }
  
  double S = tOperator->Cal_Area1( uu1 ) + tOperator->Cal_Area2( uu2 );
  
  uu1 *= sqrt(fWidth*fHeight/S)*fScale_T;
  uu2 *= sqrt(fWidth*fHeight/S)*fScale_T;
  u *= sqrt(fWidth*fHeight/S)*fScale_T;
  ua *= sqrt(fWidth*fHeight/S)*fScale_T;

  if (open_flg_tiling == 0) 
  {
    d_num_tiling = fNumDisplay;
    ++fNumDisplay;
    open_flg_tiling = 1;

    g_open((char*)"");
    g_open_window( d_num_tiling, 800, 800, (char*)"Tiling" );
    g_window( d_num_tiling, -400, -400, 400, 400 );
    //g_setFont(d_num_tiling,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num_tiling); 
    usleep(100000);
  }
  g_clearWindow(d_num_tiling);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd u1_ori_0(2*fN1);
  VectorXd u1_ori_1(2*fN1);
  VectorXd u2_ori_0(2*fN2);
  VectorXd u2_ori_1(2*fN2);

  u_ori_0 = u;
  u_ori_1 = - u_ori_0;
  u1_ori_0 = uu1.head(2*fN1);
  u1_ori_1 = - u1_ori_0;
  u2_ori_0 = uu2.head(2*fN2);
  u2_ori_1 = - u2_ori_0;

  diff(0) = u_ori_1(2*fi[3]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[3]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){ 
    u1_ori_1(2*i) -= diff(0);
    u1_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){ 
    u2_ori_1(2*i) -= diff(0);
    u2_ori_1(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[0]+1);

  double pointsD1[fN1][2];
  double pointsD2[fN2][2];  
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_0(ix1(i))+diff_0(0)*k1+diff_1(0)*k2;
	pointsD1[i][1] = u1_ori_0(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_0(ix2(i))+diff_0(0)*k1+diff_1(0)*k2;
	pointsD2[i][1] = u2_ori_0(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }

      g_setColor(d_num_tiling,10); 
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 15
      g_polygonfill( d_num_tiling, pointsD2, fN2 );      
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_1(ix1(i))+diff_0(0)*k1+diff_1(0)*k2;
	pointsD1[i][1] = u1_ori_1(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_1(ix2(i))+diff_0(0)*k1+diff_1(0)*k2;
	pointsD2[i][1] = u2_ori_1(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      g_setColor(d_num_tiling,10); // 19
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 29
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );
    }
  }

  g_flush();
}

void TDisplay::Tiling_IH5( VectorXd u, VectorXd ua, int* fn, int nv, VectorXd uu1, VectorXd uu2 )
{
   static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends
  fN = n;
  
  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 
  for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(2*i);
    if( uu1(2*i) > x_max )
      x_max = uu1(2*i);
    if( uu1(2*i+1) < y_min )
      y_min = uu1(2*i+1);
    if( uu1(2*i+1) > y_max )
      y_max = uu1(2*i+1);
  }
  for( int i = 0; i < fNN2; ++i ){
    if( uu2(2*i) < x_min )
      x_min = uu2(2*i);
    if( uu2(2*i) > x_max )
      x_max = uu2(2*i);
    if( uu2(2*i+1) < y_min )
      y_min = uu2(2*i+1);
    if( uu2(2*i+1) > y_max )
      y_max = uu2(2*i+1);
  }

  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);

  for( int i = 0; i < fNN1; ++i ){
    uu1(iix1(i)) -= x_mid;
    uu1(iiy1(i)) -= y_mid;
  }
  for( int i = 0; i < fNN2; ++i ){
    uu2(iix2(i)) -= x_mid;
    uu2(iiy2(i)) -= y_mid;
  }
  for( int i = 0; i < fN; ++i ){
    u(ix(i)) -= x_mid;
    u(iy(i)) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }
  
  double S = tOperator->Cal_Area1( uu1 ) + tOperator->Cal_Area2( uu2 );
  
  uu1 *= sqrt(fWidth*fHeight/S)*fScale_T;
  uu2 *= sqrt(fWidth*fHeight/S)*fScale_T;
  u *= sqrt(fWidth*fHeight/S)*fScale_T;
  ua *= sqrt(fWidth*fHeight/S)*fScale_T;

  if (open_flg_tiling == 0) 
  {
    d_num_tiling = fNumDisplay;
    ++fNumDisplay;
    open_flg_tiling = 1;

    g_open((char*)"");
    g_open_window( d_num_tiling, 800, 800, (char*)"Tiling" );
    g_window( d_num_tiling, -400, -400, 400, 400 );
    //g_setFont(d_num_tiling,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num_tiling); 
    usleep(100000);
  }
  g_clearWindow(d_num_tiling);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd u_ori_2(2*n);
  VectorXd u_ori_3(2*n);
  VectorXd u1_ori_0(2*fN1);
  VectorXd u1_ori_1(2*fN1);
  VectorXd u1_ori_2(2*fN1);
  VectorXd u1_ori_3(2*fN1);
  VectorXd u2_ori_0(2*fN2);
  VectorXd u2_ori_1(2*fN2);
  VectorXd u2_ori_2(2*fN2);
  VectorXd u2_ori_3(2*fN2);

  u_ori_0 = u;
  for( int i = 0; i < n; ++i ){ 
    u_ori_1(2*i) = u_ori_0(2*i);
    u_ori_1(2*i+1) = -u_ori_0(2*i+1);
    u_ori_2(2*i) =  -u_ori_0(2*i);
    u_ori_2(2*i+1) = u_ori_0(2*i+1);
    u_ori_3(2*i) = -u_ori_0(2*i);
    u_ori_3(2*i+1) = - u_ori_0(2*i+1);
  }
  u1_ori_0 = uu1.head(2*fN1);
  for( int i = 0; i < fN1; ++i ){ 
    u1_ori_1(2*i) = u1_ori_0(2*i);
    u1_ori_1(2*i+1) = -u1_ori_0(2*i+1);
    u1_ori_2(2*i) =  -u1_ori_0(2*i);
    u1_ori_2(2*i+1) = u1_ori_0(2*i+1);
    u1_ori_3(2*i) = -u1_ori_0(2*i);
    u1_ori_3(2*i+1) = -u1_ori_0(2*i+1);
  }
  u2_ori_0 = uu2.head(2*fN2);
  for( int i = 0; i < fN2; ++i ){ 
    u2_ori_1(2*i) = u2_ori_0(2*i);
    u2_ori_1(2*i+1) = -u2_ori_0(2*i+1);
    u2_ori_2(2*i) =  -u2_ori_0(2*i);
    u2_ori_2(2*i+1) = u2_ori_0(2*i+1);
    u2_ori_3(2*i) = -u2_ori_0(2*i);
    u2_ori_3(2*i+1) = -u2_ori_0(2*i+1);
  }

  diff(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[3]); 
  diff(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[3]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){
    u1_ori_1(2*i) -= diff(0);
    u1_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){
    u2_ori_1(2*i) -= diff(0);
    u2_ori_1(2*i+1) -= diff(1);
  }

  diff(0) = u_ori_2(2*fi[5]) - u_ori_1(2*fi[0]); 
  diff(1) = u_ori_2(2*fi[5]+1) - u_ori_1(2*fi[0]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){
    u1_ori_2(2*i) -= diff(0);
    u1_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){
    u2_ori_2(2*i) -= diff(0);
    u2_ori_2(2*i+1) -= diff(1);
  }

  diff(0) = u_ori_3(2*fi[2]) - u_ori_2(2*fi[1]); 
  diff(1) = u_ori_3(2*fi[2]+1) - u_ori_2(2*fi[1]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){
    u1_ori_3(2*i) -= diff(0);
    u1_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){
    u2_ori_3(2*i) -= diff(0);
    u2_ori_3(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_3(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_3(2*fi[5]+1) - u_ori_0(2*fi[0]+1);
  
  double pointsD1[fN1][2];
  double pointsD2[fN2][2];

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD1[i][1] = u1_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD2[i][1] = u2_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      
      g_setColor(d_num_tiling,10); 
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 15
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD1[i][1] = u1_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD2[i][1] = u2_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }

      ////////////////////
      g_setColor(d_num_tiling,10);  // 19
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19);  // 29 
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD1[i][1] = u1_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD2[i][1] = u2_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }

      g_setColor(d_num_tiling,10);  // 29
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19);  // 19
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD1[i][1] = u1_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD2[i][1] = u2_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }

      g_setColor(d_num_tiling,10);  // 15
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 10
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1);
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );
    }
  }
  
  g_flush();
}


void TDisplay::Tiling_IH6( VectorXd u, VectorXd ua, int* fn, int nv, VectorXd uu1, VectorXd uu2 )
{
   static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int n = u.rows()/2;
  int na = ua.rows()/2; // including both ends
  fN = n;
  
  double x_min, x_max, y_min, y_max, x_mid, y_mid, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 
  for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(2*i);
    if( uu1(2*i) > x_max )
      x_max = uu1(2*i);
    if( uu1(2*i+1) < y_min )
      y_min = uu1(2*i+1);
    if( uu1(2*i+1) > y_max )
      y_max = uu1(2*i+1);
  }
  for( int i = 0; i < fNN2; ++i ){
    if( uu2(2*i) < x_min )
      x_min = uu2(2*i);
    if( uu2(2*i) > x_max )
      x_max = uu2(2*i);
    if( uu2(2*i+1) < y_min )
      y_min = uu2(2*i+1);
    if( uu2(2*i+1) > y_max )
      y_max = uu2(2*i+1);
  }

  x_mid = 0.5*(x_min+x_max);
  y_mid = 0.5*(y_min+y_max);

  for( int i = 0; i < fNN1; ++i ){
    uu1(iix1(i)) -= x_mid;
    uu1(iiy1(i)) -= y_mid;
  }
  for( int i = 0; i < fNN2; ++i ){
    uu2(iix2(i)) -= x_mid;
    uu2(iiy2(i)) -= y_mid;
  }
  for( int i = 0; i < fN; ++i ){
    u(ix(i)) -= x_mid;
    u(iy(i)) -= y_mid;
  }
  for( int i = 0; i < na; ++i ){
    ua(2*i) -= x_mid;
    ua(2*i+1) -= y_mid;
  }
  
  double S = tOperator->Cal_Area1( uu1 ) + tOperator->Cal_Area2( uu2 );
  
  uu1 *= sqrt(fWidth*fHeight/S)*fScale_T;
  uu2 *= sqrt(fWidth*fHeight/S)*fScale_T;
  u *= sqrt(fWidth*fHeight/S)*fScale_T;
  ua *= sqrt(fWidth*fHeight/S)*fScale_T;

  if (open_flg_tiling == 0) 
  {
    d_num_tiling = fNumDisplay;
    ++fNumDisplay;
    open_flg_tiling = 1;

    g_open((char*)"");
    g_open_window( d_num_tiling, 800, 800, (char*)"Tiling" );
    g_window( d_num_tiling, -400, -400, 400, 400 );
    //g_setFont(d_num_tiling,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num_tiling); 
    usleep(100000);
  }
  g_clearWindow(d_num_tiling);


  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  
  VectorXd u_ori_0(2*n);
  VectorXd u_ori_1(2*n);
  VectorXd u_ori_2(2*n);
  VectorXd u_ori_3(2*n);
  VectorXd u1_ori_0(2*fN1);
  VectorXd u1_ori_1(2*fN1);
  VectorXd u1_ori_2(2*fN1);
  VectorXd u1_ori_3(2*fN1);
  VectorXd u2_ori_0(2*fN2);
  VectorXd u2_ori_1(2*fN2);
  VectorXd u2_ori_2(2*fN2);
  VectorXd u2_ori_3(2*fN2);

  u_ori_0 = u;
  for( int i=0; i < n; ++i ){ 
    u_ori_1(2*i) = -u_ori_0(2*i);
    u_ori_1(2*i+1) = u_ori_0(2*i+1);
    u_ori_2(2*i) =  u_ori_0(2*i);
    u_ori_2(2*i+1) = -u_ori_0(2*i+1);
    u_ori_3(2*i) = - u_ori_0(2*i);
    u_ori_3(2*i+1) = - u_ori_0(2*i+1);
  }
  u1_ori_0 = uu1.head(2*fN1);
  for( int i=0; i < fN1; ++i ){ 
    u1_ori_1(2*i) = -u1_ori_0(2*i);
    u1_ori_1(2*i+1) = u1_ori_0(2*i+1);
    u1_ori_2(2*i) =  u1_ori_0(2*i);
    u1_ori_2(2*i+1) = -u1_ori_0(2*i+1);
    u1_ori_3(2*i) = - u1_ori_0(2*i);
    u1_ori_3(2*i+1) = - u1_ori_0(2*i+1);
  }
  u2_ori_0 = uu2.head(2*fN2);
  for( int i=0; i < fN2; ++i ){ 
    u2_ori_1(2*i) = -u2_ori_0(2*i);
    u2_ori_1(2*i+1) = u2_ori_0(2*i+1);
    u2_ori_2(2*i) =  u2_ori_0(2*i);
    u2_ori_2(2*i+1) = -u2_ori_0(2*i+1);
    u2_ori_3(2*i) = - u2_ori_0(2*i);
    u2_ori_3(2*i+1) = - u2_ori_0(2*i+1);
  }

  diff(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){    
    u1_ori_1(2*i) -= diff(0);
    u1_ori_1(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){    
    u2_ori_1(2*i) -= diff(0);
    u2_ori_1(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_2(2*fi[0]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_2(2*fi[0]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){  
    u1_ori_2(2*i) -= diff(0);
    u1_ori_2(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){  
    u2_ori_2(2*i) -= diff(0);
    u2_ori_2(2*i+1) -= diff(1);
  }
  
  diff(0) = u_ori_3(2*fi[3]) - u_ori_0(2*fi[4]); 
  diff(1) = u_ori_3(2*fi[3]+1) - u_ori_0(2*fi[4]+1); 
  for( int i = 0; i < n; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN1; ++i ){
    u1_ori_3(2*i) -= diff(0);
    u1_ori_3(2*i+1) -= diff(1);
  }
  for( int i = 0; i < fN2; ++i ){
    u2_ori_3(2*i) -= diff(0);
    u2_ori_3(2*i+1) -= diff(1);
  }
  

  diff_0(0) = u_ori_3(2*fi[0]) - u_ori_0(2*fi[5]);
  diff_0(1) = u_ori_3(2*fi[0]+1) - u_ori_0(2*fi[5]+1);
  diff_1(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[5]);
  diff_1(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[5]+1);

  double pointsD1[fN1][2];
  double pointsD2[fN2][2];

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD1[i][1] = u1_ori_0(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD2[i][1] = u2_ori_0(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      
      g_setColor(d_num_tiling,10);  // 10 
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19);  // 15
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_1(ix1(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD1[i][1] = u1_ori_1(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_1(ix2(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD2[i][1] = u2_ori_1(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      
      g_setColor(d_num_tiling,10);  // 19
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 29
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_2(ix1(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD1[i][1] = u1_ori_2(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_2(ix2(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD2[i][1] = u2_ori_2(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }

      g_setColor(d_num_tiling,10);  // 29
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19);  // 19
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1); 
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );

      // polygon
      for( int i = 0; i < fN1 ; ++i ){
	pointsD1[i][0] = u1_ori_3(ix1(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD1[i][1] = u1_ori_3(iy1(i))+diff_0(1)*k1+diff_1(1)*k2;
      }
      for( int i = 0; i < fN2 ; ++i ){
	pointsD2[i][0] = u2_ori_3(ix2(i))+diff_0(0)*k1+diff_1(0)*k2; 
	pointsD2[i][1] = u2_ori_3(iy2(i))+diff_0(1)*k1+diff_1(1)*k2;
      }

      g_setColor(d_num_tiling,10); // 15
      g_polygonfill( d_num_tiling, pointsD1, fN1 );
      g_setColor(d_num_tiling,19); // 19
      g_polygonfill( d_num_tiling, pointsD2, fN2 );
      g_setColor(d_num_tiling,1);
      g_polygon( d_num_tiling, pointsD1, fN1, 4 );
      g_polygon( d_num_tiling, pointsD2, fN2, 4 );
    }
  }
  
  g_flush();
}


void TDisplay::All_tiles( VectorXd* uu1_top, VectorXd* uu2_top )
{
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 2500, 2000, (char*)"Tiling" );
    g_window( d_num, 0, 0, 2500, 2000 );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_setFont(d_num,(char*)"-*-*-20-*-*-*-*-*-*-*"); // これはエラーでない
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);

  double L = 500.0;

  VectorXd uu1, uu2;
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  
  char str[80];
  double ax = 0.5*L;
  double  ay = 0.5*L; 
  for( int s = 0; s < 20; ++s ){

    if( s % 5 == 0 && s != 0 ){
      ax = 0.5*L;
      ay += L;
    }

    uu1 = uu1_top[s];
    uu2 = uu2_top[s];    
    x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9;
    
    for( int i = 0; i < fNN1; ++i ){
    if( uu1(2*i) < x_min )
      x_min = uu1(iix1(i));
    if( uu1(2*i) > x_max )
      x_max = uu1(iix1(i));
    if( uu1(2*i+1) < y_min )
      y_min = uu1(iiy1(i));
    if( uu1(2*i+1) > y_max )
      y_max = uu1(iiy1(i));
    }
    for( int i = 0; i < fNN2; ++i ){
      if( uu2(2*i) < x_min )
	x_min = uu2(iix2(i));
      if( uu2(2*i) > x_max )
	x_max = uu2(iix2(i));
      if( uu2(2*i+1) < y_min )
	y_min = uu2(iiy2(i));
      if( uu2(2*i+1) > y_max )
	y_max = uu2(iiy2(i));
    }
    if( x_max - x_min < y_max - y_min )
      length = y_max - y_min;
    else
      length = x_max - x_min;

    for( int i = 0; i < fNN1; ++i ){
      uu1(iix1(i)) -= 0.5*(x_min+x_max);
      uu1(iiy1(i)) -= 0.5*(y_min+y_max); 
    }
    for( int i = 0; i < fNN2; ++i ){
      uu2(iix2(i)) -= 0.5*(x_min+x_max);
      uu2(iiy2(i)) -= 0.5*(y_min+y_max); 
    }

    uu1 *= L/length*0.9;
    uu2 *= L/length*0.9;

    g_setColor(d_num,2); // red
    for( int i = 0; i < fN1; ++i )
      g_line( d_num, uu1(ix1(i))+ax, uu1(iy1(i))+ay, uu1(ix1(i+1))+ax, uu1(iy1(i+1))+ay, 2 );
    for( int i = 0; i < fN2; ++i )
      g_line( d_num, uu2(ix2(i))+ax, uu2(iy2(i))+ay, uu2(ix2(i+1))+ax, uu2(iy2(i+1))+ay, 2 );

    ax += L;
  }

  g_flush();  
}
  
