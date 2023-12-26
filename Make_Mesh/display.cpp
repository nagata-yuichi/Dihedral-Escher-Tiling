/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#include "display.h"
#endif

int open_flag_tiling = 0;

TDisplay::TDisplay()
{
  fBoxL = 4;
}

TDisplay::~TDisplay()
{

}

void TDisplay::Disp_bit( VectorXd& u , VectorXd& w, int* fi, int nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
  int n = u.rows()/2;

  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < w.rows()/2; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.6*length; 
  xe = (x_min+x_max)*0.5 + 0.6*length; 
  ys = (y_min+y_max)*0.5 - 0.6*length; 
  ye = (y_min+y_max)*0.5 + 0.6*length; 

  d_num=1;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 600, 600, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xe, ye );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000); 
  }

  d_num=1;
  g_clearWindow(d_num);

  // 1:black, 2:red, 3:blue, 4:green, 9:blown, 6:pink
  cn[1] = 2; cn[2] = 2; cn[3] = 2; cn[4] = 2; cn[5] = 2; cn[6] = 2;    // all red


  // shape line (color)
  for( int h = 1; h < nv; ++h ){
    g_setColor(d_num,cn[h]); 
    for( int i = fi[h-1]; i < fi[h]; ++i )
      g_line( d_num, u(2*i), u(2*i+1), u(2*(i+1)), u(2*(i+1)+1), 4 );
  }
  g_setColor(d_num,cn[nv]); 
  for( int i = fi[nv-1]; i < fi[nv]-1; ++i )
    g_line( d_num, u(2*i), u(2*i+1), u(2*(i+1)), u(2*(i+1)+1), 4 );
  g_line( d_num, u(2*(n-1)), u(2*(n-1)+1), u(0), u(1), 4 );






  // points of the goal
    g_setColor(d_num,1); // black
  for( int i = 0; i < n; ++i ){
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
  }
  g_flush();
}



void TDisplay::Disp_AdjacentRelation( int n, VectorXd w, MatrixXi K ) 
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
 
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.6*length; 
  xe = (x_min+x_max)*0.5 + 0.6*length; 
  ys = (y_min+y_max)*0.5 - 0.6*length; 
  ye = (y_min+y_max)*0.5 + 0.6*length; 

  d_num=2;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 600, 600, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xe, ye );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }

  // goal


  int nn = K.rows();
  //  printf("nn= %d\n", nn ); 
  //  cout << K << endl;

  g_setColor(d_num,1); // black
  for( int i = 0; i < nn-1; ++i ){
    for( int j = i+1; j < nn; ++j ){
      if( K(i,j) <= -1 ){
	g_line( d_num, w(2*i), w(2*i+1), w(2*j), w(2*j+1), 2 );
      }
    }
  }

  // goal
  for( int i = 0; i < n; ++i ){
    g_setColor(d_num,2); // red
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
    g_setColor(d_num,1); // red
    g_line( d_num, w(2*i), w(2*i+1), w(2*((i+1)%n)), w(2*((i+1)%n)+1), 2 );
  }

  // inner point
  for( int i = n; i < nn; ++i ){
    g_setColor(d_num,4); // green
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
  }

  g_flush();
}

void TDisplay::Disp_AdjacentRelation2( int n, VectorXd w, MatrixXi K ) 
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
 
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.6*length; 
  xe = (x_min+x_max)*0.5 + 0.6*length; 
  ys = (y_min+y_max)*0.5 - 0.6*length; 
  ye = (y_min+y_max)*0.5 + 0.6*length; 

  d_num=3;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 600, 600, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xe, ye );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }
  g_clearWindow(d_num);
  


  int nn = K.rows();
  //  printf("nn= %d\n", nn ); 
  //  cout << K << endl;

  g_setColor(d_num,1); // black
  for( int i = 0; i < nn-1; ++i ){
    for( int j = i+1; j < nn; ++j ){
      if( K(i,j) <= -1 ){
	g_line( d_num, w(2*i), w(2*i+1), w(2*j), w(2*j+1), 2 );
      }
    }
  }

  // goal
  for( int i = 0; i < n; ++i ){
    g_setColor(d_num,2); // red
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
    g_setColor(d_num,1); // red
    g_line( d_num, w(2*i), w(2*i+1), w(2*((i+1)%n)), w(2*((i+1)%n)+1), 2 );
  }

  // inner point
  for( int i = n; i < nn; ++i ){
    g_setColor(d_num,4); // green
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
  }

  g_flush();
}



void TDisplay::Disp_goal( VectorXd w, MatrixXd G )  
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
  int n = w.rows()/2;


  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.6*length; 
  xe = (x_min+x_max)*0.5 + 0.6*length; 
  ys = (y_min+y_max)*0.5 - 0.6*length; 
  ye = (y_min+y_max)*0.5 + 0.6*length; 

  d_num=1;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 600, 600, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xe, ye );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }

  d_num=1;
  g_clearWindow(d_num);

  // points of the goal
  int stock = fBoxL;
  fBoxL = 6;
  g_setColor(d_num,1); // black
  for( int i = 0; i < n; ++i ){
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
    if( i != n-1 )
      g_line( d_num, w(2*i), w(2*i+1), w(2*(i+1)), w(2*(i+1)+1), 2 );
    else 
      g_line( d_num, w(2*i), w(2*i+1), w(0), w(1), 2 );
  }

  ////// weighted partに色
  for( int i = 0; i < n-1; ++i ){
    for( int j = i+1; j < n; ++j ){
      if( G(i,j) < -0.000001 ){
	if( G(i,j) < -0.500001 )
	  g_setColor(d_num,4); // green
	else {
	  g_setColor(d_num,1); // black

	  if( j-i == 1 || j-i == n-1 ) 
	    g_setColor(d_num,2); // red
	  else
	    g_setColor(d_num,1); // black

	}
	// g_line( d_num, w(2*i), w(2*i+1), w(2*j), w(2*j+1), 2 );
      }
    }
  }

  for( int i = 0; i < n-1; ++i ){
    for( int j = i+1; j < n; ++j ){
      G(i,i) += G(i,j);
      G(j,j) += G(j,i);
    }
  }

  for( int i = 0; i < n; ++i ){
    if( G(i,i) > 1.000001 ){ 
      g_setColor(d_num,4); // green
      g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 3 );
    }
  }
  //////


  g_flush();
  fBoxL = stock;
}

void TDisplay::Disp_tile( VectorXd& w, VectorXd& u )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
  int n = w.rows()/2;


  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < n; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.6*length; 
  xe = (x_min+x_max)*0.5 + 0.6*length; 
  ys = (y_min+y_max)*0.5 - 0.6*length; 
  ye = (y_min+y_max)*0.5 + 0.6*length; 

  d_num=3;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 600, 600, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xe, ye );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }
    g_window( d_num, xs, ys, xe, ye );

  d_num=3;
  g_clearWindow(d_num);

  g_setColor(d_num,2); 
  for( int i = 0; i < n-1; ++i ){
    g_line( d_num, u(2*i), u(2*i+1), u(2*(i+1)), u(2*(i+1)+1), 4 );
  }
  g_line( d_num, u(0), u(1), u(2*(n-1)), u(2*(n-1)+1), 4 );

   // points of the goal
  g_setColor(d_num,1); // black
  for( int i = 0; i < n; ++i ){
    g_box( d_num, w(2*i)-fBoxL, w(2*i+1)-fBoxL, w(2*i)+fBoxL, w(2*i+1)+fBoxL, 2 );
  }

  g_flush();
}


// Size:4*5
void TDisplay::Disp_all( VectorXd* u , VectorXd& w, int** fi, int* nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
  int n = u[0].rows()/2;

  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < w.rows()/2; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.55*length; 
  xe = (x_min+x_max)*0.5 + 0.55*length; 
  ys = (y_min+y_max)*0.5 - 0.55*length; 
  ye = (y_min+y_max)*0.5 + 0.55*length; 

  d_num=1;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 1750, 1400, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xs+5.0*length*1.1, ys+4.0*length*1.1 );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }

  d_num=1;
  g_clearWindow(d_num);


  // 1:black, 2:red, 3:blue, 4:green, 9:blown, 6:pink
  // cn[1] = 2; cn[2] = 4; cn[3] = 3; cn[4] = 2; cn[5] = 26; cn[6] = 24;  // IH4
  // cn[1] = 2; cn[2] = 4; cn[3] = 4; cn[4] = 2; cn[5] = 3; cn[6] = 26;  // IH5
  // cn[1] = 2; cn[2] = 4; cn[3] = 2; cn[4] = 3; cn[5] = 4; cn[6] = 24;  // IH6
  cn[1] = 2; cn[2] = 2; cn[3] = 2; cn[4] = 2; cn[5] = 2; cn[6] = 2;    // all red
  // cn[1] = 2; cn[2] = 2; cn[3] = 4; cn[4] = 4; cn[5] = 3; cn[6] = 2;    // all red


  char str[80];
  double ax = 0.0;
  double  ay = 3.0*(ye-ys);
  for( int s = 0; s < 20; ++s ){
    if( s % 5 == 0 && s != 0 ){
      ax = 0.0;
      ay -= length*1.1;
    }

    

    // shape line (color)
    for( int h = 1; h < nv[s]; ++h ){
      g_setColor(d_num,cn[h]); 
      for( int i = fi[s][h-1]; i < fi[s][h]; ++i )
	g_line( d_num, u[s](2*i)+ax, u[s](2*i+1)+ay, u[s](2*(i+1))+ax, u[s](2*(i+1)+1)+ay, 3 );
    }
    g_setColor(d_num,cn[nv[s]]); 
    for( int i = fi[s][nv[s]-1]; i < fi[s][nv[s]]-1; ++i )
      g_line( d_num, u[s](2*i)+ax, u[s](2*i+1)+ay, u[s](2*(i+1))+ax, u[s](2*(i+1)+1)+ay, 3 );
    g_line( d_num, u[s](2*(n-1))+ax, u[s](2*(n-1)+1)+ay, u[s](0)+ax, u[s](1)+ay, 3 );

    // points of the goal
    g_setColor(d_num,1); // black
    for( int i = 0; i < n; ++i ){
      g_box( d_num, w(2*i)-fBoxL+ax, w(2*i+1)-fBoxL+ay, w(2*i)+fBoxL+ax, w(2*i+1)+fBoxL+ay, 3 );
    }

    // number
    sprintf( str, "%d", s );
    g_string( d_num, xs+ax+100, ys+ay+100, str );
    ax += length*1.1;
  }
  g_flush();
}

// For supplemental file

void TDisplay::Disp_all_sup( VectorXd* u , VectorXd& w, int** fi, int* nv )
{
  static int open_flg = 0;
  static int d_num = 0;
  int ai,aj;
  int color_num;
  int cn[7];
  int n = u[0].rows()/2;

  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < w.rows()/2; ++i ){
    if( w(2*i) < x_min )
      x_min = w(2*i);
    if( w(2*i) > x_max )
      x_max = w(2*i);
    if( w(2*i+1) < y_min )
      y_min = w(2*i+1);
    if( w(2*i+1) > y_max )
      y_max = w(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.55*length; 
  xe = (x_min+x_max)*0.5 + 0.55*length; 
  ys = (y_min+y_max)*0.5 - 0.55*length; 
  ye = (y_min+y_max)*0.5 + 0.55*length; 

  d_num=1;
  /* windowを開く------------------- */
  if (open_flg == 0) 
  {
    open_flg = 1;
    /*---ウインドウの名前を付ける--------------- */
    g_open((char*)"");
    /*---ウインドウの番号と大きさを指定----------*/
    g_open_window( d_num, 1400, 1750, (char*)"Tiling" );
    /*---ウインドウ中の座標系を定義する----------*/
    g_window( d_num, xs, ys, xs+4.0*length*1.1, ys+5.0*length*1.1 );
    // g_window( d_num, 0, 1000, 1000, 0 );
    /*---ウインドウ中で使用する文字フォントを定義*/
    //g_setFont(d_num,(char*)"-*-*-*-*-*-*-*-*-*-*-*-*-*-*");
    g_setFont(d_num,(char*)"-*-*-*-r-*-*-20-*-*-*-*-*-*-*");
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    /* g_set_pixpat(d_num); */
    usleep(200000);
  }

  d_num=1;
  g_clearWindow(d_num);


  // 1:black, 2:red, 3:blue, 4:green, 9:blown, 6:pink
  // cn[1] = 2; cn[2] = 4; cn[3] = 3; cn[4] = 2; cn[5] = 26; cn[6] = 24;  // IH4
  // cn[1] = 2; cn[2] = 4; cn[3] = 4; cn[4] = 2; cn[5] = 3; cn[6] = 26;  // IH5
  // cn[1] = 2; cn[2] = 4; cn[3] = 2; cn[4] = 3; cn[5] = 4; cn[6] = 24;  // IH6
  cn[1] = 2; cn[2] = 2; cn[3] = 2; cn[4] = 2; cn[5] = 2; cn[6] = 2;    // all red
  // cn[1] = 2; cn[2] = 2; cn[3] = 4; cn[4] = 4; cn[5] = 3; cn[6] = 2;    // all red


  char str[80];
  double ax = 0.0;
  double  ay = 5.0*(ye-ys);

  /*
  // points of the goal
  g_setColor(d_num,1); // black
  for( int i = 0; i < n; ++i ){
    g_box( d_num, w(2*i)-fBoxL+ax, w(2*i+1)-fBoxL+ay, w(2*i)+fBoxL+ax, w(2*i+1)+fBoxL+ay, 3 );
   if( i != n-1 )
      g_line( d_num, w(2*i)+ax, w(2*i+1)+ay, w(2*(i+1))+ax, w(2*(i+1)+1)+ay, 2 );
    else 
      g_line( d_num, w(2*i)+ax, w(2*i+1)+ay, w(0)+ax, w(1)+ay, 2 );
  }
  sprintf( str, "%s", "W" );
  g_string( d_num, xs+ax+50, ys+ay+50, str );
  */

  ay -= length*1.1;
  for( int s = 0; s < 20; ++s ){
    if( s % 4 == 0 && s != 0 ){
      ax = 0.0;
      ay -= length*1.1;
    }

    // shape line (color)
    for( int h = 1; h < nv[s]; ++h ){
      g_setColor(d_num,cn[h]); 
      for( int i = fi[s][h-1]; i < fi[s][h]; ++i )
	g_line( d_num, u[s](2*i)+ax, u[s](2*i+1)+ay, u[s](2*(i+1))+ax, u[s](2*(i+1)+1)+ay, 3 );
    }
    g_setColor(d_num,cn[nv[s]]); 
    for( int i = fi[s][nv[s]-1]; i < fi[s][nv[s]]-1; ++i )
      g_line( d_num, u[s](2*i)+ax, u[s](2*i+1)+ay, u[s](2*(i+1))+ax, u[s](2*(i+1)+1)+ay, 3 );
    g_line( d_num, u[s](2*(n-1))+ax, u[s](2*(n-1)+1)+ay, u[s](0)+ax, u[s](1)+ay, 3 );

    // points of the goal
    /*
    g_setColor(d_num,1); // black
    for( int i = 0; i < n; ++i ){
      g_box( d_num, w(2*i)-fBoxL+ax, w(2*i+1)-fBoxL+ay, w(2*i)+fBoxL+ax, w(2*i+1)+fBoxL+ay, 3 );
    }
    */

    // number
    g_setColor(d_num,1); // black
    sprintf( str, "%d", s+1 );
    g_string( d_num, xs+ax+50, ys+ay+50, str );
    ax += length*1.1;
  }
  g_flush();
}

