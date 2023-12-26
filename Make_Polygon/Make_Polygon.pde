/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

// Instructions of this program

// Some parameter values must be set here
int flag_image = 0; // 0 (default); If set to 1, the specified image (png,jpeg) will be displayed.
double scale = 1.0; // Scaling of the input image
String st_image="example1.png";  // File name of an input image

int flag_read_data = 0; // 0 (default); If set to 1,the specified polygon data is loaded. 
String st_data="example1_75.dat"; // File name of a polygon data

// There are three modes (0, 1, and 2) in the process of creating polygons. Mode 1 is set at the beginning. 
// The mode can be changed by activating the displayed window and entering a number (0, 1, or 2).
// mode 0: By clicking on a point, it can be erased.
// mode 1: By clicking on the window, a point at the clicked location is inserted between the first last points.
// mode 2: By clicking on the window, a point at the clicked location is inserted between two consecutive points close to it.
// mode 2: By dragging a displayed point, you can move it. 

// To outpt the current data, enter "q". The result is output to the file output.dat.

// Edges of the polygon displayed in the window are colored in three different ways:
// -black: the length of this edge is approximately equal to the average length of the edges of the polygon.  
// -red: the length of this edge is longer than the average.  
// -green: the length of this edge is smaller than the average. 
// Not all edges must be black.
// To change the threshold, search for "threshold" in this program. 


///////////// End of the setting ////////////


PImage img;
int N = 300;
int [][] xy = new int[N][2];
int mode; // 0:delete, 1:add to last, 2:add to nearest
int n;
int pointMove;

int n_ori;
int [][] xy_ori = new int[N][2];



void setup() {
  size(1500, 1000);
  if ( flag_image == 1 ) {
    img = loadImage(st_image); // 読み込み画像
    img.resize( (int)(img.width*scale), (int)(img.height*scale) );
    image(img, 0, 0); 
  }

  n =  0;
  if( flag_read_data == 1){
    readData(); // 既存データを読み込む
    drawPoints();
  }
  mode = 1;
}

void readData()
{
  String lines[] = null;
  String lin = null;
  
  lines = loadStrings(st_data); // 読み込みデータ
  n = 0;
  if ( lines != null ) {
    String[] co = split(lines[0], ' ');
    if (co.length==1)
      n = int(co[0]);
    for ( int i = 0; i < n; ++i )
    {
      co = split(lines[i+1], ' ');
      if (co.length==2) {
        xy[i][0] = int(co[0]);
        xy[i][1] = int(co[1]);
      }
    }
  }

  n_ori = n;
  for ( int i = 0; i < n; ++i )
  {
    xy_ori[i][0] = xy[i][0];
    xy_ori[i][1] = xy[i][1];
  }
  drawPoints_ori();
   
}


void writeData()
{
  PrintWriter outfile;
  outfile = createWriter("output.dat");
  // outfile = createWriter(C:\Users\Nagata\Documents\pos_after.txt);
  outfile.println(n);
  for ( int i = 0; i < n; ++i )
    outfile.println((double)xy[i][0] + " " + (double)xy[i][1]);
  outfile.flush(); //残りを出力する
  outfile.close(); // ファイルを閉じる
}

void draw() {
  // println( mouseX, " ", mouseY );
}


void mouseDragged()
{
  int pointPush;

  if ( pointMove != -1 )
    movePoint( pointMove );
  drawPoints();
}


void mousePressed()
{
  int pointPush;

  pointPush = pointPushed();

  pointMove = -1;
  if ( mode == 1 ) {
    addPointLast();
  } else if ( mode == 2 ) {
    if ( pointPush == -1 ) {
      addPointNearest();
    } else {
      pointMove = pointPush;
    }
  } else if ( mode == 0 ) {
    if ( pointPush != -1 )
      deletePoint( pointPush );
  }

  drawPoints();

  println("n=", n, " ", mouseX, " ", mouseY);
}

int pointPushed()
{
  for ( int i = 0; i < n; ++i ) {
    if ( xy[i][0]-5 < mouseX && mouseX < xy[i][0]+5 &&
      xy[i][1]-5 < mouseY && mouseY < xy[i][1]+5 )
    {
      return i;
    }
  }
  return -1;
}

void addPointLast()
{
  xy[n][0] = mouseX;
  xy[n][1] = mouseY;
  ++n;
}

void insertPoint( int pointIn )
{
  for ( int i = n; i >= pointIn+1; --i ) {
    xy[i][0] = xy[i-1][0];
    xy[i][1] = xy[i-1][1];
  }
  xy[pointIn][0] = mouseX;
  xy[pointIn][1] = mouseY;
  ++n;
}

void addPointNearest()
{
  double d;
  double d_min;
  int pointNear, pointP, pointN;
  pointNear = -1;
  d_min = 9999999.9;
  for ( int i = 0; i < n; ++i ) {
    d = (xy[i][0] - mouseX )*(xy[i][0] - mouseX )+(xy[i][1] - mouseY )*(xy[i][1] - mouseY );
    if ( d < d_min ) {
      d_min = d;
      pointNear = i;
    }
  }
  double d1, d2;
  int pointIn;

  if ( pointNear != 0 )
    pointP = pointNear-1;
  else
  pointP = n-1;
  pointN = ( pointNear + 1 ) % n;
  if ( (xy[pointP][0] - mouseX )*(xy[pointP][0] - mouseX )+(xy[pointP][1] - mouseY )*(xy[pointP][1] - mouseY ) <
    (xy[pointN][0] - mouseX )*(xy[pointN][0] - mouseX )+(xy[pointN][1] - mouseY )*(xy[pointN][1] - mouseY ) )
    pointIn = pointNear;
  else
  pointIn = pointNear+1;

  insertPoint( pointIn );
}


void movePoint( int point )
{
  xy[point][0] = mouseX;
  xy[point][1] = mouseY;
}

void deletePoint( int pointDelete )
{
  for ( int i = pointDelete; i < n-1; ++i ) {
    xy[i][0] = xy[i+1][0];
    xy[i][1] = xy[i+1][1];
  }
  --n;
}

void drawPoints()
{
  background(200);
  if ( flag_image == 1 )
    image(img, 0, 0);
  drawPoints_ori();
  rectMode(CENTER);
  stroke(255, 0, 0);
  strokeWeight(2);
  for ( int i = 0; i < n; ++i ) {
    rect(xy[i][0], xy[i][1], 8, 8);
  }

  float len, dx, dy;
  len = 0.0;
  for ( int i = 0; i < n-1; ++i ) {
    dx = xy[i+1][0]- xy[i][0];
    dy = xy[i+1][1]- xy[i][1];
    len += sqrt(dx*dx+dy*dy);
  }
  dx = xy[n-1][0]- xy[0][0];
  dy = xy[n-1][1]- xy[0][1];
  len += sqrt(dx*dx+dy*dy);
  len /=(float)n;

  double l;
  strokeWeight(1);
  stroke(0, 0, 0);
  for ( int i = 0; i < n; ++i ) {
    if ( i != n-1 ) {
      dx = xy[i+1][0]- xy[i][0];
      dy = xy[i+1][1]- xy[i][1];
    } else {
      dx = xy[n-1][0]- xy[0][0];
      dy = xy[n-1][1]- xy[0][1];
    }
    l = sqrt(dx*dx+dy*dy);
   // print(len); print(" "); print(l);print("\n");
    if ( l < 0.95*len){ // threshold
      stroke(0, 255, 0);
      strokeWeight(2);
    }
    else if ( l > 1.05*len){ // threshold
      stroke(255, 0, 0);
      strokeWeight(2);
    }
    else{
      stroke(0, 0, 0);
      strokeWeight(1);
    }

    if ( i != n-1 )
      line(xy[i][0], xy[i][1], xy[i+1][0], xy[i+1][1] );
    else
    line(xy[n-1][0], xy[n-1][1], xy[0][0], xy[0][1] );
  }
}

void drawPoints_ori()
{
  rectMode(CENTER);
  stroke(255, 255, 255);
  for ( int i = 0; i < n_ori; ++i ) {
    rect(xy_ori[i][0], xy_ori[i][1], 8, 8);
  }
}

void keyPressed() {//（操作箇所）
  if ( key == 'q' ) {
    writeData();
  }
  if ( key == '0' )
    mode = 0;
  if ( key == '1' )
    mode = 1;
  if ( key == '2' )
    mode = 2;
}
