/*
Author: Yuichi Nagata
Copyright (c) 2023, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_ADD__
#include "operator_additional.h"
#endif


TOperator_ADD::TOperator_ADD( int N1, int N1_in, int N2, int N2_in ) : TOperator_IR( N1, N1_in, N2, N2_in )
{ 
}

TOperator_ADD::~TOperator_ADD()
{

}

void TOperator_ADD::SetParameter() 
{
  TOperator_IR::SetParameter();
}

void TOperator_ADD::SetInit( VectorXd w1,  VectorXd w1_in, MatrixXi kk1i, VectorXd w2,  VectorXd w2_in, MatrixXi kk2i )
{
  TOperator_IR::SetInit( w1, w1_in, kk1i, w2, w2_in, kk2i );

  // fNei1 = MatrixXi::Zero(fNN1,20); // TOperator_Iで定義
  // fNei1_Size = VectorXi::Zero(fNN1);
  fTriangleList1 = MatrixXi(10*fNN1,3); 
  fEdgeNei1_1 = MatrixXi(fNN1,fNN1);
  fEdgeNei2_1 = MatrixXi(fNN1,fNN1);
  fCell1_Edge = new MatrixXi [20*fNN1]; 
  for( int i = 0; i < 20*fNN1; ++i )
    fCell1_Edge[i] = MatrixXi(100,2);
  fCell1_Size = VectorXi::Zero(20*fNN1);

  // fNei2 = MatrixXi::Zero(fNN2,20);
  // fNei2_Size = VectorXi::Zero(fNN2);
  fTriangleList2 = MatrixXi(10*fNN2,3); 
  fEdgeNei1_2 = MatrixXi(fNN2,fNN2);
  fEdgeNei2_2 = MatrixXi(fNN2,fNN2);
  fCell2_Edge = new MatrixXi [20*fNN2]; 
  for( int i = 0; i < 20*fNN2; ++i )
    fCell2_Edge[i] = MatrixXi(100,2);
  fCell2_Size = VectorXi::Zero(20*fNN2);

  this->SetCell1_1_ling();
  this->SetNei1_Sub();

  this->SetCell2_1_ling();
  this->SetNei2_Sub();
}


void TOperator_ADD::SetCell1_1_ling()
{
  fNumOfCell1 = 0;
  for( int i = 0; i < fNN1; ++i ){
    for( int s = 0; s < fNei1_Size(i); ++s ){
      int j = fNei1(i,s); 
      fCell1_Edge[fNumOfCell1](s,0) = i;
      fCell1_Edge[fNumOfCell1](s,1) = j;
    }
    fCell1_Size(i) = fNei1_Size(i); 
    ++fNumOfCell1;
  }
}

void TOperator_ADD::SetCell2_1_ling()
{
  fNumOfCell2 = 0;
  for( int i = 0; i < fNN2; ++i ){
    for( int s = 0; s < fNei2_Size(i); ++s ){
      int j = fNei2(i,s); 
      fCell2_Edge[fNumOfCell2](s,0) = i;
      fCell2_Edge[fNumOfCell2](s,1) = j;
    }
    fCell2_Size(i) = fNei2_Size(i); 
    ++fNumOfCell2;
  }
}


void TOperator_ADD::SetNei1_Sub()
{
  //////// set edgeNei1 and edgeNei2 /////// 
  fEdgeNei1_1 = -1 * MatrixXi::Ones(fNN1,fNN1); // edgeを含むtriangleの残りの頂点
  fEdgeNei2_1 = -1 * MatrixXi::Ones(fNN1,fNN1);
  for( int i = 0; i < fNN1-1; ++i ){
    for( int j = 0; j < i+1; ++j ){
      if( fKK1i(i,j) != 0 ){
	int count = 0;
	for( int s1 = 0; s1 < fNei1_Size(i); ++s1 ){
	  int v1 = fNei1(i,s1);
	  for( int s2 = 0; s2 < fNei1_Size(j); ++s2 ){
	    int v2 = fNei1(j,s2);
	    if( v1 == v2 ){
	      if( count == 0 ){
		fEdgeNei1_1(i,j) = v1;
		fEdgeNei1_1(j,i) = v1;
		++count;
	      }
	      else if( count == 1 ){
		fEdgeNei2_1(i,j) = v1;
		fEdgeNei2_1(j,i) = v1;
		++count;
	      }
	      else
		assert( 1 == 2 );
	    }
	  }
	}
      }
    }
  }
  MatrixXi edgeNei1 = fEdgeNei1_1;
  MatrixXi edgeNei2 = fEdgeNei2_1;

  //////// set triangleList /////// 
  fTriangle1_num = 0;               
  int v1, v2, v3;
  for( v1 = 0; v1 < fNN1-1; ++v1 ){
    for( v2 = 0; v2 < v1+1; ++v2 ){
      v3 = edgeNei1(v1,v2);
      edgeNei1(v1,v2) = -1; 
      edgeNei1(v2,v1) = -1; 
      if( v3 != -1 ){
	fTriangleList1(fTriangle1_num,0) = v1;
	fTriangleList1(fTriangle1_num,1) = v2;
	fTriangleList1(fTriangle1_num,2) = v3;
	// printf( "%d: (%d, %d, %d)\n", fTriangle_num1, v1, v2, v3 );
	++fTriangle1_num;

	if( edgeNei1(v1,v3) == v2 ){
	  edgeNei1(v1,v3) = -1;
	  edgeNei1(v3,v1) = -1;
	}
	if( edgeNei2(v1,v3) == v2 ){
	  edgeNei2(v1,v3) = -1;
	  edgeNei2(v3,v1) = -1;
	}

	if( edgeNei1(v2,v3) == v1 ){
	  edgeNei1(v2,v3) = -1;
	  edgeNei1(v3,v2) = -1;
	}
	if( edgeNei2(v2,v3) == v1 ){
	  edgeNei2(v2,v3) = -1;
	  edgeNei2(v3,v2) = -1;
	}
      }

      v3 = edgeNei2(v1,v2);
      edgeNei2(v1,v2) = -1; 
      edgeNei2(v2,v1) = -1; 
      if( v3 != -1 ){
	fTriangleList1(fTriangle1_num,0) = v1;
	fTriangleList1(fTriangle1_num,1) = v2;
	fTriangleList1(fTriangle1_num,2) = v3;
	++fTriangle1_num;

	if( edgeNei1(v1,v3) == v2 ){
	  edgeNei1(v1,v3) = -1;
	  edgeNei1(v3,v1) = -1;
	}
	if( edgeNei2(v1,v3) == v2 ){
	  edgeNei2(v1,v3) = -1;
	  edgeNei2(v3,v1) = -1;
	}

	if( edgeNei1(v2,v3) == v1 ){
	  edgeNei1(v2,v3) = -1;
	  edgeNei1(v3,v2) = -1;
	}
	if( edgeNei2(v2,v3) == v1 ){
	  edgeNei2(v2,v3) = -1;
	  edgeNei2(v3,v2) = -1;
	}
      }
    }
  }
}

void TOperator_ADD::SetNei2_Sub()
{
  //////// set edgeNei1 and edgeNei2 /////// 
  fEdgeNei1_2 = -1 * MatrixXi::Ones(fNN2,fNN2); // edgeを含むtriangleの残りの頂点
  fEdgeNei2_2 = -1 * MatrixXi::Ones(fNN2,fNN2);
  for( int i = 0; i < fNN2-1; ++i ){
    for( int j = 0; j < i+1; ++j ){
      if( fKK2i(i,j) != 0 ){
	int count = 0;
	for( int s1 = 0; s1 < fNei2_Size(i); ++s1 ){
	  int v1 = fNei2(i,s1);
	  for( int s2 = 0; s2 < fNei2_Size(j); ++s2 ){
	    int v2 = fNei2(j,s2);
	    if( v1 == v2 ){
	      if( count == 0 ){
		fEdgeNei1_2(i,j) = v1;
		fEdgeNei1_2(j,i) = v1;
		++count;
	      }
	      else if( count == 1 ){
		fEdgeNei2_2(i,j) = v1;
		fEdgeNei2_2(j,i) = v1;
		++count;
	      }
	      else 
		assert( 1 == 2 );
	    }
	  }
	}
      }
    }
  }
  MatrixXi edgeNei1 = fEdgeNei1_2;
  MatrixXi edgeNei2 = fEdgeNei2_2;

  //////// set triangleList /////// 
  fTriangle2_num = 0;               
  int v1, v2, v3;
  for( v1 = 0; v1 < fNN2-1; ++v1 ){
    for( v2 = 0; v2 < v1+1; ++v2 ){
      v3 = edgeNei1(v1,v2);
      edgeNei1(v1,v2) = -1; 
      edgeNei1(v2,v1) = -1; 
      if( v3 != -1 ){
	fTriangleList2(fTriangle2_num,0) = v1;
	fTriangleList2(fTriangle2_num,1) = v2;
	fTriangleList2(fTriangle2_num,2) = v3;
	// printf( "%d: (%d, %d, %d)\n", fTriangle_num2, v1, v2, v3 );
	++fTriangle2_num;

	if( edgeNei1(v1,v3) == v2 ){
	  edgeNei1(v1,v3) = -1;
	  edgeNei1(v3,v1) = -1;
	}
	if( edgeNei2(v1,v3) == v2 ){
	  edgeNei2(v1,v3) = -1;
	  edgeNei2(v3,v1) = -1;
	}

	if( edgeNei1(v2,v3) == v1 ){
	  edgeNei1(v2,v3) = -1;
	  edgeNei1(v3,v2) = -1;
	}
	if( edgeNei2(v2,v3) == v1 ){
	  edgeNei2(v2,v3) = -1;
	  edgeNei2(v3,v2) = -1;
	}
      }

      v3 = edgeNei2(v1,v2);
      edgeNei2(v1,v2) = -1; 
      edgeNei2(v2,v1) = -1; 
      if( v3 != -1 ){
	fTriangleList2(fTriangle2_num,0) = v1;
	fTriangleList2(fTriangle2_num,1) = v2;
	fTriangleList2(fTriangle2_num,2) = v3;
	++fTriangle2_num;

	if( edgeNei1(v1,v3) == v2 ){
	  edgeNei1(v1,v3) = -1;
	  edgeNei1(v3,v1) = -1;
	}
	if( edgeNei2(v1,v3) == v2 ){
	  edgeNei2(v1,v3) = -1;
	  edgeNei2(v3,v1) = -1;
	}

	if( edgeNei1(v2,v3) == v1 ){
	  edgeNei1(v2,v3) = -1;
	  edgeNei1(v3,v2) = -1;
	}
	if( edgeNei2(v2,v3) == v1 ){
	  edgeNei2(v2,v3) = -1;
	  edgeNei2(v3,v2) = -1;
	}
      }
    }
  } 
}


void TOperator_ADD::Find_Triangle1( double xp, double yp, VectorXd& uup, int& index_triangle, double& alpha, double& beta )
{
  int v0, v1, v2;
  double u1, u2, a1, a2, b1, b2;

  index_triangle = -1;
  for( int t = 0; t < fTriangle1_num; ++t ){
    v0 = fTriangleList1(t,0);
    v1 = fTriangleList1(t,1);
    v2 = fTriangleList1(t,2);

    u1 = xp - uup(iix1(v0));
    u2 = yp - uup(iiy1(v0));
  

    a1 = uup(iix1(v1)) - uup(iix1(v0));
    a2 = uup(iiy1(v1)) - uup(iiy1(v0));
    b1 = uup(iix1(v2)) - uup(iix1(v0));
    b2 = uup(iiy1(v2)) - uup(iiy1(v0));
    alpha = (b2*u1-b1*u2)/(a1*b2- a2*b1);
    beta = (a2*u1-a1*u2)/(b1*a2- b2*a1);

    if( 0 <= alpha && 0 <= beta && alpha + beta <= 1.0 ){
      index_triangle = t;
      break;
    }
  }
}

void TOperator_ADD::Find_Triangle2( double xp, double yp, VectorXd& uup, int& index_triangle, double& alpha, double& beta )
{
  int v0, v1, v2;
  double u1, u2, a1, a2, b1, b2;

  index_triangle = -1;
  for( int t = 0; t < fTriangle2_num; ++t ){
    v0 = fTriangleList2(t,0);
    v1 = fTriangleList2(t,1);
    v2 = fTriangleList2(t,2);

    u1 = xp - uup(iix2(v0));
    u2 = yp - uup(iiy2(v0));
  

    a1 = uup(iix2(v1)) - uup(iix2(v0));
    a2 = uup(iiy2(v1)) - uup(iiy2(v0));
    b1 = uup(iix2(v2)) - uup(iix2(v0));
    b2 = uup(iiy2(v2)) - uup(iiy2(v0));
    alpha = (b2*u1-b1*u2)/(a1*b2- a2*b1);
    beta = (a2*u1-a1*u2)/(b1*a2- b2*a1);

    if( 0 <= alpha && 0 <= beta && alpha + beta <= 1.0 ){
      index_triangle = t;
      break;
    }
  }
}


bool TOperator_ADD::Is_inside_Up1( double x, double y, VectorXd& uup )
{
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;

  x1 = x;
  y1 = y;
  x2 = x1+10000.0;
  y2 = y1+10.0;

  
  int count = 0;
  for( int i = 0; i < fN1; ++i ){
    x3 = uup(ix1(i));
    y3 = uup(iy1(i));
    x4 = uup(ix1(i+1));
    y4 = uup(iy1(i+1));

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
    return true; // ポリゴン内部 -> true
  else 
    return false; // ポリゴン外部 -> false
}

bool TOperator_ADD::Is_inside_Up2( double x, double y, VectorXd& uup )
{
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;

  x1 = x;
  y1 = y;
  x2 = x1+10000.0;
  y2 = y1+10.0;

  
  int count = 0;
  for( int i = 0; i < fN2; ++i ){
    x3 = uup(ix2(i));
    y3 = uup(iy2(i));
    x4 = uup(ix2(i+1));
    y4 = uup(iy2(i+1));

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
    return true; // ポリゴン内部 -> true
  else 
    return false; // ポリゴン外部 -> false
}


double TOperator_ADD::Cal_Area1( VectorXd& uu )
{
  int v0, v1, v2;
  double x0, x1, x2, y0, y1, y2;
  double S;
  

  S = 0.0;
  for( int t = 0; t < fTriangle1_num; ++t ){
    v0 = fTriangleList1(t,0);
    v1 = fTriangleList1(t,1);
    v2 = fTriangleList1(t,2);
    x0 = uu(iix1(v0));
    y0 = uu(iiy1(v0));
    x1 = uu(iix1(v1));
    y1 = uu(iiy1(v1));
    x2 = uu(iix1(v2));
    y2 = uu(iiy1(v2));
    S += 0.5*fabs((x0-x2)*(y1-y2)-(x1-x2)*(y0-y2));
  }

  return S;
}

double TOperator_ADD::Cal_Area2( VectorXd& uu )
{
  int v0, v1, v2;
  double x0, x1, x2, y0, y1, y2;
  double S;
  

  S = 0.0;
  for( int t = 0; t < fTriangle2_num; ++t ){
    v0 = fTriangleList2(t,0);
    v1 = fTriangleList2(t,1);
    v2 = fTriangleList2(t,2);
    x0 = uu(iix2(v0));
    y0 = uu(iiy2(v0));
    x1 = uu(iix2(v1));
    y1 = uu(iiy2(v1));
    x2 = uu(iix2(v2));
    y2 = uu(iiy2(v2));
    S += 0.5*fabs((x0-x2)*(y1-y2)-(x1-x2)*(y0-y2));
  }

  return S;
}
