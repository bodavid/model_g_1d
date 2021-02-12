#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
// BJD added <fstream>
#include <fstream>
//#include <chplot.h>
using namespace std;

typedef double real;

const real dx = 1, dy = 12, a = 14, b = 29,
  g = 0.1, p = 1, q = 1, s = 0, u = 0, w = 0;

const real G0 = (a+g*w)/(q-g*p);
const real X0 = (p*a+q*w)/(q-g*p);
const real Y0 = (s*X0*X0+b)/(X0*X0+u)*X0;

//const real writeFile();

const int Length = 100; // Size of space in system units
const int Time = 100; // Duration of time in system units
const int Nx = 8; // Number of computational cells per unit length.
const int Nt = 180; // Number of computational cells per unit time.
const int Mx = Nx*Length; // Size of space in computational units
const int Mt = Nt*Time; // Duration of time in computational units

int count = 0;

enum { G, X, Y };
const real vacuum[] = { G0, X0, Y0 };
const real pcn = real(Nx*Nx)/(2*Nt); // pre-Crank-Nicolson coefficient
const real d[] = { pcn, dx*pcn, dy*pcn }; // Crank-Nicolson coefficients
// 2: current and next iteration
// 3: phi_G, phi_X, phi_Y
// Mx+2: Mx computational cells in entire spatial dimension, plus 0 on each end
real space[2][3][Mx+2];

//ofstream myfile;
//myfile.open ("/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile.txt");

void init()
{
  for( int i=0 ; i<2 ; ++i )
    for( int j=0 ; j<3 ; ++j )
      for( int k=0 ; k<Mx+2 ; ++k ) space[i][j][k] = 0;
}

real chi( const real &x, const real &t )
{
  real xoff = x - real(Length)/2;
  real toff = t - 10;
  return -exp(-xoff*xoff/2)*exp(-toff*toff/18);
}

real reaction( const real (*pcurr)[Mx+2], int R, int j, int t )
{
  switch(R)
  {
    case G: return -q*pcurr[G][j] + g*pcurr[X][j];
    case X: return p*pcurr[G][j]
                   - (1+b+3*s*X0*X0-2*X0*Y0)*pcurr[X][j]
                   + (u+X0*X0)*pcurr[Y][j]
                   + (Y0-3*s*X0)*pcurr[X][j]*pcurr[X][j]
                   + 2*X0*pcurr[X][j]*pcurr[Y][j]
                   - s*pcurr[X][j]*pcurr[X][j]*pcurr[X][j]
                   + pcurr[X][j]*pcurr[X][j]*pcurr[Y][j]
                   + chi(real(j)/Nx,real(t)/Nt);
    case Y: return (b+3*s*X0*X0-2*X0*Y0)*pcurr[X][j]
                   - (u+X0*X0)*pcurr[Y][j]
                   - (Y0-3*s*X0)*pcurr[X][j]*pcurr[X][j]
                   - 2*X0*pcurr[X][j]*pcurr[Y][j]
                   + s*pcurr[X][j]*pcurr[X][j]*pcurr[X][j]
                   - pcurr[X][j]*pcurr[X][j]*pcurr[Y][j];
    default: throw logic_error("Invalid R-value");
  }
}

void crank_nicolson( real (*pnext)[Mx+2], const real (*pcurr)[Mx+2], int R, int t )
{
  real b[Mx], v[Mx];
  real *x = pnext[R]; // Set indices 1..Mx
  const real *y = pcurr[R];
  const real c = d[R]; // Crank-Nicolson coefficient
  *b = 2*c + 1;
  *v = c*(y[2]-2*y[1]) + y[1] + reaction(pcurr,R,1,t)/Nt;
  for( int j=1; j<Mx; ++j )
  {
    real temp = c/b[j-1];
    b[j] = *b - temp*c;
    v[j] = c*(y[j]-2*y[j+1]+y[j+2]) + y[j+1]
           + reaction(pcurr,R,j+1,t)/Nt + temp*v[j-1];
  }
  x[Mx] = v[Mx-1]/b[Mx-1];
  for( int j=Mx-2 ; j>=0; --j ) x[j+1] = ( v[j] + c*x[j+2] ) / b[j];
}

// At T=100 for 1D particle, Mathematica gives:
// Laplacians of G, X, Y: 0.675636, 55.8551, -5.33839, 
// Diffusion terms of G, X, Y: 0.675636, 55.8551, -64.0607
real laplacian( const real (*pether)[Mx+2], int R, int j )
{
  const real *y = pether[R];
  return (y[j-1]-2*y[j]+y[j+1])*Nx*Nx;
}

void next_iteration( int t /* in computational units */ )
{
  real (*pcurr)[Mx+2] = space[t&1];
  real (*pnext)[Mx+2] = space[t&1^1];

  for( int i=0 ; i<3 ; ++i ) crank_nicolson( pnext, pcurr, i, t );

  div_t tdiv = div( t, Nt );
  if( !tdiv.rem )
  {
    real lapG = laplacian( pnext, G, Mx/2 );
    real lapX = laplacian( pnext, X, Mx/2 );
    real lapY = laplacian( pnext, Y, Mx/2 );
    cerr << "time=" << tdiv.quot << ", G: " << pnext[G][Mx/2] << endl;  
    cerr << "time=" << tdiv.quot << ", X: " << pnext[X][Mx/2] << endl;  
    cerr << "time=" << tdiv.quot << ", Y: " << pnext[Y][Mx/2] << endl;  
    
    cerr << "time=" << tdiv.quot << ", lapG: " << lapG << endl;
    cerr << "time=" << tdiv.quot << ", lapX: " << lapX << endl;  
    cerr << "time=" << tdiv.quot << ", lapY: " << lapY << endl;  
    cerr << "chi(0,"<< tdiv.quot << ") = "<< chi(0,tdiv.quot) << endl;
	
	// open files and allow append to file
	std::ofstream file_G( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_G.txt", std::ios::app );
	std::ofstream file_X( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_X.txt", std::ios::app );
	std::ofstream file_Y( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_Y.txt", std::ios::app );
	std::ofstream file_lapG( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapG.txt", std::ios::app );
	std::ofstream file_lapX( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapX.txt", std::ios::app );
	std::ofstream file_lapY( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapY.txt", std::ios::app );
  	file_G << tdiv.quot << ' ' << pnext[G][Mx/2] << endl;
  	file_X << tdiv.quot << ' ' << pnext[X][Mx/2] << endl;
  	file_Y << tdiv.quot << ' ' << pnext[Y][Mx/2] << endl;
  	file_lapG << tdiv.quot << ' ' << lapG << endl;
	file_lapX << tdiv.quot << ' ' << lapX << endl;
	file_lapY << tdiv.quot << ' ' << lapY << endl;

    for( int i=0 ; i<3 ; ++i ) if( pnext[i][Mx/2] <= -vacuum[i] ) {
      cerr << "Minimum concentration value hit. Aborting." << endl;
      exit(1);
    }
  }
}

int main()
{
  //open files
  ofstream file_G;
  ofstream file_X;
  ofstream file_Y;
  ofstream file_lapG;
  ofstream file_lapX;
  ofstream file_lapY;

  ofstream file_G_r;
  ofstream file_X_r;
  ofstream file_Y_r;

  cerr << "G0 = " << G0 << endl;
  cerr << "X0 = " << X0 << endl;
  cerr << "Y0 = " << Y0 << endl;
  cerr << "Nx = " << Nx << endl;
  cerr << "Nt = " << Nt << endl;
  init();
  for( int t=0 ; t<Mt ; ++t ) next_iteration(t);
  real (*pfinal)[Mx+2] = space[Mt&1];
  for( int i=0 ; i<3 ; ++i ){
    for( int j=1 ; j<=Mx ; ++j ){

      count = count + 1;
      cout << ' ' << count << ' ';

	// open files and allow append to file
  	switch(i)
  	{
    		case 0:
		{
			std::ofstream file_G_r( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_G_r.txt", std::ios::app );
			file_G_r << j << ' ' << pfinal[i][j] << endl;
			break;
		}
    		case 1:
		{
			std::ofstream file_X_r( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_X_r.txt", std::ios::app );
			file_X_r << j << ' ' << pfinal[i][j] << endl;
			break;
		}
    		case 2:
		{
			std::ofstream file_Y_r( "/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_Y_r.txt", std::ios::app );
			file_Y_r << j << ' ' << pfinal[i][j] << endl;
			break;
		}
    		//default: throw logic_error("Invalid i-value");
  	}

      cout << pfinal[i][j] << (j<Mx?',':'\n');
    }
  }

  return 0;
}
