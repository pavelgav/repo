#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <sys/time.h>

using namespace std;

/* 3-diagonal matrix class */
template <typename T> class DSMatrix
{
private:
	int msize;
	vector<T> v1, v2;
public:
	DSMatrix() {}

	DSMatrix( int n ) 
	{
		msize = n;
		v1.assign(n,0);
		v2.resize(n-1,0);
	}
	~DSMatrix(){}

	void resize( int n, T value_d = 0, T value_s = 0)
	{
		msize = n;
		v1.resize(n, value_d);
		v2.resize(n-1, value_s);
		return;	
	}

	int size() const
	{
		return msize;
	}

	void set_d (const int i, const T obj )
	{
		v1[i] = obj;
		return;
	}
	void set_s (const int i, const T obj )
	{
		v2[i] = obj;
		return;
	}

	const T get_d ( const int i )
	{
		return v1[i];
	}

	const T get_s ( const int i )
	{
		return v2[i];
	}

	DSMatrix operator+ ( const DSMatrix<T>& r )
	{
		int N = r.size();
		DSMatrix temp(N);
		for (int i = 0; i < N-1; ++i)
		{
			temp.v1[i] = v1[i] + r.v1[i];
			temp.v2[i] = v2[i] + r.v2[i];
		}
		temp.v1[N-1] = v1[N-1] + r.v1[N-1];
		return temp;
	}

	friend DSMatrix operator- (const DSMatrix<T>& l, const DSMatrix<T>& r)
	{
		int N = l.size();
		DSMatrix temp(N);
		for (int i = 0; i < N; ++i)
		{
			temp.v1[i] = l.v1[i] - r.v1[i];
			temp.v2[i] = l.v2[i] - r.v2[i];
		}
		return temp;
	}

	friend DSMatrix operator* (const T& x, const DSMatrix<T>& m )
	{
		int N = m.size();
		DSMatrix temp(N);
		for (int i = 0; i < N-1; ++i)
		{
			temp.v1[i] = x*m.v1[i];
			temp.v2[i] = x*m.v2[i];
		}
		temp.v1[N-1] = x*m.v1[N-1];
		return temp;
	}

	vector<T> operator* ( vector<T> v )
	{
		int N = v.size();
		vector<T> temp(N);

		temp[0] = v1[0]*v[0] + v2[0]*v[1];
		for (int i = 1; i < N-1; ++i)
		{
			temp[i] = v2[i-1]*v[i-1] + v1[i]*v[i] + v2[i]*v[i+1];
		}
		temp[N-1] = v2[N-2]*v[N-2] + v1[N-1]*v[N-1];
		return temp;
	}

	DSMatrix& operator= (const DSMatrix& r)
	{
		this->resize( r.size() );
		v1 = r.v1;
		v2 = r.v2;
		return *this;
	}
};



double timer()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + 1e-6 * (double)ts.tv_usec;
}




/* Grid of a problem */
class grid
{
public:
	vector<double> x;
	int Nx;
	double l, h, mu, omega;
	DSMatrix <complex<double> > M, K;
	complex<double> f;


	grid( int Nx_, double l_, double omega_, double mu_ ) :
		Nx(Nx_), 
		l(l_), 
		h( l_/Nx_ ), 
		omega(omega_),
		mu(mu_)
	{
		for (int i = 0; i <= Nx; ++i)
		{
			x.push_back( i*h );
		}

		M.resize(Nx);
		K.resize(Nx);
		init_matrix();
	}

	/* M, D, K initialization */
	void init_matrix ()
	{
		for (int i = 0; i < Nx-1; ++i)
		{
			double gm1 = 1/pow(c( (i+0.5)*h ), 2) , gm2 = 1/pow(c( (i+1.5)*h ), 2);
			complex<double> zt1 ( zt( (i+0.5)*h ) ), zt2( zt( (i+1.5)*h ) ); 

			M.set_d( i, h*(gm1*zt1+gm2*zt2)/3. );
			M.set_s( i, h*gm2*zt2/6. );

			K.set_d( i, ( 1./zt1 + 1./zt2 )/h );
			K.set_s( i, -1./zt2/h );
		}
		double gm1 = 1/pow(c( (Nx-0.5)*h ), 2);
		complex<double> zt1( zt( (Nx-0.5)*h ) );		
		//cout << omega << endl;

		M.set_d( Nx-1, h*zt1*gm1/3. );
		K.set_d( Nx-1, 1./zt1/h );

		gm1 = 1/pow( c(0.5*h), 2 );
		zt1 = zt(0.5*h);
		f = h*gm1*zt1*pow(omega,2)/6. + 1./zt1/h;
		//f += (0.004, 0.);

		//cout << K.get_d(Nx-1) << endl;
		return;
	}

	void set_mu(double mu_)
	{
		mu = mu_;
	}

	void set_omega(double omega_)
	{
		omega = omega_;
	}


	void refine_grid(int N_)
	{
		Nx = N_;
		h = l/Nx;
	}


	/* Sound velocity */
	double c( double x )
	{
		return ( x <= 1 ) ? ( ( 7 + cos(4*M_PI*x) )/8 ) : ( 1 );
	}

	/* Sigma function */
	double sgm( double x )
	{
		return ( x <= 1 ) ? ( 0 ) : ( mu*exp(x-1) );
	}

	complex<double> zt(double x)
	{
		complex<double> im(0., 1.);
		return (omega - im*sgm(x) )/omega;
	}

	~grid(){}
	
};






/* Solver of the problem */
class solver
{
private:
	grid &g;
	double t, tau;
	vector<complex<double> > b, p_temp, r;
	int step;
	DSMatrix<complex<double> > A;
public:
	solver( grid& g_) :
		t(0.),
		step(0), 
		g(g_)
	{
		A.resize(g.Nx);
		b.assign(g.Nx, 0);
		p_temp.assign(g.Nx, 0);
		r.resize(g.Nx);
		tau = g.h/2;

		build_system();
	}

	~solver(){}

	/* Building of linear system to solve */
	void build_system()
	{

		int N = g.Nx;
		complex<double> im(0, 1.);
		A = (-pow(g.omega, 2) )*g.M + g.K;

		b[0] = g.f;

		//cout << A.get_s(N-2) << endl;
		return;
	}


	/* Progonka */
	void progonka( )
	{
		int N = g.Nx;
		vector<complex<double> > p(N), q(N);

		p[0] = -A.get_s(0)/A.get_d(0);
		q[0] = b[0]/A.get_d(0);
		//cout << q[0] << endl;
		for (int i = 1; i < N-1; ++i)
		{
			p[i] = -A.get_s(i)/( A.get_d(i) + A.get_s(i-1)*p[i-1] );
			q[i] = (-A.get_s(i-1)*q[i-1] + b[i])/( A.get_d(i) + A.get_s(i-1)*p[i-1] );
			//cout << p[i] << " " << q[i] << endl;
		}
		p_temp[N-1] = ( -A.get_s(N-2)*q[N-2] + b[N-1] )/( A.get_d(N-1) + A.get_s(N-2)*p[N-2] );
		for (int i = N-1; i > 0; i--)
		{
			p_temp[i-1] = p[i-1]*p_temp[i] + q[i-1];
		}
		//cout << p_temp[0] << endl;
		return ;
	}

	void reduction ( DSMatrix<complex<double> >& M, vector<complex<double> >& v, vector<complex<double> >& r)
	{
		int N = g.Nx;
		double p = ceil(log(N-1)/log(2.));
		int K = (int)(pow(2, p) );
		vector<complex<double> > a(K+1), b(K+1), c(K+1), d(K+1);

		M.resize(K+1, 1, 0);
		v.resize(K+1, 0);
		r.resize(K+1);

		//cout << M.get_d(K) << endl;

		for (int i = 0; i < K; ++i)
		{
			b[i] = M.get_d(i);
			c[i] = M.get_s(i);
			a[i+1] = M.get_s(i);
		}
		b[K] = A.get_d(K); d = v;

		//cout << c[K-1] << endl;
		for (int s = 1; s <= K/2; s *= 2 )
		{
			complex<double> C (-c[0]/b[s] );
			b[0] += C*a[s];
			d[0] += C*d[s];
			c[0] = C*c[s];
			complex<double> A (-a[K]/b[K-s] );
			b[K] += A*c[K-s];
			d[K] += A*d[K-s];
			a[K] = A*a[K-s];
			for (int j = 2*s; j <= K-2*s; j += 2*s)
			{
				A = -a[j]/b[j-s];
				C = -c[j]/b[j+s];
				b[j] += A*c[j-s] + C*a[j+s];
				d[j] += A*d[j-s] + C*d[j+s];
				a[j] = A*a[j-s];
				c[j] = C*c[j+s];
			}
		}

		r[0] = (d[0]*b[K] - c[0]*d[K])/(b[0]*b[K] - a[K]*c[0]);
		r[K] = (b[0]*d[K] - a[K]*d[0])/(b[0]*b[K] - a[K]*c[0]);

		for (int s = K/2; s >= 1; s = s/2)
		{
			for (int j = s; j <= K-s; j += 2*s)
			{
				r[j] = (d[j] - a[j]*r[j-s] - c[j]*r[j+s])/b[j];
			}
		}
		return;
	}

	/* Time integration */
	void integrate ( double t_stop )
	{
		//progonka();
		double t0 = timer();
		reduction(A, b, r);
		//cout << timer() - t0 << endl;

		save(g.mu);
		return;
	}


	/* Output */
	void save(double mu)
	{
		ofstream fout;
		char buf[128];
		int N = g.Nx;
		double h = g.h;
		int scale = 8;
		sprintf(buf,"./out/output%d.csv", scale);
		fout.open(buf, std::ios::out);

		fout << "x"  << "," << "Im" << endl;
		fout << 0  << "," << 0 << endl;
		complex <double> im(0, 1.);
		complex <double> w(exp(im*g.omega*10.) ); 
		//cout << w << endl;
		for (int i = 0; i < N; i += scale)
		{
			fout << g.x[i] + scale*h  << "," << r[i].imag() << '\n';
		}
		fout.close();
	}

	
};


int main(int argc, char const *argv[])
{
	grid G( 16000, 2, 5, 10 );
	solver S( G );
	S.integrate(10);
	/*for (int i = 1; i < 11; ++i)
	{
		G.set_mu(100*i);
		G.init_matrix();
		S.build_system();
		S.integrate(10);
	}*/
	

	return 0;
}