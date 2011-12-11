#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Superclass for vectors and SquareMatrix matrixes
 */

class Array {
	private:
		int _n;
	public:
		Array(int n) {
			_n = n;
		}
		int n() {
			return _n;
		}
};

/*
 * Abstract class for matrixes
 */

class Matrix : public Array {
	public:
		Matrix(int n) : Array(n) {}
		virtual double M(int i, int j) = 0;
		virtual void setM(int i, int j, double val) = 0;
		void dump() {
			for (int i = 1; i <= n(); ++i) {
				printf("\n");
				for (int j = 1; j <= n(); ++j) {
					printf("%.2lf     ", M(i, j));
				}
			}
			printf("\n");
		}
};

class PMatrix : public Matrix {
	private:
		int *_r;
	public:
		PMatrix(int n) : Matrix(n) {
			_r = (int *) calloc(n, sizeof(int));
			for (int i = 0; i < n; ++i) {
				_r[i] = i + 1;
			}
		}
		virtual double M(int i, int j) {
			if (i < 1 || j < 1 || i > n() || j > n()) {
				throw "Wrong index!";
			}
			return _r[i - 1] == j ? 1 : 0;
		}
		virtual void setM(int i, int j, double val) {
			throw "Unsupported action!";
		};
		void swapRows(int i1, int i2) {
			int tmp = _r[i1 - 1];
			_r[i1 -1] = _r[i2 -1];
			_r[i2 -1] = tmp;
		}
};

class Vector : public Array {
	private:
		double *_v;
	public:
		Vector() : Array(0) {}
		Vector(int n) : Array(n) {
			_v = (double *) calloc(n, sizeof(double));
		}
		~Vector() {
			free(_v);
		}
		double V(int i) {
			return _v[i - 1];
		}
		void setV(int i, double val) {
			_v[i - 1] = val;
		}
		void dump() {
			for (int i = 1; i <= n(); ++i) {
				printf("%.4lf    ", V(i));
			}
			printf("\n");
		}
		void leftMultiply(Matrix *A, Vector *res) {
			for (int i = 1; i <= n(); ++i) {
				double a = 0;
				for (int j = 1; j <= n(); ++j) {
					a += A -> M(i, j) * V(j);
					res -> setV(i, a);
				}
			}
		}
		void copyTo(Vector *that) {
			for (int i = 1; i <= n(); ++i) {
				that -> setV(i, V(i));
			}
		}
		void print() {
			for (int i = 1; i <= n(); ++i) {
				if (i > 1) {
					printf(" ");
				}
				printf("%.8lf", V(i));
			}
		}
		bool nullVector() {
			for (int i = 1; i <= n(); ++i) {
				if (fabs(V(i)) > 1e-30) return false;
			}
			return true;
		}
		double norm2() {
			double a = 0;
			for (int i = 1; i <= n(); ++i) {
				a += V(i) * V(i);
			}
			return sqrt(a);
		}
		double normInf() {
			double ret = fabs(V(1));
			for (int i = 2; i <= n(); ++i) {
				if (ret < fabs(V(i))) ret = fabs(V(i));
			}
			return ret;
		}
		void multiplyScalar(double c) {
			for (int i = 1; i <= n(); ++i) {
				setV(i, V(i) * c);
			}
		}
		double scalarProduct(Vector *w) {
			double a = 0;
			for (int i = 1; i <= n(); ++i) {
				a += V(i) * w -> V(i);
			}
			return a;
		}
		void minus(Vector * w, Vector *res) {
			for (int i = 1; i <= n(); ++i) {
				res -> setV(i, V(i) - w -> V(i));
			}
		}
};

class SquareMatrix : public Matrix {
	private:
		double **_m;
		bool _reg;
	public:
		SquareMatrix(int n) : Matrix(n) {
			_m = (double **) calloc(n, sizeof(double *));
			for (int i = 0; i < n; ++i) {
				_m[i] = (double *) calloc(n, sizeof(double));
			}
		}
		~SquareMatrix() {
			for (int i = 0; i < n(); ++i) {
				free(_m[i]);
			}
			free(_m);
		}
		virtual double M(int i, int j) {
			return _m[i - 1][j - 1];
		}
		virtual void setM(int i, int j, double val) {
			if (i < 1 || j < 1 || i > n() || j > n()) {
				throw "Wrong index!";
			}
			_m[i - 1][j - 1] = val;
		}
		void swapRows(int i1, int i2) {
			for (int i = 1; i <= n(); ++i) {
				double tmp = M(i1, i);
				setM(i1, i, M(i2, i));
				setM(i2, i, tmp);
			}
		}
		void PLUDecomposite(PMatrix *P) {
			for (int k = 1; k <= n()-1; ++k) {
				int m = k;
				for (int i = k+1; i <= n(); ++i) {
					if (fabs(M(i, k)) > fabs(M(m, k))) {
						m = i;
					}
				}
				if (fabs(M(m, k) + 0.0) < 1e-15) {
					_reg = false;
					return;
				}
				if (m != k) {
					swapRows(m, k);
					P -> swapRows(m, k);
				}
				for (int i = k+1; i <= n(); ++i) {
					setM(i, k, M(i, k) / M(k, k));
					for (int j = k+1; j <= n(); ++j) {
						setM(i, j, M(i, j) - M(i, k) * M (k, j));
					}
				}
			}
			_reg = fabs(M(n(), n())) >= 1e-15;
		}
		bool regular() {
			return _reg;
		}
		void solve(Vector *b) {
			for (int i = 1; i <= n(); ++i)
				for (int j = 1; j <= i-1; ++j) {
					b -> setV(i, b -> V(i) - M(i, j) * b -> V(j));
				}
			for (int i = n(); i >= 1; --i) {
				for (int j = i+1; j <= n(); ++j) {
					b -> setV(i, b -> V(i) - M(i, j) * b -> V(j));
				}
				b -> setV(i, b -> V(i) / M(i, i));
			}
		}
		void copyTo(SquareMatrix *that) {
			for (int i = 1; i <= n(); ++i)
				for (int j = 1; j <= n(); ++j) {
				that -> setM(i, j, M(i, j));
			}
		}
};

class LMatrix : public Matrix {
	private:
		SquareMatrix *_M;
	public:
		LMatrix(SquareMatrix *M) : Matrix(M -> n()) {
			_M = M;
		}
		virtual double M(int i, int j) {
			if (i == j)	{
				return 1;
			}
			if (i < j) {
				return 0;
			}
			return _M -> M(i, j);
		}
		virtual void setM(int i, int j, double val) {
			if (i <= j) {
				throw "I wouldn't try that!";
			}
			_M -> setM(i, j, val);
		}
};

class UMatrix : public Matrix {
	private:
		SquareMatrix *_M;
	public:
		UMatrix(SquareMatrix *M) : Matrix(M -> n()) {
			_M = M;
		}
		virtual double M(int i, int j) {
			if (i > j) {
				return 0;
			}
			return _M -> M(i, j);
		}
		virtual void setM(int i, int j, double val) {
			if (i > j) {
				throw "I wouldn't try that!";
			}
			_M -> setM(i, j, val);
		}
};

class Equation {
	private:
		int num;
	public:
		Equation(int num) {
			this -> num = num;
		}
		int dim() {
			switch (num) {
				case 1: return 3;
				case 2: return 2;
				case 3: return 2;
			}
			return -1;
		}
		void F(Vector *x, Vector *y) {
			switch (num) {
				case 1:
					y -> setV(1, - x->V(1) * x->V(1) + x->V(3) + 3);
					y -> setV(2, - x->V(1) + 2 * x->V(2) * x->V(2) - x->V(3) * x->V(3) - 3);
					y -> setV(3, x->V(2) - 3 * x->V(3) * x->V(3) + 2);
					break;
				case 2:
					y -> setV(1, 2 * x->V(1) * x->V(1) - x->V(2) - 1);
					y -> setV(2, - x->V(1) + 2 * x->V(2) * x->V(2) - 1);
					break;
				case 3:
					y -> setV(1, -4.0 * x->V(1) + cos(2.0 * x->V(1) - x->V(2)) - 3.0);
					y -> setV(2, sin(x->V(1)) - 3.0 * x->V(2) - 2.0);
					break;
			}
		}
		void jacobi(Vector *x, Matrix *J) {
			switch (num) {
				case 1:
					J -> setM(1, 1, -2 * x->V(1)); J -> setM(1, 2, 0           ); J -> setM(1, 3, 1           );
					J -> setM(2, 1, -1          ); J -> setM(2, 2, 4 * x->V(2) ); J -> setM(2, 3, -2 * x->V(3));
					J -> setM(3, 1, 0           ); J -> setM(3, 2, 1           ); J -> setM(3, 3, -6 * x->V(3));
					break;
				case 2:
					J -> setM(1, 1, 4 * x->V(1)); J -> setM(1, 2, -1         );
					J -> setM(2, 1, -1         ); J -> setM(2, 2, 4 * x->V(2));
					break;
				case 3:
					J -> setM(1, 1, -4.0 - 2.0 * sin(2.0 * x->V(1) - x->V(2))); J -> setM(1, 2, sin(2.0 * x->V(1) - x->V(2)));
					J -> setM(2, 1, cos(x->V(1))                       ); J -> setM(2, 2, -3.0                        );
					break;
			}
		}
};

int main(int argc, char* argv[]) {

	int N;
	scanf("%d", &N);

	for (int s0 = 0; s0 < N; ++s0) {
		int num;
		scanf("%d", &num);
		int maxit;
		scanf("%d", &maxit);
		double eps;
		scanf("%lf", &eps);
		Equation eq(num);
		int n = eq.dim();
		Vector x(n);
		for (int i = 1; i <= n; ++i) {
			double num;
			scanf("%lf", &num);
			x.setV(i, num);
		}
		Vector y0(n);
		eq.F(&x, &y0);
		double t = 1;
		double norm0 = y0.normInf();
		for (int k = 1; k <= maxit; ++k) {
			eq.F(&x, &y0);
			double normy = y0.normInf();
			SquareMatrix J(n);
			eq.jacobi(&x, &J);
			PMatrix P(n);
			J.PLUDecomposite(&P);
			if (!J.regular()) {
				printf("szingularis");
				for (int i = 1; i <= n; ++i) {
					printf(" %.8lf", x.V(i));
				}
				printf(" %.8lf\n", normy);
				break;
			}
			if (normy <= eps * (1 + norm0)) {
				printf("siker");
				for (int i = 1; i <= n; ++i) {
					printf(" %.8lf", x.V(i));
				}
				printf(" %.8lf %d\n", y0.normInf(), k - 1);
				break;
			}
			Vector dx(n);
			for (int i = 1; i <= n; ++i) {
				y0.setV(i, -y0.V(i));
			}
			y0.leftMultiply(&P, &dx);
			J.solve(&dx);
			double xnorm = x.normInf();
			bool br = false;
			for (int l = 0; l < 9; ++l) {
				Vector y(n);
				for (int i = 1; i <= n; ++i) {
					y.setV(i, x.V(i) + t * dx.V(i));
				}
				Vector fy(n);
				eq.F(&y, &fy);
				if (fy.normInf() < normy) {
					if (l == 0) {
						t = (1.5 * t) > 1 ? 1 : (1.5 * t);
					}
					y.copyTo(&x);
					break;
				}
				t = t / 2;
				if (t < 1e-3) {
					br = true;
					break;
				}
				if (l == 8) {
					br = true;
					break;
				}
			}
			if (br) {
				printf("sikertelen");
				for (int i = 1; i <= n; ++i) {
					printf(" %.8lf", x.V(i));
				}
				printf(" %.8lf\n", normy);
				break;
			}
			if (k == maxit) {
				printf("maxit\n");
			}
		}
	}

	return 0;
}