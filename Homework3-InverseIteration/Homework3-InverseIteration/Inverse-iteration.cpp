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

int main(int argc, char* argv[]) {

	int N;
	scanf("%d", &N);

	for (int s0 = 0; s0 < N; ++s0) {
		int n;
		scanf("%d", &n);
		double c;
		scanf("%lf", &c);
		int maxit;
		scanf("%d", &maxit);
		double eps;
		scanf("%lf", &eps);
		SquareMatrix A(n);
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; j <= n; ++j) {
				double num;
				scanf("%lf", &num);
				A.setM(i, j, num);
			}
		}
		Vector x(n);
		for (int i = 1; i <= n; ++i) {
			double num;
			scanf("%lf", &num);
			x.setV(i, num);
		}
		SquareMatrix A1(n);
		A.copyTo(&A1);
		for (int i = 1; i <= n; ++i) {
			A1.setM(i, i, A.M(i, i) - c);
		}

		PMatrix P(n);
		A1.PLUDecomposite(&P);

		if (!A1.regular()) {
			printf("%.8lf\n", c);
			continue;
		}
		if (x.nullVector()) {
			printf("kezdovektor\n", c);
			continue;
		}

		x.multiplyScalar(1 / x.norm2());
		
		Vector Ax(n);
		x.leftMultiply(&A, &Ax);
		double l0 = x.scalarProduct(&Ax);
		double l1;

		int m;
		for (m = 1; m <= maxit; ++m) {
			Vector y(n);
			x.leftMultiply(&P, &y);
			A1.solve(&y);
			y.copyTo(&x);
			x.multiplyScalar(1 / x.norm2());
			x.leftMultiply(&A, &Ax);
			l1 = x.scalarProduct(&Ax);
			if (fabs(l1 - l0) < eps * (1 + fabs(l1)))	break;
			l0 = l1;
		}
		
		if (m > maxit) {
			printf("maxit\n");
			continue;
		}

		double d = 0;

		for (int i = 1; i <= n; ++i) {
			d += (Ax.V(i) - l1 * x.V(i)) * (Ax.V(i) - l1 * x.V(i));
		}


		printf(d <= eps ? "siker" : "sikertelen");
		printf(" %.8lf", l1);
		for (int i = 1; i <= n; ++i) {
			printf(" %.8lf", x.V(i));
		}
		printf(" %.8lf", d);
		printf("\n");
	}

	return 0;
}