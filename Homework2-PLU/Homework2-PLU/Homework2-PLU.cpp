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
				printf("%.2lf    ", V(i));
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
		void copy(Vector *that) {
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
					printf("%d, %d : %.310lf\n", m, k, M(m, k));
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

#define LN 1000000

int main(int argc, char* argv[]) {

	double res[LN];
	bool newline[LN];
	bool sing[LN];
	int lines = 0;

	for (int i = 0; i < LN; ++i) {
		res[i] = 0;
		newline[i] = false;
		sing[i] = false;
	}

	while (1 == 1) {
		int n, m;
		scanf("%d", &n);
		if (n == 0) break;
		SquareMatrix A(n);
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; j <= n; ++j) {
				double num;
				scanf("%lf", &num);
				A.setM(i, j, num);
			}
		}

		PMatrix P(n);

		A.PLUDecomposite(&P);
		
		if (!A.regular()) {
			++lines;
			newline[lines - 1] = true;
			sing[lines - 1] = true;
		}
		scanf("%d", &m);
		for (int k = 1; k <= m; ++k) {
			Vector b(n);
			for (int j = 1; j <= n; ++j) {
				double num;
				scanf("%lf", &num);
				b.setV(j, num);
			}
			if (A.regular()) {
				Vector pb(n);
				b.leftMultiply(&P, &pb);
				A.solve(&pb);
				for (int i = 1; i <= n; ++i) {
					++lines;
					if (i == 1) {
						newline[lines - 1] = true;
					}
					res[lines - 1] = pb.V(i);
				}
			}
		}
	}

	for (int i = 0; i < lines; ++i) {
		if (newline[i] && i > 0) {
			printf("\n");
		}
		if (sing[i]) {
			printf("szingularis");
		} else {
			if (!newline[i]) {
				printf(" ");
			}
			printf("%.8lf", res[i]);
		}
	}
	printf("\n");
	//system("PAUSE");
	return 0;
}
