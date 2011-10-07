#include <stdio.h>
#include <iostream>

#define DIM 100

/*
 * Superclass for vectors and square matrixes
 * n: Dimension
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
		Matrix(int n) : Array(n) {};
		virtual double M(int i, int j) = 0;
		virtual void setM(int i, int j, double val) = 0;
};

class SimpleMatrix : public Matrix {
	private:
		double _m[DIM][DIM];
	public:
		SimpleMatrix(int n) : Matrix(n) {};
		virtual double M(int i, int j) {
			return _m[i - 1][j - 1];
		}
		virtual void setM(int i, int j, double val) {
			if (i < 1 || j < 1 || i > n() || j > n())
				throw "Wrong index!";
			_m[i - 1][j - 1] = val;
		}
};

class LMatrix : public Matrix {
	private:
		SimpleMatrix *_M;
	public:
		LMatrix(SimpleMatrix *M) : Matrix(M -> n()) {
			_M = M;
		}
		virtual double M(int i, int j) {
			if (i == j)
				return 1;
			if (i < j)
				return 0;
			return _M -> M(i, j);
		}
		virtual void setM(int i, int j, double val) {
			if (i <= j)
				throw "I wouldn't try that!";
			_M -> setM(i, j, val);
		}
};

class UMatrix : public Matrix {
	private:
		SimpleMatrix *_M;
	public:
		UMatrix(SimpleMatrix *M) : Matrix(M -> n()) {
			_M = M;
		}
		virtual double M(int i, int j) {
			if (i > j)
				return 0;
			return _M -> M(i, j);
		}
		virtual void setM(int i, int j, double val) {
			if (i > j)
				throw "I wouldn't try that!";
			_M -> setM(i, j, val);
		}
};



class Vector : public Array {
	private:
		double _v[DIM];
	public:
		Vector(int n) : Array(n) {};
		double V(int i) {
			return _v[i - 1];
		}
		void setV(int i, double val) {
			_v[i - 1] = val;
		}
};

int main(int argc, char* argv[]) {

	while (1 == 1) {
		int n, m, i, j;
		double num;
		scanf("%d", &n);
		if (n == 0) break;
		SimpleMatrix A(n);
		for (i = 1; i <= n; ++i) {
			for (j = 1; j <= n; ++j) {
				scanf("%lf", &num);
				A.setM(i, j, num);
			}
		}
		scanf("%d", &m);
		Vector b(n);
		for (i = 1; i <= m; ++i) {
			for (j = 1; j <= n; ++j) {
				scanf("%lf", &num);
				b.setV(j, num);
			}
		}

		LMatrix L(&A);
		UMatrix U(&A);

		for (i = 1; i <= A.n(); ++i) {
			printf("\n");
			for (j = 1; j <= A.n(); ++j)
				printf("%.3lf / ", A.M(i, j));
		}
		printf("\n");
		
		for (i = 1; i <= L.n(); ++i) {
			printf("\n");
			for (j = 1; j <= L.n(); ++j)
				printf("%.3lf / ", L.M(i, j));
		}
		printf("\n");
		
		for (i = 1; i <= U.n(); ++i) {
			printf("\n");
			for (j = 1; j <= U.n(); ++j)
				printf("%.3lf / ", U.M(i, j));
		}
		printf("\n");
		
	} 

	//system("PAUSE");
	return 0;
}
