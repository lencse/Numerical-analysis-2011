#include <math.h>
#include <stdio.h>

#define EE4 10000

int main(int argc, char* argv[]) {
	int n;
	double pi = 2 * acos(0.0);
	double res[EE4];
	scanf("%d", &n);
	for (int i = 0; i < n; ++i) {
		double v[EE4];
		double w[EE4];
		int dim;
		scanf("%d", &dim);
		for (int j = 0; j < dim; ++j) scanf("%lf", &v[j]);
		for (int j = 0; j < dim; ++j) scanf("%lf", &w[j]);
		double sp = 0;
		for (int j = 0; j < dim; ++j) sp += v[j] * w[j];
		double vnorm = 0;
		for (int j = 0; j < dim; ++j) vnorm += v[j] * v[j];
		vnorm = sqrt(vnorm);
		double wnorm = 0;
		for (int j = 0; j < dim; ++j) wnorm += w[j] * w[j];
		wnorm = sqrt(wnorm);
		res[i] = acos(sp / (vnorm * wnorm)) * 180 / pi;
	}
	for (int i = 0; i < n; ++i)
		printf("%.9lf\n", res[i]);
	return 0;
}
