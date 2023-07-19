#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <config.h>
#include <math.h>
#include <windows.h>
#include <ppl.h>
#include <Bluestein.h>
#include <complex>
#include <Eigen/../unsupported/Eigen/FFT>
#include <fftw3.h>

#include <vector>


// Inspired by code found in https://gist.github.com/endolith/2783807

Eigen::MatrixXcd chirpZ(Eigen::MatrixXcd x,  std::complex<double> A,  std::complex<double> W, int M) {
	const std::complex<double> i(0.0, 1.0);
	int N_size = x.size();
	int N_rows = x.rows();
	int N_cols = x.cols();
	Eigen::FFT<double> fft;

	int L = pow(std::ceil(log2(M + N_size - 1)), 2);

	Eigen::VectorXcd yA(L);
	Eigen::VectorXcd yW(L);
	Eigen::VectorXcd Y(L);
	Eigen::VectorXcd V(L);
	//std::unique_ptr<Eigen::VectorXcd> yA(new Eigen::VectorXcd(L));
	//std::unique_ptr<Eigen::VectorXcd> yW(new Eigen::VectorXcd(L));
	//std::unique_ptr<Eigen::VectorXcd> Y(new Eigen::VectorXcd(L));
	//std::unique_ptr<Eigen::VectorXcd> V(new Eigen::VectorXcd(L));

	Eigen::VectorXcd yTotal(L);
	for (int j = 0; j < L; ++j) {
		yA[j] = pow(A, -j);
		yW[j] = pow(W, pow(j, 2) / 2.);
	}
	yTotal = yA.array() * yW.array() * x.array();
	fft.fwd(Y, yTotal);

	Eigen::VectorXcd v(L);
	for (int m = 0; m < M; ++m) {
		v[m] = pow(W, pow(-1. * m, 2) / 2.);
	}
	for (int p = L - M + 1; p < L; ++p) {
		v[p] = pow(W, pow(-1. * (N_size - p - 1.0), 2) / 2.);
	}
	fft.fwd(V, v);

	Eigen::VectorXcd VY(M);
	VY = Y.array() * V.array();
	Eigen::VectorXcd g(M);
	fft.inv(g, VY);

	for (int m = 0; m < M; ++m) {
		g[m] *= pow(W, pow(m, 2) / 2);
	}

	return g;

}

