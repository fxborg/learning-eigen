#include "ssa.h"
EXPORT CSSA * __stdcall Create(const unsigned int size, const unsigned int length, const  double omega0, const unsigned int max_et, const double c0step) {
	std::cout << "OK!";
	return new CSSA(size, length, omega0, max_et, c0step);
}

EXPORT void __stdcall Destroy(CSSA * instance)
{
	delete instance;
}

EXPORT int __stdcall Push(CSSA * instance, const int x, const double y, const time_t t0, const time_t t1)
{
	return instance->push(x, y, t0, t1);
}
EXPORT int __stdcall Calculate(CSSA* instance)
{
	return instance->calculate();
}

CSSA::CSSA(const unsigned int size, const unsigned int length, const  double omega0, const unsigned int max_et, const double c0step) :
	m_size(minmax(size, 50, 2000)),
	m_L(minmax(length, 5, m_size)),
	m_omega0(minmax(omega0, 0.001, 0.5)),
	m_max_et(minmax(max_et, 2, int(0.5*m_L))),
	m_c0min(0.5),
	m_c0max(1.0),
	m_c0step(c0step),
	m_rdelta_threshold(0.05),
	m_c0eps(0),
	m_series(CSeries(m_size + 1))
{}


int CSSA::push(const int x, const double y, const time_t t0, const time_t t1)
{
	int result = 0;
	try
	{
		result = m_series.push(x, y, t0, t1);
	}
	catch (...)
	{
		result = -9999;
	}

	return result;
}

int CSSA::calculate()
{
	std::cout << "calc" << std::endl;
	if (!m_series.is_adding())return false;
	const std::deque<double> & series = m_series.get_stats_series();
	std::cout << "!!!" << std::endl;

	double v = calc_cmax(series);
	return 0;
}


void CSSA::ssa(const std::deque<double> & series, cv::Mat1d &S, cv::Mat1d &V, cv::Mat1d &U)
{
	int k = m_size - m_L + 1;
	int L = int(m_L);
	std::cout << "m_size=" << m_size << std::endl;
	std::cout << "size=" << series.size() << std::endl;
	std::cout << "L=" << L << std::endl;
	std::cout << "k=" << k << std::endl;

	std::vector<double> vec;
	vec.reserve(L*k);
	int adj = series.size() - m_size;
	auto it = series.cbegin() + adj;

	for (int i = 0; i < L; i++)
	{
		std::copy(it, it + k, std::back_inserter(vec));
		++it;
	}

	cv::Mat1d X(L, k, vec.data());

	std::cout << "X=" << X << std::endl << std::endl;
	std::cout << "series=" << std::endl;
	for (auto it = series.begin(); it != series.end(); ++it)
		std::cout << *it << ",";
	std::cout << std::endl;

	// SVD
	cv::SVD::compute(X, S, U, V);
	cv::Mat1b nonZeroSingularValues = S > 0.0001;
	// rank
	int M = cv::countNonZero(nonZeroSingularValues);
	if (M < L)
	{
		S = S.rowRange(0, M);
		U = U.colRange(0, M);
		V = V.colRange(0, M);
	}
	std::cout << M;

	std::cout << "------------------S---------------" << std::endl;
	std::cout << S.size() << std::endl;


	std::cout << "-----------------U--------------" << std::endl;
	std::cout << U.size() << std::endl;
	std::cout << "------------------V--------------" << std::endl;
	std::cout << V.size() << std::endl;

}



double CSSA::calc_cmax(const std::deque<double> & series)
{
	cv::Mat1d sing_values, U, V, trend;
	ssa(series, sing_values,U,V);
	int d = std::min(sing_values.rows, int(m_max_et));

	//  Calculate LowFreq values for all ETs
	std::vector<double>cval_forET(d + 1);
	for (int et = 1; et <= d; et++)
	{
		cval_forET[et] = LFvalue_vect(U, m_omega0, et);
	}
	//---
	int c0s_amount = floor((m_c0max - m_c0min) / m_c0step) + 1;
	int c0_i =1;
	for (double c0 = m_c0min; c0 <= m_c0max; c0 += m_c0step)
	{
		std::vector<int> trend_ETs;
		for (int et = 1; et <= d; et++)
		{
			if (cval_forET[et] > c0) trend_ETs.push_back(et);
		}
		reconstruct(sing_values, U, V, trend_ETs, trend);
		if (c0 > 0)
		{
		}
		else
		{
		}
	}


	return 0.0;
}

unsigned int CSSA::minmax(const unsigned int n, const unsigned int min, const unsigned int max)
{
	return std::max(min, std::min(max, n));
}

double CSSA::minmax(const double n, const double min, const double max)
{
	return std::max(min, std::min(max, n));
}

double CSSA::periodogram(const std::valarray<double>& F, const int k)
{

	int N = F.size();
	double w = (k - 1) / N;
	std::valarray<double> T(N);
	for (int i = 0; i < N; ++i) T[i] = (double)i;


	std::valarray<double> c = (F * cos(M_PI_2*T*w));
	std::valarray<double> s = (F * sin(M_PI_2*T*w));
	double csum = c.sum();
	double ssum = s.sum();
	double Pi = csum*csum + ssum*ssum;

	Pi = Pi * 2.0 / N;

	if (k == 1 || k == floor(N / 2) + 1)
		Pi = Pi / 2;
	return Pi;
}

double CSSA::LFvalue_vect(const cv::Mat1d & M, const double omega0, const int ET)
{
	int L = M.cols;
	double nrm = cv::norm(M.rowRange(ET, ET));
	if (nrm == 0.0)
	{
		return -0.0;
	}

	std::valarray<double> vect(std::vector<double>(M.rowRange(ET, ET)).data(), L);
	vect = vect / nrm;
	// calculate sum of periodogram values for low frequencies
	int k = 1;
	double w = 0;
	double val = 0;
	while (w < omega0)
	{
		val += periodogram(vect, k);
		k = k + 1;
		w = (k - 1) / L;
	}


	return 0.0;
}

void CSSA::reconstruct(const cv::Mat1d & sing_values, const cv::Mat1d & U, const cv::Mat1d & V, const std::vector<int> ETs,  cv::Mat1d & F)
{
	int L = U.size().height;
	int K = V.size().height;
	int N = K + L - 1;
	cv::Mat1d X(L, K);
	int sz = ETs.size();
	for (int i = 0; i < sz; i++)
	{
		int ET = ETs[i];
		const cv::Mat1d sing_value{ sing_values.rowRange(ET, ET) };
		const cv::Mat1d u{ U.rowRange(ET, ET) };
		const cv::Mat1d v{ V.rowRange(ET, ET)};
		const cv::Mat1d Xi{ sing_value * u * v.t() };
		X += Xi;
	}
	cv::Mat1d X(N, 1);
	double L0 = std::min(L, K);
	double K0 = std::max(L, K);
 
	cv::Mat1d Y = (L<K) ? X : X.t();


	for (int k = 1; k < L0; k++)
	{
		int d = 0;
		for (int m = 1; m <= k; m++)
		{
			d +=  Y(m, k - m + 1);
		}
		F(k) = 1 / k * d;
	}
	for (int k = L0; k <= K0; k++)
	{

		int d = 0;
		for (int m = 1; m <= L0; m++)
		{
			d += Y(m, k - m + 1);
		}
		F(k) = 1 / k * d;
	}
	for (int k = K0 + 1; k <= N; k++)
	{
		int d = 0;
		for (int m = k - K0; m <= N - K - +1; m++)
		{
			d += Y(m, k - m + 1);
		}
		F(k) = 1 / (N-k+1) * d;
	}

}

