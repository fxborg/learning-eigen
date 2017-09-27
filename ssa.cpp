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
{
	std::cout << m_omega0 << std::endl;
}


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
	if (!m_series.is_adding())return 0;
	const std::deque<double> & series = m_series.get_stats_series();
	if (series.size() != m_size+1)return 0;
	double v = calc_cmax(std::vector<double>(series.cbegin(), series.cend()- 1));
	return 0;
}


void CSSA::ssa(const std::vector<double> & series, cv::Mat1d &S, cv::Mat1d &V, cv::Mat1d &U)
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
	auto it = series.cbegin();

	for (int i = 0; i < L; i++)
	{
		std::copy(it, it + k, std::back_inserter(vec));
		++it;
	}

	cv::Mat1d X(L, k, vec.data());


	for (auto it = series.begin(); it != series.end(); ++it)
		std::cout << *it << ",";
	std::cout << std::endl;

	// SVD
	cv::SVD::compute(X, S, U, V);
	cv::Mat1d nonZeroSingularValues = S > 0.0001;
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
	std::cout << "height " << V.size().height << std::endl;
	std::cout << "wodth " << V.size().width << std::endl;
	std::cout << "cols " << V.cols << std::endl;
	std::cout << "rows " << V.rows << std::endl;

}



double CSSA::calc_cmax(const std::vector<double> & series)
{

	cv::Mat1d sing_values, U, V, trend;
	ssa(series, sing_values, V, U);
	int d = std::min(sing_values.rows, int(m_max_et));
	std::cout << V.size() << std::endl;
	//std::cout << U.rowRange(0,5).colRange(0,5) << std::endl;
	//  Calculate LowFreq values for all ETs
	std::vector<double>cval_forET(d);
	for (int et = 0; et < d; et++)
	{
		cval_forET[et] = LFvalue_vect(U, m_omega0, et);
	}
	//---
	int c0s_amount = floor((m_c0max - m_c0min) / m_c0step) + 1;
	for (double c0 = m_c0min; c0 <= m_c0max; c0 += m_c0step)
	{
		std::vector<int> trend_ETs;
		for (int et = 0; et < d; et++)
		{
			if (cval_forET[et] > c0) trend_ETs.push_back(et);
		}
		// reconstruct
		reconstruct(sing_values, U, V, trend_ETs, trend);

		std::vector<double> crit_LFres;

		const cv::Mat1d ts{ series };
		
		double crit = (c0 > 0) ? LF_criteria(ts, trend) : 0;
		crit_LFres.push_back(crit);

	}
	//---
	std::deque<double> c0s{ 0.0 };
	for (double c0 = m_c0min; c0 < m_c0max; c0 += m_c0step)
	{
		c0s.push_back(c0);
	}

	std::valarray<double> tmp1{ std::vector<double>(c0s).data() };
	c0s.
	std::valarray<double> tmp2{ 0.0 };

	tmp(2:length(crit_LFres) + 1) = crit_LFres;
	C0Deltas = [crit_LFres, 0] - tmp;
	C0Deltas = C0Deltas(2:length(C0Deltas) - 1);

	Cmax = 1;
	CmaxDelta = 0;



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


double CSSA::LFvalue_vect(const cv::Mat1d & M, const double omega0, const int ET)
{
	int L = M.cols;

	double nrm = cv::norm(M.col(ET));
	if (nrm == 0.0)
	{
		return -1;
	}

	cv::Mat1d vect(M.col(ET).t() / nrm);
	// calculate sum of periodogram values for low frequencies

	int k = 1;
	double w = 0.0;
	double val = 0;

	while (w < omega0)
	{
		val += periodogram(std::valarray<double>(std::vector<double>(vect).data(), vect.cols), k);
		k = k + 1;
		w = (k - 1.0) / L;
	}

	return val;
}

double CSSA::periodogram(const std::valarray<double>& F, const int k)
{

	int N = F.size();
	double w = (k - 1.0) / N;
	std::valarray<double> T(N);
	for (int i = 0; i < N; ++i) T[i] = (double)i;

	double c = (F * cos(2.0*M_PI*T*w)).sum();
	double s = (F * sin(2.0*M_PI*T*w)).sum();
	double Pi = c*c + s*s;
	Pi = Pi * 2.0 / N;

	if (k == 1 || k == floor(N / 2) + 1) Pi = Pi / 2;
	return Pi;
}

void CSSA::reconstruct(const cv::Mat1d & sing_values, const cv::Mat1d & U, const cv::Mat1d & V, const std::vector<int> ETs, cv::Mat1d & F)
{
	int L = U.cols;
	int K = V.cols;
	int N = K + L - 1;
	cv::Mat1d X{ cv::Mat1d::zeros(L, K) };
	int sz = ETs.size();
	for (int i = 0; i < sz; i++)
	{
		int ET = ETs[i];
		const double sing_value{ sing_values(ET,0) };
		const cv::Mat1d u{ U.col(ET) };
		const cv::Mat1d v{ V.row(ET) };
		const cv::Mat1d Xi{ sing_value * u * v };

		X += Xi;
	}
	//	std::cout << "X:" << std::endl;
	//	std::cout << X << std::endl;
	cv::Mat1d Freq(1, N);
	double L0 = std::min(L, K);
	double K0 = std::max(L, K);

	cv::Mat1d Y = (L<K) ? X : X.t();


	for (int k = 1; k < L0; k++)
	{
		int d = 0;
		for (int m = 1; m <= k; m++)
		{
			d += Y(m, k - m + 1);
		}
		Freq(k - 1) = 1.0 / k * d;
	}
	for (int k = L0; k <= K0; k++)
	{

		int d = 0;
		for (int m = 1; m <= L0; m++)
		{
			d += Y(m, k - m + 1);
		}
		Freq(k - 1) = 1.0 / k * d;
	}
	for (int k = K0 + 1; k <= N; k++)
	{
		int d = 0;
		for (int m = k - K0; m <= N - K - +1; m++)
		{
			d += Y(m, k - m + 1);
		}
		Freq(k - 1) = 1.0 / (N - k + 1) * d;
	}
	F = std::move(Freq);
	std::cout << F;
}

double CSSA::LF_criteria(const cv::Mat1d & ts, const cv::Mat1d & trend)
{

	double lf_ts = LFvalue_vect(ts, m_omega0, 0);
	double lf_res = LFvalue_vect((ts - trend), m_omega0, 0);
	if (lf_ts == -1 || lf_res == -1) return 0.0;

	return (lf_res / lf_ts);
}
