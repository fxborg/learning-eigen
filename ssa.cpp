//+------------------------------------------------------------------+
//|                                                          ssa.cpp |
//| SSA Trend                                 Copyright 2017, fxborg |
//| this is a c++ reimplementation of AutoSSA Matlab package         |
//| AutoSSA(http://www.pdmi.ras.ru/~theo/autossa/english/soft.htm)   |
//|                                   http://fxborg-labo.hateblo.jp/ |
//+------------------------------------------------------------------+

#include "ssa.h"
EXPORT CSSA * __stdcall Create(const size_t size, const size_t length, const  double omega0, const size_t max_et, const double c0step) {

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
EXPORT double __stdcall Slope(CSSA* instance, const int period)
{
	return instance->slope(period);
}


EXPORT bool __stdcall GetResults(CSSA * instance, const size_t idx, double &y)
{
	return instance->get_results(idx, y);
}

bool CSSA::get_results(const size_t idx, double &y)
{
	if (m_results.size() <= idx) return false;
	y = m_results[idx];
	return true;
}

CSSA::CSSA(const size_t size, const size_t length, const  double omega0, const size_t max_et, const double c0step) :
	m_size(minmax(size, 50, 2000)),
	m_L(minmax(length, 5, m_size)),
	m_omega0(minmax(omega0, 0.001, 0.5)),
	m_max_et(minmax(max_et, 2, int(0.5*m_L))),
	m_c0min(0.5),
	m_c0max(1.0),
	m_c0step(c0step),
	m_rdelta_threshold(0.05),
	m_c0eps(0),
	m_hist_size(20),
	m_hist_coef(0.8),
	m_series(CSeries(m_size + 1)),
	m_hist_ETs({})
{
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
	if (series.size() != m_size + 1)return 0;
	// 結果を初期化
	m_results.erase(m_results.begin(), m_results.end());
	// メモを初期化
	m_cval_for_ET.erase(m_cval_for_ET.begin(), m_cval_for_ET.end());
	m_memo_criteria.erase(m_memo_criteria.begin(), m_memo_criteria.end());



	const std::vector<double> ts(series.cbegin(), series.cend() - 1);




	cv::Mat1d sing_values, V, U;

	timer.emplace_back(std::chrono::system_clock::now());
	//--- SSA
	ssa(ts, sing_values, V, U);
	//--- cmax

	timer.emplace_back(std::chrono::system_clock::now());
	double c0maxval = calc_cmax(ts, sing_values, V, U);
	timer.emplace_back(std::chrono::system_clock::now());
	for (int i = 1; i < timer.size(); i++)
	{
		std::cout << "time" << (i - 1) << " -> time" << i << " = "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(timer[i] - timer[i - 1]).count()
			<< " msec."
			<< std::endl;
	}


	cv::Mat1d s, u, v;

	if (sing_values.rows > int(m_max_et))
	{

		s = sing_values.rowRange(0, m_max_et);
		u = U.colRange(0, m_max_et);
		v = V.colRange(0, m_max_et);

	}
	else
	{
		s = sing_values;
		u = U;
		v = V;

	}

	// calc trend
	const std::vector<int>trend_ETs{ calc_trend(ts, c0maxval, s, v, u) };

	m_hist_ETs.push_back(std::move(trend_ETs));
	// history
	if (m_hist_ETs.size() > m_hist_size)
	{

		size_t delsz = m_hist_ETs.size() - m_hist_size;
		m_hist_ETs.erase(m_hist_ETs.begin(), m_hist_ETs.begin() + (delsz));

	}

	std::vector<int> certETs{ certainly_ETs() };

	//for (int i = m_hist_ETs.size()-1; i >= 0; i--)
	//{
	//	if (m_hist_ETs[i].size()>0)
	//		std::cout << i <<":"<< join_ETs(m_hist_ETs[i]) << std::endl;
	//}
	//std::cout << "CertETs : " << join_ETs(certETs) << std::endl;

	const std::vector<double> res{ reconstruct(sing_values, U, V, certETs) };
	m_results = std::move(res);
	timer.emplace_back(std::chrono::system_clock::now());

	return m_results.size();
}



std::vector<int> CSSA::calc_trend(const std::vector<double> & series, const double c0, const cv::Mat1d &sing_values, const cv::Mat1d &V, const cv::Mat1d &U)
{
	int maxET = std::min((int)m_cval_for_ET.size(), U.cols);

	std::vector<int> trend_ETs;

	for (int et = 0; et < maxET; et++)
	{
		if (m_cval_for_ET[et] > c0)trend_ETs.push_back(et);
	}
	return trend_ETs;
}

double CSSA::slope(const int period)
{
	int sz = m_results.size();
	if (sz < period) return 0.0;

	Stats stat = Stats();
	for (int i = sz - period; i < sz; i++)
	{
		stat += Stats(1, i, m_results[i]);
	}
	return stat.slope();
}

void CSSA::ssa(const std::vector<double> & series, cv::Mat1d &S, cv::Mat1d &V, cv::Mat1d &U)
{

	int k = m_size - m_L + 1;
	int L = int(m_L);

	timer.emplace_back(std::chrono::system_clock::now());
	std::vector<double> vec;
	vec.reserve(L*k);
	int adj = series.size() - m_size;
	auto it = series.cbegin();

	for (int i = 0; i < L; i++)
	{
		std::copy(it, it + k, std::back_inserter(vec));
		++it;
	}



	timer.emplace_back(std::chrono::system_clock::now());


	//	typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

	Eigen::MatrixXd A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(&vec[0], L, k);

	cv::Mat1d X(L, k, vec.data());




	// opencv SVD
	timer.emplace_back(std::chrono::system_clock::now());
	cv::SVD::compute(X, S, U, V);
	timer.emplace_back(std::chrono::system_clock::now());

	// eigen3 SVD 
	timer.emplace_back(std::chrono::system_clock::now());
	//	JacobiSVD< MatrixXd> svd(A, ComputeThinU | ComputeThinV);	timer.emplace_back(std::chrono::system_clock::now());


	// RedSVD
	timer.emplace_back(std::chrono::system_clock::now());
	RedSVD::RedSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	timer.emplace_back(std::chrono::system_clock::now());
	cv::Mat1d U2,S2,V2;
	cv::eigen2cv(svd.matrixU(), U2);
	cv::eigen2cv(svd.singularValues(), S2);
	cv::eigen2cv(svd.matrixV(), V2);



	cv::Mat1d nonZeroSingularValues = S > 0.0001;


	// rank
	int M = cv::countNonZero(nonZeroSingularValues);
	if (M < L)
	{
		S = S.rowRange(0, M);
		U = U.colRange(0, M);
	}

	V = V.rowRange(0, M).t();


}


double CSSA::calc_cmax(const std::vector<double> & series, const cv::Mat1d & sing_values, const cv::Mat1d & V, const cv::Mat1d &U)
{
	//	timer.emplace_back(std::chrono::system_clock::now());

	cv::Mat1d trend;
	int d = std::min(sing_values.rows, int(m_max_et));
	//  Calculate LowFreq values for all ETs
	for (int et = 0; et < d; et++)
	{
		m_cval_for_ET.push_back(LFvalue_vect(U, m_omega0, et));
	}

	//---
	int c0s_amount = (int)floor((m_c0max - m_c0min) / m_c0step) + 1;


	const cv::Mat1d ts{ series };
	double lf_ts = LFvalue_vect(ts, m_omega0, 0);

	std::vector<double> crit_LFres;
	std::vector<double> c0s;

	for (int n = 0; n < c0s_amount; n++)
	{
		double c0 = m_c0min + m_c0step*n;
		std::vector<int> trend_ETs;
		for (int et = 0; et < d; et++)
		{
			if (m_cval_for_ET[et] > c0)
			{

				trend_ETs.push_back(et);
			}
		}


		// reconstruct

		//--- test cache
		const std::string key{ join_ETs(trend_ETs) };

		double crit;
		if (m_memo_criteria.count(key) == 0)
		{
			//---
			const cv::Mat1d trend{ reconstruct(sing_values, U, V, trend_ETs) };

			double lf_res = LFvalue_vect((ts - trend), m_omega0, 0);
			crit = (lf_ts > 0 && lf_res > 0) ? (lf_res / lf_ts) : 0.0;
			m_memo_criteria[key] = crit;
			//---

		}
		else
		{
			crit = m_memo_criteria[key];
		}
		c0s.push_back(c0);
		crit_LFres.push_back(crit);

	}
	//---
	int sz = crit_LFres.size();

	//int n = 0;
	//for (auto & v : crit_LFres)
	//{
	//	std::cout << n++ << ":" << v << std::endl;
	//}

	const std::valarray<double> tmp(crit_LFres.data(), sz);
	const std::valarray<double> c0dtl = tmp[std::slice(1, sz - 1, 1)] - tmp[std::slice(0, sz - 1, 1)];



	double cmax = 1.0;
	for (int i = 0; i < (int)c0dtl.size(); i++)
	{
		if (c0dtl[i] > m_rdelta_threshold)
		{
			cmax = c0s[i];
			break;
		}
	}


	//	timer.emplace_back(std::chrono::system_clock::now());
	return cmax;
}

size_t CSSA::minmax(const size_t n, const size_t min, const size_t max)
{
	return std::max(min, std::min(max, n));
}

double CSSA::minmax(const double n, const double min, const double max)
{
	return std::max(min, std::min(max, n));
}


double CSSA::LFvalue_vect(const cv::Mat1d & M, const double omega0, const int ET)
{
	int L = M.rows;


	double nrm = cv::norm(M.col(ET), cv::NORM_L2);
	if (nrm == 0.0)
	{
		return -1;
	}

	cv::Mat1d vect(M.col(ET).t() / nrm);
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

cv::Mat1d CSSA::reconstruct(const cv::Mat1d & sing_values, const cv::Mat1d & U, const cv::Mat1d & V, const std::vector<int> & ETs)
{
	int L = U.rows;
	int K = V.rows;
	int N = K + L - 1;
	cv::Mat1d X{ cv::Mat1d::zeros(L, K) };
	for (size_t i = 0; i < ETs.size(); i++)
	{
		int ET = ETs[i];
		const double sing_value{ sing_values(ET,0) };
		const cv::Mat1d u{ U.col(ET) };
		const cv::Mat1d v{ V.col(ET) };
		const cv::Mat1d Xi{ sing_value * u * v.t() };

		X += Xi;
	}
	cv::Mat1d Freq{ cv::Mat1d::zeros(N, 1) };
	int L0 = std::min(L, K);
	int K0 = std::max(L, K);
	cv::Mat1d Y;
	if (L<K)
		Y = X;
	else
		Y = X.t();


	//--- 斜めに加算 前半
	for (int k = 1; k < L0; k++)
	{
		double d = 0.0;
		for (int m = 1; m <= k; m++)
		{
			d += Y(m - 1, k - m);
		}
		Freq(k - 1, 0) = 1.0 / k * d;
	}

	//---
	for (int k = L0; k <= K0; k++)
	{

		double d = 0.0;
		for (int m = 1; m <= L0; m++)
		{
			d += Y(m - 1, k - m);
		}
		Freq(k - 1, 0) = 1.0 / L0 * d;
	}
	//---
	for (int k = K0 + 1; k <= N; k++)
	{
		double d = 0.0;
		for (int m = k - K0 + 1; m <= N - K0 + 1; m++)
		{
			d += Y(m - 1, k - m);
		}

		Freq(k - 1, 0) = 1.0 / (N - k + 1) * d;
	}
	return  Freq;
}


std::string CSSA::join_ETs(const std::vector<int> & ETs)
{
	if (ETs.empty()) return "";

	std::string s = std::to_string(ETs[0]);
	for (size_t i = 1; i < ETs.size(); ++i) s += "," + std::to_string(ETs[i]);

	return s;
}

std::vector<int> CSSA::certainly_ETs()
{
	std::map<int, int> work;
	int sz = m_hist_ETs.size();
	if (sz == 0) return {};

	for (auto &vec : m_hist_ETs)
	{
		for (auto &et : vec)
		{
			if (work.count(et) == 0)	work[et] = 1;
			else						work[et]++;
		}
	}
	std::vector<int> res;
	//	std::cout << "m_hist_coef=" << (m_hist_coef) << std::endl;
	for (auto &kv : work)
	{
		//		std::cout << kv.first << "->" << kv.second;
		//		std::cout << " = " << ((double)kv.second / sz);
		if (((double)kv.second / sz) >= m_hist_coef)
		{
			//			std::cout << "*";
			res.push_back(kv.first);
		}
		//		std::cout << std::endl;
	}
	return res;
}