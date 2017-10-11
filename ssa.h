//+------------------------------------------------------------------+
//|                                                          ssa.h   |
//| SSA Trend                                 Copyright 2017, fxborg |
//| this is a c++ reimplementation of AutoSSA Matlab package         |
//| AutoSSA(http://www.pdmi.ras.ru/~theo/autossa/english/soft.htm)   |
//|                                   http://fxborg-labo.hateblo.jp/ |
//+------------------------------------------------------------------+
#define EXPORT  extern "C" // __declspec(dllexport)  
#define _USE_MATH_DEFINES
#include <deque> 
#include <vector> 
#include <unordered_map> 
#include <valarray> 
#include <numeric>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <chrono>
#include "series.h"
#include "stats.h"
#include "opencv2/core/core.hpp"
#include "RedSVD-h"
class CSSA
{
public:
	bool get_results(const size_t idx, double & y);
	// コンストラクタ
	explicit CSSA(const size_t size, const size_t length, const  double omega0, const size_t max_et, const double c0step);
	int push(const int x, const double y, const time_t t0, const time_t t1);
	int calculate();
	double slope(const int period);


private:
	std::vector<int> CSSA::calc_trend(const std::vector<double> & series, const double c0, const cv::Mat1d &S, const cv::Mat1d &V, const cv::Mat1d &U);

	double calc_cmax(const std::vector<double> & series, const cv::Mat1d & sing_values, const cv::Mat1d & V, const cv::Mat1d &U);
	void ssa(const std::vector<double> & series, cv::Mat1d &S, cv::Mat1d &V, cv::Mat1d &U);
	size_t minmax(const size_t n, const size_t min, const size_t max);
	double minmax(const double n, const double min, const double max);
	double periodogram(const std::valarray<double> & F, const int k);
	double LFvalue_vect(const cv::Mat1d & M, const double omega0, const int ET);

	cv::Mat1d reconstruct(const cv::Mat1d & sing_values, const cv::Mat1d & U, const cv::Mat1d & V, const std::vector<int> & ETs);
	std::string join_ETs(const std::vector<int>& ETs);
	std::vector<int> certainly_ETs();


	const size_t m_size;
	const size_t m_L;
	const double m_omega0;
	const size_t m_max_et;
	const double m_c0min;
	const double m_c0max;
	const double m_c0step;
	const double m_rdelta_threshold;
	const double m_c0eps;
	const size_t m_hist_size;
	const double m_hist_coef;
	CSeries m_series;
	std::vector<double>m_results;
	std::vector<double> m_cval_for_ET;
	std::unordered_map<std::string, double> m_memo_criteria;
	std::deque<std::vector<int>> m_hist_ETs;
	std::vector<std::chrono::time_point<std::chrono::system_clock>> timer;

};

//--- インスタンスを生成
EXPORT CSSA * __stdcall Create(const size_t size, const size_t length, const  double omega0, const size_t max_et, const double c0step);
//--- インスタンスを破棄
EXPORT void __stdcall Destroy(CSSA* instance);
//--- インスタンス経由でpushメソッドをコール
EXPORT int __stdcall Push(CSSA* instance, const int x, const double y, const time_t t0, const time_t t1);
//--- SSAを計算
EXPORT int __stdcall Calculate(CSSA* instance);


EXPORT bool __stdcall GetResults(CSSA* instance, const size_t idx, double &y);

EXPORT double __stdcall Slope(CSSA* instance, const int period);

