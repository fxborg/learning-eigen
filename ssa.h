#define EXPORT  extern "C" // __declspec(dllexport)  
#define _USE_MATH_DEFINES
#include <deque> 
#include <vector> 
#include <valarray> 
#include <numeric>
#include <algorithm>
#include <iterator>
#include <math.h>
#include "series.h"
#include "opencv2/core/core.hpp"
class CSSA
{
public:
	// コンストラクタ
	explicit CSSA(const unsigned int size, const unsigned int length, const  double omega0, const unsigned int max_et, const double c0step);
	int push(const int x, const double y, const time_t t0, const time_t t1);
	int calculate();



private:
	double calc_cmax(const std::vector<double> & series);
	void ssa(const std::vector<double> & series, cv::Mat1d &S, cv::Mat1d &V, cv::Mat1d &U);
	unsigned int minmax(const unsigned int n, const unsigned int min, const unsigned int max);
	double minmax(const double n, const double min, const double max);
	double periodogram(const std::valarray<double> & F, const int k);
	double LFvalue_vect(const cv::Mat1d & M, const double omega0, const int ET);
	void reconstruct(const cv::Mat1d & sing_values, const cv::Mat1d & U, const cv::Mat1d & V, const std::vector<int> ET, cv::Mat1d & F);
	double LF_criteria(const cv::Mat1d & ts, const  cv::Mat1d & trend);
	const unsigned int m_size;
	const unsigned int m_L;
	const double m_omega0;
	const unsigned int m_max_et;
	const double m_c0min;
	const double m_c0max;
	const double m_c0step;
	const double m_rdelta_threshold;
	const double m_c0eps;
	CSeries m_series;
};

//--- インスタンスを生成
EXPORT CSSA * __stdcall Create(const unsigned int size, const unsigned int length, const  double omega0, const unsigned int max_et, const double c0step);
//--- インスタンスを破棄
EXPORT void __stdcall Destroy(CSSA* instance);
//--- インスタンス経由でpushメソッドをコール
EXPORT int __stdcall Push(CSSA* instance, const int x, const double y, const time_t t0, const time_t t1);
//--- 予測値を計算
EXPORT int __stdcall Calculate(CSSA* instance);
