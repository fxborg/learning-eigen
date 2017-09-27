
#include <Windows.h>
#include <tchar.h>
#include <iostream>

#include <fstream>
#include <string>
#include <vector>


typedef int(__stdcall *Fn_ExpClass_Create)(int, int, double ,int,  double);
typedef int(__stdcall *Fn_ExpClass_Push)(int, int, double, time_t, time_t);
typedef int(__stdcall *Fn_ExpClass_Calc)(int);
typedef void(__stdcall *Fn_ExpClass_Destroy)(int);
int main()
{
	HMODULE hDll = ::LoadLibrary(_T("ssa.dll"));

	Fn_ExpClass_Create ExpClass_Create = (Fn_ExpClass_Create)GetProcAddress(hDll, "Create");
	Fn_ExpClass_Destroy ExpClass_Destroy = (Fn_ExpClass_Destroy)GetProcAddress(hDll, "Destroy");
	Fn_ExpClass_Push ExpClass_Push = (Fn_ExpClass_Push)GetProcAddress(hDll, "Push");
	Fn_ExpClass_Calc ExpClass_Calc = (Fn_ExpClass_Calc)GetProcAddress(hDll, "Calculate");

	int ptr = ExpClass_Create(100, 48, 0.001,24, 0.01);
	time_t t0 = 1000;
	time_t t1 = 1000;
	t1 = t0;

	std::string src_path="C:\\Users\\fxb\\Desktop\\AutoSSA-matlab_22102007\\data\\USDJPY4H_2.dat";
	std::ifstream fin(src_path);
	if (fin.bad()) {
		std::cout << "File Open Error:" << src_path << std::endl;
		return 1;
	}
	std::cout << "File Open!:" << src_path << std::endl;
	std::string str;
	int n = 0;
	double price;
	std::vector<std::vector<double>>sample;
	while (std::getline(fin, str))
	{
		int res = sscanf_s(str.c_str(), "%lf", &price);

		if (res != 1)continue;
		n++;
		t1 = t0++;
		ExpClass_Push(ptr, n, price, t0, t1);
	}
	fin.close();
	n++;
	t1 = t0++;
	ExpClass_Push(ptr, n, (double)100.0, t0, t1);

	//
	ExpClass_Calc(ptr);
	ExpClass_Destroy(ptr);

	//dllで作成したクラスインスタンスを作成する
	system("pause");
	return 0;
}

