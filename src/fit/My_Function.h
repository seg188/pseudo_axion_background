#include <iostream>
#ifndef MYFUNCTION_H
#define MYFUNCTION_H



class My_Function {
	
public:

My_Function();

static TF1* TFit_Function;
static void Fit_Function(TF1* func) { TFit_Function = func; }

static std::vector<TFitResultPtr> Fit_Results;
static void Set_Fit_Results(std::vector<TFitResultPtr> res){Fit_Results = res; }

static std::vector<std::vector<double>> pt_params;
static void Set_pt_params(std::vector<std::vector<double>> val){ pt_params = val;}

static std::vector<double> bin_centers;
static void Set_bin_centers(std::vector<double> fbin_centers){bin_centers = fbin_centers;}

static double math_pow(double x, double pow); //pow > 0
static double Pnomial(std::vector<double> params, double x);
static double background_pt_function(double x, double mass_omg);

static std::vector<TFitResultPtr> First_Fit(std::vector<TH1*>, bool, double, double);

static std::vector<TH1*> data_histo_vector;
static void Set_data_histo_vector(std::vector<TH1*> vec){ data_histo_vector = vec; }

static bool s_opt_draw;
static void Set_opt_draw(bool opt){ s_opt_draw = opt; }

static double s_fit_min;
static double s_fit_max;
static void Set_Range(double min, double max){ s_fit_min = min; s_fit_max = max; }

static TString Pnomial_String(int order);

static double LagrangePol(std::vector<double> x_vals, std::vector<double> y_vals, double x);

static std::vector<std::vector<double>> Get_Fit_Params(std::vector<TFitResultPtr>);
static std::vector<double> GetBinCenters();

static void SetSecondaryFunction(){ Set_bin_centers(GetBinCenters()) ; Set_pt_params(Get_Fit_Params(First_Fit(data_histo_vector, s_opt_draw, s_fit_min, s_fit_max))); }
static void Retry_Fit();

static TFitResultPtr Fit_Fixed(int, int);
static TFitResultPtr Fit_Fixed2(int, int);

static void FF_All();

static void New_Fit();

static double GetMaxBin(TH1*);

static void HistogramizeParameters();

static TF1 *TFit_First_Function;
static TF1 *TFit_Second_Function;
static TF1Convolution *f_conv;

static std::vector<TFitResultPtr> pol_fit(std::vector<TH1*>, TString);



};

#endif