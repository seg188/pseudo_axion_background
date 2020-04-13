
#include <iostream>
#include "My_Function.h"

#ifndef MYFUNCTION_CLASS
#define MYFUNCTION_CLASS


std::vector<TH1*> My_Function::data_histo_vector({0});
bool My_Function::s_opt_draw(true);
double My_Function::s_fit_min(0.3);
double My_Function::s_fit_max(5.0);
std::vector<double> My_Function::bin_centers(0);
std::vector<TFitResultPtr> My_Function::Fit_Results(0);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void My_Function::Retry_Fit()
{
	std::vector<std::vector<double>> good_params;
	std::vector<double> good_params_centers;
	int l = 0;


	for (TFitResultPtr res : My_Function::Fit_Results)
	{
		if (res->IsValid())
		{
			good_params.push_back(res->Parameters());
			good_params_centers.push_back(My_Function::bin_centers[l++]);
		}
		else {l++;}
	}

	std::vector<std::vector<double>> l_params;
	l_params.resize(good_params[0].size());

	for (int i = 0; i < good_params[0].size(); i++)
	{
		for (int j = 0; j < good_params.size(); j++)
		{
			//std::cout << "beep boop" << std::endl;
			l_params[i].push_back(good_params[j][i]);
		}
	}

	std::vector<TFitResultPtr> new_results;

	for (int k = 0; k < My_Function::Fit_Results.size(); k++)
	{
		if (! My_Function::Fit_Results[k]->IsValid())
		{
			for (int i = 0; i < l_params[0].size(); i++)
			{
				TFit_Function->SetParameter(i, My_Function::LagrangePol(good_params_centers, good_params[i], My_Function::bin_centers[k]));
			}

			double new_fit_max = My_Function::GetMaxBin(data_histo_vector[k]);

			new_results.push_back(data_histo_vector[k]->Fit("TFit_Function", "S", "M", s_fit_min, new_fit_max));
		}
		else {new_results.push_back(My_Function::Fit_Results[k]);}
	}

	My_Function::Set_Fit_Results(new_results);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TF1* My_Function::TFit_First_Function = new TF1("my_gaussian", "[0] * exp(-1 / [1] * TMath::Power(x - [2], 2))", 0.3, 4.0);
//TF1* My_Function::TFit_Second_Function = new TF1("my_laundau", "TMath::Landau([0] * x - [1]  )", 0.3, 4.0);
//TF1Convolution *My_Function::f_conv = new TF1Convolution(My_Function::TFit_First_Function, My_Function::TFit_Second_Function);
//TF1* My_Function::TFit_Function(new TF1("TFit_Function", *My_Function::f_conv , 0.5, 4.0, My_Function::f_conv->GetNpar()));

TF1* My_Function::TFit_Function(new TF1("TFit_Function", "[0] * exp(-1 / [3] * TMath::Power(x - [4], 2)) * TMath::Landau([1] * x - [2])" , 0.5, 4.5));
//TF1* My_Function::TFit_Function(new TF1("TFit_Function", "[0] * exp([1] * (log(x - [3]) - [2] * x - [4] * TMath::Power(x, 2)))" , 0.5, 4.5));
//
// * ) * exp(-1*[5]*x)"
//* 

std::vector<std::vector<double>> My_Function::pt_params{{0}};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double My_Function::math_pow(double x, double pow) //pow > 0
{
   double val = x;
   for (int i = 1; i < pow; i++){val = val * x;}
   return val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double My_Function::background_pt_function(double x, double mass_omg) 
{
   double val = 1;

   val = val * My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[0], x);
   val = val * TMath::Landau(My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[1], x) * mass_omg - My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[2], x));
   val = val * exp(-1 / My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[3], x) * My_Function::math_pow( mass_omg - My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[4], x), 2));
  // val = val * exp(-1 * My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[5], x) * mass_omg);

   return val;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//using LagrangePol function of form [0] + [1] * x + [2] * X^2 + [3] * x^3 ... up to [6] or above
//THE GOAL OF THIS FUNCTION IS TO RETURN SOMETHING OF THE FOLLOWING FORM: (stored in pt_params)
//{ [O] , [1] , [2] , [3] ...} FOR FIRST PARAMETER IN DATA SET 
//{ [O] , [1] , [2] , [3] ...} FOR SECOND PARAMETER IN DATA SET ... ETC..... FOR PARAMETER [0], INDEX BY pt_params[0][x]

std::vector<TFitResultPtr> My_Function::First_Fit(std::vector<TH1*> histo_vector, bool opt_draw = true, double min = 0.5, double max = 4.)
{

	std::vector<TFitResultPtr> ffit_results;
	std::vector<TCanvas*> canvas_vector;
	std::vector<std::vector<double>> fit_params;
	int j = 0; 
	double fit_min = min; double fit_max = max;

	///////////////////////////////////////////////////////////////////////////////
	///FINDING MAX/MIN OF RANGE TO USE IN FUNCTION, OTHERWISE USING DEFAULT MAX, MIN

  	for (TH1* histo : histo_vector){ 

  		for (int i = 1; i < histo->GetNbinsX(); i++) { 
  			if (histo->GetBinContent(i -1) < 0.001 and histo->GetBinContent(i) < 0.001 and histo->GetBinContent(i + 1) > 0.001 and i < 10){
 	   			fit_min = histo->GetBinCenter(i);
    			if (histo->GetBinContent(i+1) > 0.01) { break; } 
    		//	else if (histo->GetBinContent(i -1) > 0.001 and histo->GetBinContent(i) > 0.001 and  histo->GetBinContent(i + 1) < 0.001) {fit_max = histo->GetBinCenter(i);}
    			} } 


    canvas_vector.push_back(new TCanvas(TString(histo->GetName()) + TString("_canvas"), histo->GetName()));
    My_Function::TFit_Function->SetParameters(1, 3, 4, 1, -3, -1);
    //My_Function::TFit_Function->SetParameters(3, 0, 0, 10, 10);
    My_Function::TFit_Function->SetParLimits(0, 0, 2);
    My_Function::TFit_Function->SetParLimits(2, 3., 6.);
    My_Function::TFit_Function->SetParLimits(3, 0., 6.);
   // My_Function::TFit_Function->FixParameter(4, 1.25);
   // My_Function::TFit_Function->FixParameter(2, 4.);
     //  My_Function::TFit_Function->SetParLimits(4, -10., 10.);
     //  My_Function::TFit_Function->SetParLimits(2, -10., 10.);
     //	  My_Function::TFit_Function->SetParLimits(5, 0, 100.);

    
    ffit_results.push_back(histo->Fit("TFit_Function", "S", "M", fit_min, fit_max));
}


My_Function::Set_Fit_Results(ffit_results);

My_Function::New_Fit();

//My_Function::New_Fit();



/*
 for (int i = 0; i < fit_results.size(); i++)
 {
 	if (!fit_results[i]->IsValid())
 	{
 		for (int j = 0; j < fit_results[i]->Parameters().size(); j++ ) {My_Function::TFit_Function->SetParameter(j, fit_results[j]->Parameter(j));}

 		fit_results[i] = histo_vector[i]->Fit("TFit_Function", "S", "S", fit_min, fit_max);
 	}
 	
 }
*/


////IN THIS SECTION, WE WILL CHECK IF A FIT IS VALID. IF IT IS, WE WILL KEEP IT. IF NOT, WE WILL DISCARD IT, AND USE THE OTHER FIT
//VALUES TO GUESS THE VALUE OF THE PARAMETERS TO DO A NEW FIT/////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
for (int times = 0; times < 3; times++)
{
	std::cout << "\n\n\n  IMPROVED FIT RESULTS *********************************************** \n\n" << std::endl;
	My_Function::Set_Fit_Results(ffit_results);
	My_Function::Retry_Fit(fit_min, fit_max);


}
*/



for (TH1* histo : histo_vector){ 
	if (opt_draw) {
	//	canvas_vector[j]->cd(); 
//		histo->Draw();

		//canvas_vector[j]->Print(TString("../../Plots/") + TString(canvas_vector[j]->GetName()) + TString(".pdf"), ".pdf" );
		j++;
	}
}


return My_Function::Fit_Results; 

}


double My_Function::LagrangePol(std::vector<double> real_x, std::vector<double> real_y, double x)
{
	double total_val = 0;

	std::vector<double> x_vals = real_x;
	std::vector<double> y_vals = real_y;

	for (int ii =0; ii < real_x.size() - 1; ii++){

		x_vals.push_back((x_vals[ii] + x_vals[ii + 1]) / 2 );
		y_vals.push_back((y_vals[ii] + y_vals[ii + 1]) / 2 );
	}


	for (int i = 0; i < x_vals.size(); i++)
	{
		double val = y_vals[i] * TMath::Power(x_vals[i] / x, 2);

		for (int j = 0; j < x_vals.size(); j++) {
			if (i != j) {

				val = val * (x - x_vals[j]) / (x_vals[i] - x_vals[j]);
			}
		}

		
		total_val += val;
		
		
	}

	return total_val;
}

std::vector<std::vector<double>> My_Function::Get_Fit_Params(std::vector<TFitResultPtr> results_v)
{

	std::vector<std::vector<double>> og_params, pol_params;

	for (TFitResultPtr res : results_v){ og_params.push_back(res->Parameters()); }

	pol_params.resize(og_params[0].size());

	for (int i = 0; i < og_params[0].size(); i++){ 
		for (int j = 0; j < og_params.size(); j++) { 

			pol_params[i].push_back(og_params[j][i]);
		}
	}

	return pol_params;
}

std::vector<double> My_Function::GetBinCenters()
{
	std::vector<double> v;

	v.push_back(30);
	v.push_back(50);
	v.push_back(70);
	v.push_back(100);
	v.push_back(140);
	v.push_back(200);
	v.push_back(270);
	v.push_back(380);
	v.push_back(580);

	return v;
}

TFitResultPtr My_Function::Fit_Fixed(int n, int func_n){

	for (int k = 0; k < n; k++){My_Function::TFit_Function->ReleaseParameter(k);}

	if (n < My_Function::TFit_Function->GetNpar()) {My_Function::TFit_Function->FixParameter(n, (My_Function::Fit_Results[func_n]->Parameters())[n]);}
	else { 
		double new_maxx = My_Function::GetMaxBin(data_histo_vector[func_n]);
	//	My_Function::TFit_Function->SetParLimits(0, 0, 2);
	//	My_Function::TFit_Function->SetParLimits(3, -2.5, 3.);
   //		 My_Function::TFit_Function->SetParLimits(4, -2., 3.);
		return data_histo_vector[func_n]->Fit("TFit_Function", "S", "M", My_Function::s_fit_min, new_maxx); }
	
	return My_Function::Fit_Fixed(n + 1, func_n);

}

TFitResultPtr My_Function::Fit_Fixed2(int n, int func_n){

	for (int k = 0; k < n; k++){My_Function::TFit_Function->ReleaseParameter(k);}

	if (n < My_Function::TFit_Function->GetNpar() and n >= 0) {My_Function::TFit_Function->FixParameter(n, (My_Function::Fit_Results[func_n]->Parameters())[n]);}
	else { 
		double new_maxx = My_Function::GetMaxBin(data_histo_vector[func_n]);
//		My_Function::TFit_Function->SetParLimits(0, 0, 2);
//		My_Function::TFit_Function->SetParLimits(3, 0, 5.);
//		    My_Function::TFit_Function->SetParLimits(3, -2.5, 3.);
    My_Function::TFit_Function->SetParLimits(4, -2., 3.);
		return data_histo_vector[func_n]->Fit("TFit_Function", "S", "M", My_Function::s_fit_min, new_maxx); }
	
	return My_Function::Fit_Fixed(n - 1, func_n);

}

void My_Function::FF_All(){ //can only be called after first fit or within first fit

	for (int jj = 0; jj < My_Function::Fit_Results.size(); jj++){

		if (My_Function::Fit_Results[jj]->IsValid() == 1){

		} else {

			My_Function::Fit_Results[jj] = My_Function::Fit_Fixed(0, jj);
			My_Function::Fit_Results[jj] = My_Function::Fit_Fixed2(My_Function::Fit_Results.size() - 1, jj);
			My_Function::Fit_Results[jj] = My_Function::Fit_Fixed(0, jj);
		}

	}



}

void My_Function::New_Fit(){

	//for (int jj = 0; jj < 5; jj++){

		My_Function::FF_All();
	//}


	for (TFitResultPtr res : My_Function::Fit_Results)
	{

	for (int j = 0; j < 3; j++){
		if (res->IsValid() == 0 or (res->Parameters())[0] < 0 )
		{
		//	std::cout << "Reinterpolating Parameters..........................................................................\n\n" << std::endl;
		//	My_Function::Retry_Fit();

		//	for (int jj = 0; jj < 3; jj++){

				std::cout << "\n\n\n  IMPROVED FIT RESULTS *********************************************** \n\n" << std::endl;

				My_Function::FF_All();
		//	}

			//break;	

		}
		}
	}


	
}

double My_Function::GetMaxBin(TH1* histo){
	double maxx = My_Function::s_fit_max;

	for (int jj = 1; jj <= histo->GetNbinsX(); jj++){
		if (histo->GetBinContent(jj) == 0 and histo->GetBinContent(jj - 1) != 0 and histo->GetBinContent(jj + 1) == 0)
		{
			maxx = histo->GetBinCenter(jj);
		}
	}

	return maxx;
}


TString getName(int j)
{
  std::ostringstream strs_name;
  strs_name << j;

  return ( TString("PT_") + TString(strs_name.str()) );
}


void My_Function::HistogramizeParameters(){ //CALL THIS FUNCTION TO DRAW THE PARAMETER VALUES AND THE FUNCTIONS DESIGNED TO FIT HTEM

double bin_edges[10]{20, 40, 60, 80, 120, 160, 240, 300, 460, 700};
std::vector<TH1*> histogram_vector;
int j = 0;
for (double par : My_Function::Fit_Results[0]->Parameters()){

	std::ostringstream strs;
	strs << j++;
	histogram_vector.push_back(new TH1D(TString("par_") + TString(strs.str()), TString("par_") + TString(strs.str()), My_Function::bin_centers.size(), bin_edges ));

}


for (int ii = 0; ii < My_Function::pt_params.size(); ii++){

	for ( int jj = 0; jj < My_Function::pt_params[ii].size() ; jj++){
		
		histogram_vector[ii]->SetBinContent(jj + 1, My_Function::pt_params[ii][jj]);
	}
}

std::vector<TF1*> func_vector;
std::vector<TCanvas*> canvas_v;

std::vector<double> x1 = My_Function::bin_centers;

j = 0;

for (TH1* histo : histogram_vector){

	std::ostringstream strs;
	strs << j;
	func_vector.push_back(new TF1(histo->GetName() + TString("_function"), TString("My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[") + TString(strs.str()) + TString("], x)"), 20, 500 ));

	canvas_v.push_back(new TCanvas(TString("canvas_par_") + TString(strs.str()), TString("canvas_par_") + TString(strs.str()) ));

	canvas_v[j]->cd();
	histo->Draw();

	func_vector[j]->Draw("SAME"); 


	//std::cout << My_Function::LagrangePol(My_Function::bin_centers, My_Function::pt_params[j], 580) << std::endl;
	j++;
	
}


}





std::vector<TFitResultPtr> My_Function::pol_fit(std::vector<TH1*> ratio_vector, TString tag){
	
	// 
	//+ x * x * x * x * x * [5] + x * x * x * x * x * x * [6]
	//+  x * x * x * [3] + x * x * x * x * [4] + x * x * x * x * x * [5] + x * x * x * x * x * x * [6]
	//+ [2] * x * x +  x * x * x * [3] 
	
		
	std::vector<TFitResultPtr> res;


	return res;


}




#endif


