#include <iostream>
#include <math.h>
#include <dirent.h>

#ifndef CISO_STATS
#define CISO_STATS

namespace analysis::helper{

/*
    For each observed number in the table subtract the corresponding expected number (O — E).
    Square the difference [ (O —E)2 ].
    Divide the squares obtained for each cell in the table by the expected number for that cell [ (O - E)2 / E ].
    Sum all the values for (O - E)2 / E. This is the chi square statistic. 
*/

//chi sqaured test for expectation value always equals 1, with error bars, then normalized by DOF

double getError(double x, double y, double ex, double ey) //for quanity x/y
{
   double error = ((ey * ey)/(y * y) + (ex * ex)/(x * x));
   error = sqrt(error);
   error = (x / y)*error;
   return error;
}

std::vector<std::string> getFiles(const char* dir_name){
 	std::vector<std::string> f_Files;
  	std::string path_to_dir = dir_name;
  	std::vector<std::string> error1;
  	error1.push_back("error... could not open directory");

  auto dir = opendir(path_to_dir.c_str());
  if (dir == NULL)
  {
    std::cout << "Could not open directory" << std::endl;
    return error1;
  }
  auto entity = readdir(dir);
  while (entity != NULL)
  {
    if(entity->d_type == DT_REG) 
    {
      TString file_name = TString(entity->d_name);
      if (file_name.EndsWith(".root")) f_Files.push_back(std::string(entity->d_name));
    }

    entity = readdir(dir);
  }

  closedir(dir);

  std::vector<std::string> f_Files2;
  f_Files2.push_back("ttb_output.root");
  f_Files2.push_back("dy_output.root");
  f_Files2.push_back("wj_output.root");
  f_Files2.push_back("qcd_output.root");
  f_Files2.push_back("sm_output.root");

  return f_Files2;
}

std::vector<EColor> color_v{kRed, kBlue, kOrange , kGreen,  kViolet};

EColor getColor(int index, int max_vals = 6){

	return color_v[index % max_vals];

}

TString getFileName(std::string name){
	if (name == "dy_output.root") return TString("DY mc");
	else if (name == "qcd_output.root") return TString("QCD mc");
	else if (name == "wj_output.root") return TString("wjets mc");
	else if (name == "ttb_output.root") return TString("ttbar mc");
	else return TString("unknown");
}

TString getFileNameFull(std::string name){
	if (name == "/home/stephen/hex/wjets/output_files/dy_output.root") return TString("DY_");
	else if (name == "/home/stephen/hex/wjets/output_files/qcd_output.root") return TString("QCD_mc_");
	else if (name == "/home/stephen/hex/wjets/output_files/wj_output.root") return TString("wjets_mc_");
	else if (name == "/home/stephen/hex/wjets/output_files/sm_output.root") return TString("sm_data_");
	else if (name == "/home/stephen/hex/wjets/output_files/ttb_output.root") return TString("ttbar_mc_");
	else return TString("unknown");
}

void Scale(TH1* histo, double scale_factor){

	//histo->SetBinErrorOption(TH1::kPoisson2);
	std::vector<double> errup, errdown;
	for (int j = 1; j <= histo->GetNbinsX(); j++){
		if (histo->GetBinContent(j) * scale_factor < 0.033){
			errup.push_back(histo->GetBinErrorUp(j));
			errdown.push_back(histo->GetBinErrorLow(j));
			std::cout << errup[j] * scale_factor << " " << errdown[j] * scale_factor << std::endl;
		} else {
			errup.push_back(-100);
			errdown.push_back(-100);
		}
		
	}

	histo->SetBinErrorOption(TH1::kNormal);

	for (int j = 1; j <= histo->GetNbinsX(); j++){
		if (errup[j] > -1) histo->SetBinError(j, (errup[j-1] + errdown[j - 1]) / 2 * scale_factor );
		if (errup[j] < -1) histo->SetBinError(j, histo->GetBinError(j) * scale_factor);
		histo->SetBinContent(j, histo->GetBinContent(j) * scale_factor);
	}	


}

double ChiSqaured(std::vector<double> observed_data, std::vector<double> error) 
{
	double chi_sq(0);
	double expectation(1); 
	double error_val(0);
	double val(0);
	int val_nan(0);


	for (int iii = 0; iii < observed_data.size(); iii ++)
	{
		error_val = error[iii];

		if (error_val > 0 and !isnan(error_val) and !isnan(observed_data[iii])){
			val = (observed_data[iii] - expectation);
			val = val / error_val;
			val = val * val;
			chi_sq = chi_sq + val;
		
		}
		else{

			val_nan++;
		}
		
		
		
	}

	return chi_sq/(observed_data.size() - val_nan - 2);
} 

TString getSciNotation(double num){ //get scientific notation for order of magntitude, no decimal places, 3 sig figs

	bool go(true);
	int magntitude(0);
	while (go){

		if (num / TMath::Power(10, magntitude) < 10){

			go = false;
		} else {

			magntitude++;
		}
	}

	std::ostringstream strs, strsmag;
	strs << num;

	TString prelim = TString(strs.str());

	strs.str("");
	strs.clear();

	strs << prelim(0);
	strs << ".";
	strs << prelim(1, 2);
	strsmag << magntitude;

	return TString(strs.str()) + TString("e") + TString(strsmag.str());

}

TString getName(TString cat, double num1, double num2){
	std::ostringstream strs1, strs2;
	strs1 << num1;
	strs2 << num2;

	return cat + TString(strs1.str()) + TString("_to_") + TString(strs2.str());
}


std::vector<TFitResultPtr> pol_fit(std::vector<TH1*> ratio_vector, TString tag);
double getFitMin(TH1*);
double getFitMax(TH1*);


std::vector<int> PT_BOUNDS2{20, 80, 180, 300, 500, 1000};
std::vector<int> PT_BOUNDS1{20, 40, 60, 80, 120, 160, 240, 300, 460, 700};
std::vector<int> PT_BOUNDS0{20, 40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, 130, 140, 160, 170, 180, 200, 240, 320, 460};

TH1* getRatios(TH1* tight, TH1* loose, bool draw, TString tag){
	int nonzero = 0;
	TH1* ratio = new TH1D(loose->GetName() + TString("ratio") + tag, loose->GetName() + TString("ratio") + tag, 30, 0.0, 5.0 );
	for (int k = 1; k <= tight->GetNbinsX(); k++){
		double t = tight->GetBinContent(k);
		double l = loose->GetBinContent(k);
		double el = loose->GetBinError(k);
		double et = tight->GetBinError(k);

		if ( getError(t, l, et, el) < 5.0 and l > 0.001) {//  getError(t, l, et, el) < 3.0 and l
			ratio->SetBinContent(k, t / l );
			ratio->SetBinError(k, getError(t, l, et, el)); //for quanity x/y
			nonzero++;
		}

	}

	ratio->SetLineColor(kViolet);
	//ratio->SetMarkerColorAlpha(kViolet, 0);

	return ratio;
	
}

void Scale(TH1* plot){
	std::cout << "\n\n\n******************" << std::endl;
	std::cout << plot->GetName() << std::endl;
		//Scale(MassPi0_PT_proj[pt_index]);

		for (int j = 1; j < plot->GetNbinsX(); j++){
			std::cout << j << ": " << "err: " << plot->GetBinError(j)  <<  std::endl;
		}
	
	
		plot->SetBinErrorOption(TH1::kNormal);

		for (int j = 1; j < plot->GetNbinsX(); j++){
			std::cout << j << ": " << "err: " << plot->GetBinError(j)  <<  std::endl;
			//std::cout << j << ": " << "err: " << plot->GetBinErrorUp(j)  <<  std::endl;
			//std::cout << j << ": " << "err: " << plot->GetBinErrorLow(j)  <<  std::endl;
		}
		std::cout << "******************\n\n\n" << std::endl;


	std::vector<double> error_up;
	std::vector<double> error_down;

	for (int j = 1; j <= plot->GetNbinsX(); j++){
		error_down.push_back(plot->GetBinErrorLow(j));
		error_up.push_back(plot->GetBinErrorUp(j));
	}
	plot->SetBinErrorOption(TH1::kNormal);
	double integral = plot->Integral();
	plot->Scale(1 / plot->Integral());
	for (int j = 1; j <= plot->GetNbinsX(); j++){
		if (plot->GetBinContent(j) != 0) {
			plot->SetBinError(j, (error_up[j - 1] + error_down[j -1 ]) / (2 * integral) );
		
			std::cout << j << ": " << "err: " << (error_down[j-1] + error_up[j - 1])/ (2 * integral) <<  std::endl;
			
		}

			
	}

}
TString getPTRange(int pt_index, std::vector<int> pt_bounds){
   std::ostringstream strs1, strs2;

   strs1 << pt_bounds[pt_index];
   strs2 << pt_bounds[pt_index + 1];

   return TString(strs1.str()) + TString("_to_") + TString(strs2.str());
}


void FgScale(TH1* plot){
	std::vector<double> vals;
	std::vector<double> errs;
	double integral = plot->Integral();

	for (int x = 0; x < plot->GetNbinsX(); x++){
		vals.push_back(plot->GetBinContent(x));
		errs.push_back(plot->GetBinError(x));
	}
//	plot->SetBinErrorOption(TH1::kNormal);

	for (int x = 0; x < plot->GetNbinsX(); x++){
		//std::cout << "bin val: " << vals[x] << std::endl;
		//std::cout << "bin err: " << errs[x] << std::endl;
		plot->SetBinContent(x + 1, vals[x] / integral);
		plot->SetBinError(x + 1, errs[x] / integral);
	}
}

}; //NAMESPACE helper


#endif