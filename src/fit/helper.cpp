#include "My_Function.h"
#include "My_Function.cpp"

namespace fit::helper{

TString getName(int number, TString tag = "pt_"){
        std::ostringstream strs_num;
        strs_num << number;

        return tag + TString(strs_num.str());
}

TLegend* getLegend(TString title){
	TLegend* legend = new TLegend(0.71, 0.77, .98, .98);
	legend->SetHeader(title);
	return legend;
}



TString getName1(TString cat, int num1){
	std::ostringstream strs1;
	strs1 << num1;

	return cat + TString(strs1.str());
}

double PolChangeSquareUp(std::vector<double> par_changes, std::vector<double> real_params, double x){
	double return_val = 0;
	for (int k = 0; k < real_params.size(); k++){
		return_val += TMath::Power(x, k) * (par_changes[k]);
		//std::cout << par_changes[k] << std::endl;
		//std::cout << par_changes[k] + real_params[k] << std::endl;
	}

	return return_val * return_val;
}

double PolChangeSquareLow(std::vector<double> par_changes, std::vector<double> real_params, double x){
	double return_val = 0;
	for (int k = 0; k< real_params.size(); k++){
		return_val += TMath::Power(x, k) * par_changes[k];
	}

	return return_val * return_val;
}

double GetFitValue(double x, std::vector<double> real_params){
	double return_val = 0;
	for (int k = 0; k< real_params.size(); k++){
		return_val += TMath::Power(x, k) * (real_params[k]);
	}
	return return_val;
}

void DrawRatioFitEnvelopes(std::vector<std::vector<std::vector<double>>> par_low_mean_up){


}

double polynomial(std::vector<double> params, double x){
	double val = 0;
	for (int k = 0; k < params.size(); k++){
		val += params[k] * TMath::Power(x, k);
	}
	return val;
}

double deriv_polynomial(std::vector<double> params, double x){
	double val = 0;
	double power = 1.00;
	for (int k = 1; k < params.size(); k++){
		val += power * params[k] * TMath::Power(x, k - 1);
		power += 1.00;
	}
	return val;
}

std::vector<double> get_Fk(TMatrixD eigen_vector_matrix, int k){
	std::vector<double> parameters;
	for (int row = 0; row < eigen_vector_matrix.GetNrows(); row++){
		parameters.push_back(eigen_vector_matrix[row][k]);
		std::cout << "FK: " << eigen_vector_matrix[row][k] << std::endl;
	}

	return parameters;
}

double dot_product_2(std::vector<double> one, std::vector<double> two){
	double val = 0;
	for (int k  = 0; k < one.size(); k++){
		val += one[k] * two[k];
	}

	return val;
}

double get_error(TMatrixD cov_matrix, double x){
	double val = 0;
	for (int col = 0; col < cov_matrix.GetNcols(); col++){
		for (int row = 0; row < cov_matrix.GetNrows(); row++){
			val += cov_matrix[row][col] * TMath::Power(x, row + col); 
		}
	}

	return TMath::Sqrt(val);
}

double get_x_max(TH1* plot){
 	int num_bins = plot->GetNbinsX();
 	int max_x = 5.0;
 	for (int bin = 10; bin <= num_bins; bin++){ 
 		if (plot->GetBinContent(bin) < 0.005 and plot->GetBinContent(bin + 1) < 0.005) {
 			max_x = plot->GetBinCenter(bin);
 			break;		
 		}
 	}

 	return max_x;
 }

double get_x_min(TH1* plot){
 	int num_bins = plot->GetNbinsX();
 	int min_x = 0.0;
 	for (int bin = 1; bin <= num_bins; bin++){ 
 		if (plot->GetBinContent(bin) > 0.005 and plot->GetBinContent(bin + 1) > 0.005){
 			min_x = plot->GetBinCenter(bin);		
 			break;
 		} 
 	}

 	return min_x;
 }


double getFitMin(TH1* plot){
	int nbins = plot->GetNbinsX();
	int begin_bin = 1;
	for (int j = 1; j <= nbins; j++){
		if (plot->GetBinContent(j) > 0.04 and plot->GetBinContent(j + 1) > 0.04){
			begin_bin = j;
			break;
		}
	}

	return plot->GetBinCenter(begin_bin);
}

double getFitMax(TH1* plot){

	int nbins = plot->GetNbinsX();
	int end_bin = nbins;
	for (int j = 7; j < nbins; j++){
		if (plot->GetBinContent(j) < 0.04 and plot->GetBinContent(j + 1) < 0.04 and plot->GetBinContent(j + 2) < 0.04){
			end_bin = j;
			break;
		}
	}

	return plot->GetBinCenter(end_bin);
	
}
	//  + x * x * x * x * x * x * x* [7]
int get_max_index(std::vector<double> vec){
	double max = -100000.0;
	int max_index = -1;
	for (int k = 0; k < vec.size(); k++){
		if (vec[k] > max){
			max = vec[k];
			max_index = k;
		} 
	}

	return max_index;
}
int F_Test(TH1* plot){
	std::cout << "IN F TEST***************************" << std::endl;
	if (plot == nullptr) return -2;
	std::vector<TFitResultPtr> fit_results;
	double FIT_MIN = getFitMin(plot);
	double FIT_MAX = getFitMax(plot);
	int n_data_points = 0;

	for (int bin = 1; bin <= plot->GetNbinsX(); bin++){
		if (plot->GetBinContent(bin) > 0.001 and plot->GetBinCenter(bin) > FIT_MIN and plot->GetBinCenter(bin) < FIT_MAX) n_data_points++;
	}

	std::cout << "N DATA POINTS: " << n_data_points << std::endl;
	if (n_data_points < 2) return -2;

	fit_results.push_back(plot->Fit("pol_fit_function1", "S", "Q", FIT_MIN, FIT_MAX));
	fit_results.push_back(plot->Fit("pol_fit_function2", "S", "Q", FIT_MIN, FIT_MAX));
	fit_results.push_back(plot->Fit("pol_fit_function3", "S", "Q", FIT_MIN, FIT_MAX));
	fit_results.push_back(plot->Fit("pol_fit_function4", "S", "Q", FIT_MIN, FIT_MAX));

	std::vector<double> chi_2_v, npar_v, ndof_v;
	int npar = 1;
	for (TFitResultPtr res : fit_results){
		chi_2_v.push_back(res->Chi2());
		npar_v.push_back(npar);
		ndof_v.push_back(res->Ndf());
		npar++;
	}

	std::vector<double> F_statistics;

	for (int k = 0; k <  fit_results.size() - 1; k++){
		double F_stat = chi_2_v[0] - chi_2_v[k + 1];
		F_stat = F_stat / (npar_v[k + 1] - 1);
		F_stat = F_stat / chi_2_v[k + 1];
		F_stat = F_stat * (n_data_points - npar_v[k + 1]);
		F_stat = TMath::FDistI(F_stat, 1, n_data_points - npar_v[k + 1] );
		F_statistics.push_back(F_stat);
	}

	bool to_continue = false;
	int good_index = -1;
	
	if (F_statistics[get_max_index(F_statistics)] > 0.90){
		good_index = get_max_index(F_statistics);
		to_continue = true;
	}
	

	if (!to_continue) {
		std::cout << "RETURNING INDEX: " << good_index << std::endl;
		return good_index;
	}

	std::vector<double> F_statistics2;

	for (int k = good_index; k <  fit_results.size() - 1; k++){
		double F_stat = chi_2_v[good_index] - chi_2_v[k + 1];
		F_stat = F_stat / (npar_v[k + 1] - npar_v[good_index]);
		F_stat = F_stat / chi_2_v[k + 1];
		F_stat = F_stat * (n_data_points - npar_v[k + 1]);
		F_stat = TMath::FDistI(F_stat, 1, n_data_points - npar_v[k + 1] );
		F_statistics2.push_back(F_stat);
	}

	if (F_statistics2[get_max_index(F_statistics2)] > 0.90){
		good_index = get_max_index(F_statistics2);
	}

	std::cout << "RETURNING INDEX: " << good_index << std::endl;

	return good_index;





}



TH1* get_harmonic_average(std::vector<TH1*> plots, TString name){ //FOR PLOTS THAT ALL HAVE SAME NUMBER OF BINS
	if (plots.size() == 0) return nullptr;
	TH1* avg_plot = new TH1D(name, name, plots[0]->GetNbinsX(), plots[0]->GetXaxis()->GetXmin(), plots[0]->GetXaxis()->GetXmax() );
	for (int x_bin = 1; x_bin <= plots[0]->GetNbinsX(); x_bin++){
		double bin_val = 0.0;
		double inv_bin_err = 0.0;

		for (TH1* plot : plots){
			double val = plot->GetBinContent(x_bin);
			double err = plot->GetBinError(x_bin);

			if (err == 0 or val == 0) continue;

			bin_val += val / (err * err);
			inv_bin_err += 1.0 / (err * err);

		}

		if (bin_val == 0 or inv_bin_err == 0) continue;

		avg_plot->SetBinContent(x_bin, bin_val / inv_bin_err );
		avg_plot->SetBinError(x_bin, TMath::Sqrt(1.0 / TMath::Sqrt(inv_bin_err)));

	}

	return avg_plot;
}

double get_min(std::vector<double> vals){
	double min = 99999999.0;
	for (double val : vals){
		if (val < min) min = val;
	}

	return min;
}

double get_diff_ratio(TH1* plot, TH1* check_plot, int bin){

	//return 1.0;
	double plot_val = plot->GetBinContent(bin);
	double chek_val = check_plot->GetBinContent(bin);
	double plot_err = plot->GetBinError(bin);
	double chek_err = check_plot->GetBinError(bin);


	double difference = TMath::Abs(plot_val - chek_val);

	if (difference < plot_err or difference < chek_err) return 0.999999;
	else{
		std::vector<double> vals;
		vals.push_back(TMath::Abs( plot_val - (chek_val + chek_err)  ));
		vals.push_back(TMath::Abs( plot_val - (chek_val - chek_err)  ));
		vals.push_back(TMath::Abs( (chek_val + chek_err) - plot_val));
		return get_min(vals) / plot_err;
	}
}

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

  return f_Files;
}

std::vector<EColor> color_v{kRed, kOrange, kGreen, kBlue, kViolet, kBlack};

EColor getColor(int index, int max_vals = 6){

	return color_v[index % max_vals];

}

TString getFileName(std::string name){
	if (name == "dy_output.root") return TString("DY + ttbar");
	else if (name == "qcd_output.root") return TString("QCD mc");
	else if (name == "wj_output.root") return TString("wjets mc");
	else return TString("unknown");
}

TString getFileNameFull(std::string name){
	if (name == "/home/stephen/hex/wjets/output_files/dy_output.root") return TString("DY_ttbar_");
	else if (name == "/home/stephen/hex/wjets/output_files/qcd_output.root") return TString("QCD_mc_");
	else if (name == "/home/stephen/hex/wjets/output_files/wj_output.root") return TString("wjets_mc_");
	else if (name == "/home/stephen/hex/wjets/output_files/sm_output.root") return TString("sm_data_");
	else return TString("unknown");
}

void Scale(TH1* histo, double scale_factor){

	histo->SetBinErrorOption(TH1::kPoisson2);
	std::vector<double> errup, errdown;
	for (int j = 1; j <= histo->GetNbinsX(); j++){
		if (histo->GetBinContent(j) * scale_factor < 0.033){
			errup.push_back(histo->GetBinErrorUp(j));
			errdown.push_back(histo->GetBinErrorLow(j));
			std::cout << "ERRORS: " << errup[j] * scale_factor << " " << errdown[j] * scale_factor << std::endl;
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





}// NAMESPACE helper