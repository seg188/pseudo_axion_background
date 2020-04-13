#include "GenerateAnalysis.cpp"
#include "../fit/GPR/GP1D_Fitter.hh"
#include "../fit/GPR/GP1D_Fitter.cpp"
#include "../fit/GPR/test_gpr.cpp"
#include "../DataSet.hh"
#include "../../include/UserIn.cpp"
#include "../fit/drawCombined.cpp"


///////////////////////////////////////////////////////////////////////////////////////////////////////
//CALLING GetBackGroundEstimate() returns a std::vector of "function" class
//To use the class, the operator () is overloaded for x vales (w-mass values)
//method sigma() returns 1-std deviation error


// STEPS FOR GETTING BACKGROUND ESTIMATE ----------------------------------------------------------------------------------------------------------------------
class fit_result_function{
public:
	int j =0;
	std::vector<double> _x_axis;
	std::vector<double> _y_axis;
	std::vector<double> _errors;
	bool void_plot = false;

	fit_result_function(TH1* plot){
		std::cout << "resizing!!" << std::endl;
		if (!plot) void_plot = true;
		else{
			_x_axis.resize(plot->GetNbinsX());
			_y_axis.resize(plot->GetNbinsX());
			_errors.resize(plot->GetNbinsX());

			std::cout << "stting points!!" << std::endl;
			for (int j = 0; j < plot->GetNbinsX(); j++){
			_x_axis[j] = plot->GetBinCenter(j+1);
			_y_axis[j] = plot->GetBinContent(j+1);
			_errors[j] = plot->GetBinError(j+1);
		}
		}
		
		
	} 

	double operator()(double x){
		//std::cout << "call number: " << std::endl;
		//std::cout << j++ << std::endl;
		//std::cout << "brr" << std::endl;
		if (void_plot) return 0.00;
		for (int k = 0; k < _x_axis.size(); k++){
			if (x > _x_axis[k] and x < _x_axis[k+1]){
				return _y_axis[k];
			}
		}

		return 0.0;
	}

	double sigma(double x){
		for (int k = 0; k < _x_axis.size(); k++){
			if (x > _x_axis[k] and x < _x_axis[k+1]){
				return _errors[k];
			}
		}

		return 10000.00;
		
	}

};


class product_function{
public:
	fit::GP1D_Fitter* _fit;
	fit_result_function* _ratio;

	product_function(fit::GP1D_Fitter* fit, fit_result_function* ratio){
		_fit = fit;
		_ratio = ratio;
	}

	double operator()(double x){
		return (*_ratio)(x) * _fit->mean(x);
	}

	double sigma(double x){
		return TMath::Sqrt(TMath::Power(_fit->error(x) * (*_ratio)(x), 2.0) + TMath::Power(_fit->mean(x) * _ratio->sigma(x), 2.0));
	}
};

class sum_function{
public:
	std::vector<product_function> _functions;
	std::vector<double> _weights;
	std::vector<double> _weight_sigma;

	sum_function(std::vector<product_function> functions, std::vector<double> weights, std::vector<double> weight_sigma){
		_functions = functions;
		_weights = weights;
		_weight_sigma = weight_sigma;
	}

	double operator()(double x){
		double value = 0;
		for (int n =0; n < _weights.size(); n++){
			value += _weights[n] * (_functions[n])(x);
		}
		return value;
	}

	double sigma(double x){
		double sigma2 = 0;
		for (int k = 0; k < _functions.size(); k++){
			sigma2 += TMath::Power(_functions[k](x) * _weight_sigma[k], 2.00) + TMath::Power(_functions[k].sigma(x) * _weights[k], 2.00);
		}
		return TMath::Sqrt(sigma2);
	}

		

};


auto GetBackGroundEstimate(const char* file_folder = "~/hex/"){

//IMPORTING OR READING DATA TO BE USED IN BACKGROUND ESTIMATE/////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////

//GETTING RELEVANT DATA FROM GJETS MC TO MAKE PREDICTION

	TFile* gj_file = TFile::Open(file_folder + TString("gjets2016_mc.root"));
	auto Phi_Mass_PT = (TH2*) gj_file->Get("fPhiMass_PT_signal");

	std::vector<TH1*> gj_ratios;

	for (int k = 0; k < 9; k++){
		std::ostringstream strs;   
		strs << k;
		gj_ratios.push_back((TH1*) gj_file->Get(TString("gjets2016_mc_r_MassPi0_PT_") + TString(strs.str()) + TString("_ratio") ));
		
	}

												


//GETING WJETS MC AND DATA TO FIT RATIO PLOTS FOR SIDEBAND CORRECTION

	TFile* wj_file = TFile::Open(file_folder + TString("wjets2016_mc.root"));
	TString pre_ = TString("wjets2016_mc_r_MassPi0_PT_");
	TString _suf = TString("_ratio");

	std::vector<TH1*> wj_ratios;
	for (int k = 0; k < 9; k++){
		std::ostringstream strs;
		strs << k;
		wj_ratios.push_back((TH1*)wj_file->Get(pre_ + TString(strs.str()) + _suf ));
	}


	TFile* sm_file = TFile::Open(file_folder + TString("wjets_data.root"));
	std::vector<TH1*> sm_ratios;
	
	for (int k = 0; k < 9; k++){
		std::ostringstream strs;
		strs << k;
		sm_ratios.push_back((TH1*)sm_file->Get(TString("wjets_data_r_MassPi0_PT_") + TString(strs.str()) + TString("_ratio") ));
	}
                                                       

	std::vector<std::vector<TH1*>> ratios = {gj_ratios, wj_ratios, sm_ratios};
	std::cout << "beep boop" << std::endl;
	auto combined_ratios_th1 = fit::combine::drawCombined(ratios);
	std::cout << "making result functions" << std::endl;
	std::vector<fit_result_function> combined_ratios;
	for (auto fit : combined_ratios_th1){
		std::cout << fit->GetName() << std::endl;
		combined_ratios.push_back(fit_result_function(fit));
	}

//GETTING SIDE BANDS//////////////////////////////////////

	TFile* data_file = TFile::Open(file_folder + TString("sp_data.root"));
	std::vector<TH1*> plots;
		
	for (int k = 0; k < 9; k++){
		std::ostringstream strs;
		strs << k;
		plots.push_back((TH1*)data_file->Get(TString("sp_data_MassPi0_PT_") + TString(strs.str()) + TString("_side") ));

		
	}
	
	auto fits = fit_vector(plots, false);
	std::cout <<"FINISHED FITS" << std::endl;

	int plt = 0;
	std::vector<TCanvas*> canv;
	for (auto fit : fits){
		std::ostringstream strs;
		strs << plt;
		canv.push_back(new TCanvas(TString("canv") + TString(strs.str())) );
		plots[plt]->Draw();
	
		fit->_sigma_1_graph->Draw("SAME f");
		fit->_mean_graph->Draw("SAME");
		canv[plt]->Print(TString("plots/") + TString("sideband_fit_") + TString(strs.str() + TString(".pdf")), ".pdf");
		plt++;
	}

	


	std::vector<sum_function> predictions_by_phi;
	for (int phi_val =0; phi_val < histogram::PHI_MASS_NUMBINS; phi_val++){
		std::vector<double> pt_weights;
		std::vector<double> pt_weights_sigma;

		for (int pt_b = 0; pt_b < user_input::PT_BOUNDS1.size()-1; pt_b++){
			pt_weights.push_back(Phi_Mass_PT->GetBinContent(phi_val+1, pt_b+1));
			pt_weights_sigma.push_back(Phi_Mass_PT->GetBinError(phi_val+1, pt_b+1));

			//pt_weights.push_back(1.000);

		}

		std::vector<product_function> corrected_sideband;

		for (int pt_b = 0; pt_b < user_input::PT_BOUNDS1.size()-1; pt_b++){
			corrected_sideband.push_back(product_function(fits[pt_b], &(combined_ratios[pt_b])) );
		}

		predictions_by_phi.push_back(sum_function(corrected_sideband, pt_weights, pt_weights_sigma)); 

	}

	std::vector<TCanvas*> canvas_vector;
	std::vector<TGraph*> graphs;
	std::vector<TGraph*> err_graphs;
	int NPLOTPOINTS = 300;
	double x_min = 0.00;
	double x_max = 5.00;

	int j = 0;
	for (auto func : predictions_by_phi){
		std::ostringstream strs;
		strs << j;
		canvas_vector.push_back(new TCanvas("predictions_phi_" + TString(strs.str())));
		graphs.push_back(new TGraph(NPLOTPOINTS));
		err_graphs.push_back(new TGraph(2*NPLOTPOINTS));

		for (int xb = 0; xb < NPLOTPOINTS; xb++){
			double x = static_cast<double>(xb) / static_cast<double>(NPLOTPOINTS) *(x_max - x_min) + x_min;
			double x2 = static_cast<double>(NPLOTPOINTS - xb - 1) / static_cast<double>(NPLOTPOINTS) *(x_max - x_min) + x_min;
			graphs[j]->SetPoint(xb, x, func(x));

			err_graphs[j]->SetPoint(xb, x, func(x) + func.sigma(x));
			err_graphs[j]->SetPoint(NPLOTPOINTS + xb, x2, func(x2) - func.sigma(x2));
		}
		canvas_vector[j]->cd();
		graphs[j]->SetTitle(TString("Prediction_Phi_Region_") + TString(strs.str()));
		graphs[j]->SetLineWidth(1);
		graphs[j]->Draw();
		err_graphs[j]->SetFillStyle(3002);
    	err_graphs[j]->SetFillColor(2);
    	err_graphs[j]->Draw("SAME f");
		canvas_vector[j]->Print(TString("plots/") + TString(graphs[j]->GetTitle()) + TString(".pdf"), ".pdf");
		

		j++;
	}

	return predictions_by_phi;




	//scale by ratios
	//badabing badaboom


	//CREATING 

   	








}