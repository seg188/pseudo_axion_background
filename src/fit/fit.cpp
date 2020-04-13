#include "helper.cpp"
#include "../analysis/GenerateAnalysis.cpp"

namespace fit{
using namespace helper;

TF1* pol_fit_function1 = new TF1("pol_fit_function1", " [0] + [1] * x " , 1.0, 5.);
TF1* pol_fit_function2 = new TF1("pol_fit_function2", " [0] + [1] * x + [2] * x * x  " , 1.0, 5.);
TF1* pol_fit_function3 = new TF1("pol_fit_function3", " [0] + [1] * x + [2] * x * x +  x * x * x * [3] " , 1.0, 5.);
TF1* pol_fit_function4 = new TF1("pol_fit_function4", " [0] + [1] * x + [2] * x * x +  x * x * x * [3] + x * x * x * x * [4] " , 1.0, 5.);

void expand_coverage(TH1* plot, std::vector<TH1*> check_plots){
	double largest_expand_factor = 1.00;

	for (int x_bin = 1; x_bin <= plot->GetNbinsX(); x_bin++){
		
		if (plot->GetBinContent(x_bin) == 0 or plot->GetBinError(x_bin) == 0) continue;

		for (TH1* check_plot : check_plots){
			//case for avg > plot
			if (check_plot->GetBinContent(x_bin) == 0 or check_plot->GetBinError(x_bin) == 0) continue;
			double diff = get_diff_ratio(plot, check_plot, x_bin);
			if ( diff > largest_expand_factor ) largest_expand_factor = diff;

		}

		//plot->SetBinError(x_bin, plot->GetBinError(x_bin) * largest_expand_factor);

	}

	for (int x_bin = 0; x_bin < plot->GetNbinsX(); x_bin++){
		if (plot->GetBinContent(x_bin) == 0 or plot->GetBinError(x_bin) == 0) continue;
		else (plot->SetBinError(x_bin, plot->GetBinError(x_bin) * largest_expand_factor));
	}

}

template <typename Data_Set>
std::vector<TH1*> GetRatioPlots(Data_Set data){


	return analysis::GenerateAnalysis(data);


}

template <typename Data_Set>
std::vector<TH1*> FitRatios(Data_Set data){
	
	TString tag = data._tag;
	bool isdata = data._is_data;
	std::vector<TH1*> ratio_vector_sm = GetRatioPlots(data);

	std::cout << "BEEP BOOP" << std::endl;
	//std::vector<TFitResultPtr> sm_fit_results = My_Function::pol_fit(ratio_vector_sm, tag);

	std::vector<TFitResultPtr> ratio_fit_results;

	pol_fit_function1->SetParameters(1, 0, -1, -1);
	pol_fit_function2->SetParameters(1, 0, -1, -1);
	pol_fit_function3->SetParameters(1, 0, -1, -1);
	//pol_fit_function2->SetParLimits(2, -10.0, 1.00);
	std::vector<TCanvas*> canvas_v;
	std::vector<TH1*> fitted_ratios_with_errors;
	int fit_num = 0;

	for (TH1* ratio_plot : ratio_vector_sm){
		//ratio_plot->Rebin(2);

		std::cout << "WORKING ON: " << tag << " " << ratio_plot->GetName() << std::endl;
		int which_fit = F_Test(ratio_plot);
		std::cout << "out of F TEST***************************" << std::endl;
		double FIT_MIN = getFitMin(ratio_plot);
		double FIT_MAX = getFitMax(ratio_plot);
		canvas_v.push_back(new TCanvas(ratio_plot->GetName() + TString("canvas") + tag , tag));
		int num_data_points = 0;
		ratio_plot->Draw();
		for (int bin = 1; bin <= ratio_plot->GetNbinsX(); bin++){
			if (ratio_plot->GetBinContent(bin) > 0.001 and ratio_plot->GetBinCenter(bin) > FIT_MIN and ratio_plot->GetBinCenter(bin) < FIT_MAX) num_data_points++;
		}
		
		std::cout << "min: " << FIT_MIN << std::endl;
		std::cout << "max: " << FIT_MAX << std::endl;
		std::cout << "num points: " << num_data_points << std::endl;

		TFitResultPtr fit_result;
		//add F test right here
		if (which_fit == 2) fit_result = ratio_plot->Fit("pol_fit_function4", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == 1) fit_result = ratio_plot->Fit("pol_fit_function3", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == 0) fit_result = ratio_plot->Fit("pol_fit_function2", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == -1) fit_result = ratio_plot->Fit("pol_fit_function1", "S", "S", FIT_MIN, FIT_MAX);
		else {
			std::cout << "WHICH FIT: " << which_fit << std::endl;
			std::cout << "NOT DOING FIT" << std::endl;
			std::cout << "THIS IS NULLLLLL*****************************************************************" << std::endl;
			fitted_ratios_with_errors.push_back( new TH1D(ratio_plot->GetName() + TString("_fitted"), ratio_plot->GetName() + TString("_fitted"), 120, 0, 5.0) );
			fit_num++;
			continue;
		}

		fitted_ratios_with_errors.push_back( new TH1D(ratio_plot->GetName() + TString("_fitted"), ratio_plot->GetName() + TString("_fitted"), 500, 0, 5.0) );
		
		std::cout << "getting cov matrix" << std::endl;
		TMatrixD covv_matrix_v(fit_result->GetCovarianceMatrix());
		
		std::cout << "about to fill plots " << std::endl;
		std::vector<double> y_params = fit_result->Parameters();
		int x_bins =  fitted_ratios_with_errors[fit_num]->GetNbinsX();

		std::cout << "filling plots" << std::endl;

		for (int i = 1; i <= x_bins; i++){    
			double x = fitted_ratios_with_errors[fit_num]->GetBinCenter(i);
			double total_error_x = get_error(covv_matrix_v, x);


			if (x > FIT_MIN and x < FIT_MAX){

				fitted_ratios_with_errors[fit_num]->SetBinContent(i, GetFitValue(x, y_params));
				std::cout << "SETTING val: " <<  GetFitValue(x, y_params) << std::endl;
				fitted_ratios_with_errors[fit_num]->SetBinError(i, total_error_x);
			}

		}

		fit_num++;

	} //THIS IS THE IMPORTANT ONE

	std::cout << "returning plots" << std::endl;
				

	return fitted_ratios_with_errors;


}


std::vector<TH1*> FitRatios(std::vector<TH1*> plots, TString tag = TString("comb")){
	TF1* pol_fit_function1 = new TF1("pol_fit_function1", " [0] + [1] * x " , 1.0, 5.);
	TF1* pol_fit_function2 = new TF1("pol_fit_function2", " [0] + [1] * x + [2] * x * x  " , 1.0, 5.);
	TF1* pol_fit_function3 = new TF1("pol_fit_function3", " [0] + [1] * x + [2] * x * x +  x * x * x * [3] " , 1.0, 5.);
	TF1* pol_fit_function4 = new TF1("pol_fit_function4", " [0] + [1] * x + [2] * x * x +  x * x * x * [3] + x * x * x * x * [4] " , 1.0, 5.);

	
	std::vector<TH1*> ratio_vector_sm = plots;

	std::cout << "BEEP BOOP" << std::endl;
	//std::vector<TFitResultPtr> sm_fit_results = My_Function::pol_fit(ratio_vector_sm, tag);

	std::vector<TFitResultPtr> ratio_fit_results;
	std::cout << pol_fit_function1 << std::endl;
	pol_fit_function1->SetParameters(1, 1);
	pol_fit_function2->SetParameters(1, 0, -1);
	pol_fit_function3->SetParameters(1, 0, -1, -1);
	//pol_fit_function2->SetParLimits(2, -10.0, 1.00);
	std::vector<TCanvas*> canvas_v;
	std::vector<TH1*> fitted_ratios_with_errors;
	int fit_num = 0;

	for (TH1* ratio_plot : ratio_vector_sm){
		//ratio_plot->Rebin(2);

		std::cout << "WORKING ON: " << tag << " " << ratio_plot->GetName() << std::endl;
		int which_fit = F_Test(ratio_plot);
		std::cout << "out of F TEST***************************" << std::endl;
		double FIT_MIN = getFitMin(ratio_plot);
		double FIT_MAX = getFitMax(ratio_plot);
		canvas_v.push_back(new TCanvas(ratio_plot->GetName() + TString("canvas") + tag , tag));
		int num_data_points = 0;
		ratio_plot->Draw();
		for (int bin = 1; bin <= ratio_plot->GetNbinsX(); bin++){
			if (ratio_plot->GetBinContent(bin) > 0.001 and ratio_plot->GetBinCenter(bin) > FIT_MIN and ratio_plot->GetBinCenter(bin) < FIT_MAX) num_data_points++;
		}
		
		std::cout << "min: " << FIT_MIN << std::endl;
		std::cout << "max: " << FIT_MAX << std::endl;
		std::cout << "num points: " << num_data_points << std::endl;

		TFitResultPtr fit_result;
		//add F test right here
		if (which_fit == 2) fit_result = ratio_plot->Fit("pol_fit_function4", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == 1) fit_result = ratio_plot->Fit("pol_fit_function3", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == 0) fit_result = ratio_plot->Fit("pol_fit_function2", "S", "S", FIT_MIN, FIT_MAX);
		else if (which_fit == -1) fit_result = ratio_plot->Fit("pol_fit_function1", "S", "S", FIT_MIN, FIT_MAX);
		else {
			std::cout << "WHICH FIT: " << which_fit << std::endl;
			std::cout << "NOT DOING FIT" << std::endl;
			std::cout << "THIS IS NULLLLLL*****************************************************************" << std::endl;
			fitted_ratios_with_errors.push_back( new TH1D(ratio_plot->GetName() + TString("_fitted"), ratio_plot->GetName() + TString("_fitted"), 120, 0, 5.0) );
			fit_num++;
			continue;
		}

		fitted_ratios_with_errors.push_back( new TH1D(ratio_plot->GetName() + TString("_fitted"), ratio_plot->GetName() + TString("_fitted"), 500, 0, 5.0) );
		
		std::cout << "getting cov matrix" << std::endl;
		TMatrixD covv_matrix_v(fit_result->GetCovarianceMatrix());
		
		std::cout << "about to fill plots " << std::endl;
		std::vector<double> y_params = fit_result->Parameters();
		int x_bins =  fitted_ratios_with_errors[fit_num]->GetNbinsX();

		std::cout << "filling plots" << std::endl;

		for (int i = 1; i <= x_bins; i++){    
			double x = fitted_ratios_with_errors[fit_num]->GetBinCenter(i);
			double total_error_x = get_error(covv_matrix_v, x);


			if (x > FIT_MIN and x < FIT_MAX){

				fitted_ratios_with_errors[fit_num]->SetBinContent(i, GetFitValue(x, y_params));
				std::cout << "SETTING val: " <<  GetFitValue(x, y_params) << std::endl;
				fitted_ratios_with_errors[fit_num]->SetBinError(i, total_error_x);
			}

		}

		fit_num++;

	} //THIS IS THE IMPORTANT ONE



	std::cout << "returning plots" << std::endl;
				

	return fitted_ratios_with_errors;




}


template <typename set1, typename set2, typename set3>
std::vector<TH1*> drawRatios(set1 data1, set2 data2, set3 data3){

	std::vector<std::vector<TH1*>> ratio_histos;
	ratio_histos.push_back(GetRatioPlots(data1));
	ratio_histos.push_back(GetRatioPlots(data2));
	ratio_histos.push_back(GetRatioPlots(data3));

	std::vector<TH1*> avg_plots;
	avg_plots.resize(ratio_histos[0].size());


	for (int k = 0; k < avg_plots.size(); k++){
		std::vector<TH1*> combine = {ratio_histos[0][k], ratio_histos[1][k], ratio_histos[2][k]};
		avg_plots[k] = get_harmonic_average(combine, ratio_histos[0][k]->GetName() + TString("avg"));
	}

		std::vector<TH1*> avg_plots_fit = FitRatios(avg_plots, "_avg");

        std::vector<TH1*> ratio_plots = FitRatios(ratio_histos[0], data1._tag);

        //TCanvas* ratio_canvas = new TCanvas("ratios", "ratios");
        //ratio_canvas->cd();
        //ratio_plots[0]->Draw();
        //
        TLegend* all_legend = new TLegend(0.71, 0.77, .98, .98);
        all_legend->SetHeader("Sample", "C");
        std::vector<TCanvas*> canvas_v;
        int num_plots =  ratio_plots.size();
        for (int j = 0; j < num_plots; j++){

                canvas_v.push_back(new TCanvas(getName(j, "canvas_"), ratio_plots[j]->GetName()));
                canvas_v[j]->cd();
                ratio_plots[j]->SetLineColor(kRed);
                ratio_plots[j]->SetMaximum(2.);
                ratio_plots[j]->SetLineWidth(6);
                ratio_plots[j]->SetMinimum(0.0);
                ratio_plots[j]->SetMaximum(3.5);
                //ratio_plots[j]->SetMarkerStyle(kFullSquare);
                //ratio_plots[j]->SetMarkerColor(kBlack);
                ratio_plots[j]->SetFillColor(kRed);
         	//	ratio_plots[j]->Rebin(2);
                if (j == 0) all_legend->AddEntry(ratio_plots[j], "wjets mc");
                ratio_plots[j]->Draw();
        }

      
        std::vector<TH1*> ratio_plots1 = FitRatios(ratio_histos[1], data2._tag);


        for (int j = 0; j < num_plots; j++){
                canvas_v[j]->cd();
                ratio_plots1[j]->SetLineColor(kBlue);
                ratio_plots1[j]->SetMaximum(2.);
                ratio_plots1[j]->SetLineWidth(6);
                //ratio_plots[j]->SetMarkerStyle(kFullDotMedium);
                //ratio_plots[j]->SetMarkerColor(kBlack);
                ratio_plots1[j]->SetFillColor(kBlue);
            //    ratio_plots[j]->Rebin(2);
                if (j == 0) all_legend->AddEntry(ratio_plots1[j], "wjets data");
                ratio_plots1[j]->Draw("SAME");
        }

       
         std::vector<TH1*> ratio_plots2 = FitRatios(ratio_histos[2], data3._tag);


        for (int j = 0; j < num_plots; j++){
                canvas_v[j]->cd();
                ratio_plots2[j]->SetLineColor(kGreen);
                ratio_plots2[j]->SetMaximum(2.);
                ratio_plots2[j]->SetLineWidth(6);
                ratio_plots2[j]->SetFillColor(kGreen);
                if (j == 0) all_legend->AddEntry(ratio_plots2[j], "gjets mc");
                ratio_plots2[j]->Draw("SAME");
                all_legend->Draw("SAME");
                avg_plots_fit[j]->SetLineColor(kViolet);
                avg_plots_fit[j]->Draw("SAME");
         }


         return avg_plots_fit;
}	



TH1* getRatio(TH1* top, TH1* bot, TString name){

	int num_x_bins = top->GetNbinsX();
	double min_x = 0.0;
	double max_x = 5.0;

	TH1* ratio_plot = new TH1D(name, name, num_x_bins, min_x, max_x);

	for (int k = 1; k <= num_x_bins; k++){
		if (bot->GetBinContent(k) < 0.001) continue;
		ratio_plot->SetBinContent(k, top->GetBinContent(k) / bot->GetBinContent(k));
		ratio_plot->SetBinError(k, getError(top->GetBinContent(k), bot->GetBinContent(k), top->GetBinError(k), bot->GetBinError(k)));
	}

	return ratio_plot;

}

std::vector<TH1*> getVectorRatio(std::vector<TH1*> top_v, std::vector<TH1*> bot_v, TString id){
	std::vector<TH1*> ratio_vector;

	for (int k = 0; k < top_v.size(); k++) {
		ratio_vector.push_back(getRatio(top_v[k], bot_v[k], getName1(TString("pt_") + TString(id), k)));
	}

	return ratio_vector;
}



}; // NAMESPACE fit