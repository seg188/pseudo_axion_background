#include "fit.cpp"
namespace fit{
namespace combine{


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



int expand_errors(TH1* plot, double val = 1.050){
	for (int k =1; k <= plot->GetNbinsX(); k++){
		if (plot->GetBinError(k) < 0.001 or plot->GetBinContent(k) < 0.001) continue;
		//std::cout <<  "setting " << plot->GetBinError(k) <<" to " << val * plot->GetBinError(k) << std::endl;
		plot->SetBinError(k, val * plot->GetBinError(k));
		if (plot->GetBinError(k) > 999999999.9) return -1;
	}

	return 0;
}
double chi_2(TH1* plot1, TH1* plot2){
   double val = 0.000;
   double dof = -2.0;
   for (int x = 1; x <= plot2->GetNbinsX(); x++){
   	  
   	  int x2 = x * plot1->GetNbinsX() / plot2->GetNbinsX();
      double diff = TMath::Abs(plot2->GetBinContent(x) - plot1->GetBinContent(x2));
      double err = TMath::Power(plot1->GetBinError(x2), 2.00) + TMath::Power(plot2->GetBinError(x), 2.00);
      double v1 = plot1->GetBinContent(x2);
      double v2 = plot2->GetBinContent(x);
      if (err < 0.001 or v1 < 0.0001 or v2 < 0.0001) continue;
    //  if (plot1->GetBinContent(x) < 0.000001 or plot2->GetBinContent(x) < 0.000001) continue;
   	  	  val += diff / TMath::Sqrt(err);
    	  dof += 1.0;
   }

   if (dof < 1) return 0;
   return val / dof;
}

void expand_cov(TH1* exp, TH1* ref){
	TH1* copy = (TH1*) exp->Clone();
	double chi_2_last = -1.0;
	while (chi_2(exp, ref) > 1.00){
		//if (chi_2_last == chi_2(exp, ref)) break;
		//std::cout << chi_2(exp, ref) << std::endl;
		expand_errors(exp);
		chi_2_last = chi_2(exp, ref);


	}

}



std::vector<TH1*> drawCombined(std::vector<std::vector<TH1*>> ratio_histos, int draw_opt = 0, int opt =1){
	
	std::vector<TH1*> avg_plots;
	avg_plots.resize(ratio_histos[0].size());


	for (int k = 0; k < avg_plots.size(); k++){
		std::vector<TH1*> combine = {ratio_histos[0][k], ratio_histos[1][k], ratio_histos[2][k]};
		avg_plots[k] = get_harmonic_average(combine, ratio_histos[0][k]->GetName() + TString("avg"));
	}

	std::vector<TH1*> avg_plots_fit = FitRatios(avg_plots);

	std::vector<TCanvas*> canvas_v;
	
	std::vector<TH1*> avg_plots_fit_copy;


	for (int k = 0; k < avg_plots_fit.size(); k++){
		canvas_v.push_back(new TCanvas(TString(avg_plots_fit[k]->GetName()) + TString("canv"), "canv"));
		avg_plots_fit_copy.push_back( (TH1 *) avg_plots_fit[k]->Clone());
		expand_cov(avg_plots_fit[k], ratio_histos[0][k]);
		expand_cov(avg_plots_fit[k], ratio_histos[1][k]);
		expand_cov(avg_plots_fit[k], ratio_histos[2][k]);

		canvas_v[k]->Divide(2, 2);
	
		avg_plots[k]->SetLineColor(kBlack);
		avg_plots[k]->SetMarkerColor(kBlack);
		avg_plots[k]->SetMarkerStyle(kFullSquare);
		avg_plots_fit[k]->SetLineColor(kCyan);
		avg_plots_fit[k]->SetLineWidth(6);

		ratio_histos[0][k]->SetLineColor(kRed);
		ratio_histos[0][k]->SetMarkerColor(kRed);
		ratio_histos[0][k]->SetMarkerStyle(kFullSquare);
		ratio_histos[1][k]->SetLineColor(kBlue);
		ratio_histos[1][k]->SetMarkerColor(kBlue);
		ratio_histos[1][k]->SetMarkerStyle(kFullSquare);
		ratio_histos[2][k]->SetLineColor(kGreen);
		ratio_histos[2][k]->SetMarkerColor(kGreen);
		ratio_histos[2][k]->SetMarkerStyle(kFullSquare);

		avg_plots_fit[k]->SetMaximum(3.5);
		avg_plots_fit[k]->SetMinimum(0.0);
		canvas_v[k]->cd(1);
		avg_plots_fit_copy[k]->SetLineColor(kViolet);
		avg_plots_fit_copy[k]->SetLineWidth(6);


		avg_plots_fit[k]->Draw();
		std::ostringstream strs;
		strs << chi_2(avg_plots_fit[k], avg_plots[k]);
		avg_plots_fit_copy[k]->Draw("SAME");
		avg_plots[k]->Draw("SAME");
		getLegend(TString("chi2=") + TString(strs.str()))->Draw("SAME");

		canvas_v[k]->cd(2);
		std::ostringstream strs1;
		strs1 << chi_2(avg_plots_fit[k], ratio_histos[0][k]);
		avg_plots_fit[k]->Draw();
		avg_plots_fit_copy[k]->Draw("SAME");
		ratio_histos[0][k]->Draw("SAME");
		getLegend(TString("chi2=") + TString(strs1.str()))->Draw("SAME");

		std::ostringstream strs2;
		strs2 << chi_2(avg_plots_fit[k], ratio_histos[1][k]);
		
		canvas_v[k]->cd(3);
		avg_plots_fit[k]->Draw();
		avg_plots_fit_copy[k]->Draw("SAME");
		ratio_histos[1][k]->Draw("SAME");
		getLegend(TString("chi2=") + TString(strs2.str()))->Draw("SAME");
		std::ostringstream strs3;
		strs3 << chi_2(avg_plots_fit[k], ratio_histos[2][k]);
		
		canvas_v[k]->cd(4);
		avg_plots_fit[k]->Draw();
		avg_plots_fit_copy[k]->Draw("SAME");
		ratio_histos[2][k]->Draw("SAME");
	
		getLegend(TString("chi2=") + TString(strs3.str()))->Draw("SAME");


	
	}

	for (TCanvas* canvas : canvas_v){
      
          canvas->Print(TString("/home/stephen/hex/wjets/fit/plots_mar11/") + canvas->GetName() + TString(".pdf"), "pdf");
    } 
    
   

    return avg_plots_fit;
}
}; //namepsce combine
}; //namespace fit