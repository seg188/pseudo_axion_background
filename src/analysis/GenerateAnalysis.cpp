#include "../DataSet.hh"
#include "../../include/UserIn.cpp"
#include "helper.cpp"



namespace analysis{

using namespace helper;

template <typename Data_Set>
std::vector<TH1*> GenerateAnalysis(Data_Set data, bool draw = false){

   std::vector<int> PT_BOUNDS = user_input::PT_BOUNDS1;
   TString tag = data._tag + TString("_r");

//GETTING DATA FROM OUTPUT FILE **********************************************************************

TH2* PT_Mass_Pi0_01 = data._PT_Mass_Pi0_01;
TH2* PT_Mass_Pi0_02 = data._PT_Mass_Pi0_02;
TH2* PT_Mass_Pi0_03 = data._PT_Mass_Pi0_03;
TH2* PT_Mass_Pi0_04 = data._PT_Mass_Pi0_04;
TH2* PT_Mass_Pi0_05 = data._PT_Mass_Pi0_05;
TH2* PT_Mass_Pi0_06 = data._PT_Mass_Pi0_06;
TH2* PT_Mass_Pi0_07 = data._PT_Mass_Pi0_07;
TH2* PT_Mass_Pi0_08 = data._PT_Mass_Pi0_08;
TH2* PT_Mass_Pi0_09 = data._PT_Mass_Pi0_09;
TH2* PT_Mass_Pi0_10 = data._PT_Mass_Pi0_10;


TH2* PT_Mass_Pi0_01_side = data._PT_Mass_Pi0_01_side;
TH2* PT_Mass_Pi0_02_side = data._PT_Mass_Pi0_02_side;
TH2* PT_Mass_Pi0_03_side = data._PT_Mass_Pi0_03_side;
TH2* PT_Mass_Pi0_04_side = data._PT_Mass_Pi0_04_side;
TH2* PT_Mass_Pi0_05_side = data._PT_Mass_Pi0_05_side;
TH2* PT_Mass_Pi0_06_side = data._PT_Mass_Pi0_06_side;
TH2* PT_Mass_Pi0_07_side = data._PT_Mass_Pi0_07_side;
TH2* PT_Mass_Pi0_08_side = data._PT_Mass_Pi0_08_side;
TH2* PT_Mass_Pi0_09_side = data._PT_Mass_Pi0_09_side;
TH2* PT_Mass_Pi0_10_side = data._PT_Mass_Pi0_10_side;

std::vector<TH2*> ptmass_v, ptmass_v_side;

ptmass_v.push_back(PT_Mass_Pi0_01);
ptmass_v.push_back(PT_Mass_Pi0_02);
ptmass_v.push_back(PT_Mass_Pi0_03);
ptmass_v.push_back(PT_Mass_Pi0_04);
ptmass_v.push_back(PT_Mass_Pi0_05);
ptmass_v.push_back(PT_Mass_Pi0_06);
ptmass_v.push_back(PT_Mass_Pi0_07);
ptmass_v.push_back(PT_Mass_Pi0_08);
ptmass_v.push_back(PT_Mass_Pi0_09);
ptmass_v.push_back(PT_Mass_Pi0_10);

ptmass_v_side.push_back(PT_Mass_Pi0_01_side);
ptmass_v_side.push_back(PT_Mass_Pi0_02_side);
ptmass_v_side.push_back(PT_Mass_Pi0_03_side);
ptmass_v_side.push_back(PT_Mass_Pi0_04_side);
ptmass_v_side.push_back(PT_Mass_Pi0_05_side);
ptmass_v_side.push_back(PT_Mass_Pi0_06_side);
ptmass_v_side.push_back(PT_Mass_Pi0_07_side);
ptmass_v_side.push_back(PT_Mass_Pi0_08_side);
ptmass_v_side.push_back(PT_Mass_Pi0_09_side);
ptmass_v_side.push_back(PT_Mass_Pi0_10_side);

TH1* sample_weights = data._sample_weights;
//**********************************************************************************************************
TH1* dummy_hist = new TH1D("dummy_hist", "dummy", 1, 0., 1.);

   std::vector<TLegend*> PT_legend;
   for (int pt_val : PT_BOUNDS){
      PT_legend.push_back(new TLegend(0.71, 0.77, .98, .98));
   }
   TLine* horz_line = new TLine(0, 1, 5, 1);
    
   //PT PROJECTIONS
   std::vector<TH1*> MassPi0_PT_proj, MassPi0_PT_proj_side;
   std::vector<TCanvas*> PT_canvas_v;
   std::vector<TH1*> ratio_plots;
  

   if (draw) TCanvas* all_canvas = new TCanvas("all_c", "Pt_bins");
   for (int pt_index = 0; pt_index < PT_BOUNDS.size(); pt_index++){

       std::ostringstream name_strs;
       name_strs << pt_index;

      TString name = tag + TString("_MassPi0_PT_") + TString(name_strs.str());
      TString title = TString("MassPi0_PT_") + TString(getPTRange(pt_index, PT_BOUNDS1)) + TString("_w/weights");

      MassPi0_PT_proj.push_back(new TH1D(name, title, 30, 0.0, 5.0));
      MassPi0_PT_proj_side.push_back(new TH1D(name + TString("_side"), title, 30, 0.0, 5.0));
      if (draw) PT_canvas_v.push_back(new TCanvas(getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]), getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]) + "side" ));
      if (draw) PT_canvas_v[pt_index]->Divide(1, 2);
      if (draw) PT_canvas_v[pt_index]->cd(1);

      TH2* plot = ptmass_v_side[pt_index];  //////////////////////////DEFINE THIS
      TH1* out_plot = MassPi0_PT_proj_side[pt_index];

      for (int x_bin = 1; x_bin < PT_Mass_Pi0_01->GetNbinsX(); x_bin++){

         double bin_val = 0.;
         double sum_up_weights = 0.;
         double sum_lo_weights = 0.;

         for (int w_bin = 1; w_bin <= PT_Mass_Pi0_01->GetNbinsY(); w_bin++){
         	double bin_weight = sample_weights->GetBinContent(w_bin);
         	double tbin_val = plot->GetBinContent(x_bin, w_bin);
      //   	if (bin_weight != 0.)
      //   	{
         //		std::cout << "N EVENTS: " << plot->GetBinContent(x_bin, w_bin) << std::endl;
         //		std::cout << "bin weight: " << bin_weight << std::endl;

       //  	} 
         	if (tbin_val < 1) continue;

            bin_val += tbin_val * bin_weight ;
            dummy_hist->SetBinContent(1, tbin_val);
            dummy_hist->SetBinErrorOption(TH1::kPoisson);

            if (bin_val > 0.0){


           		// std::cout << "BIN ERROR UP: " << dummy_hist->GetBinErrorUp(1) << std::endl;
           		// std::cout << "BIN val     : " << dummy_hist->GetBinContent(1) << std::endl;
            }

            sum_up_weights += dummy_hist->GetBinErrorUp(1) * dummy_hist->GetBinErrorUp(1) * bin_weight * bin_weight;
            sum_lo_weights += dummy_hist->GetBinErrorLow(1) * dummy_hist->GetBinErrorLow(1) * bin_weight * bin_weight;

         }

         sum_up_weights = TMath::Sqrt(sum_up_weights);
         sum_lo_weights = TMath::Sqrt(sum_lo_weights);

         //std::cout << sum_up_weights << std::endl;

        out_plot->SetBinContent(x_bin, bin_val);
        out_plot->SetBinError(x_bin, 0.50 * (sum_up_weights + sum_lo_weights));

      }

      plot = ptmass_v[pt_index];
      out_plot = MassPi0_PT_proj[pt_index];

      for (int x_bin = 1; x_bin < PT_Mass_Pi0_01->GetNbinsX(); x_bin++){

         double bin_val = 0.;
         double sum_up_weights = 0.;
         double sum_lo_weights = 0.;

         for (int w_bin = 1; w_bin <= PT_Mass_Pi0_01->GetNbinsY(); w_bin++){
         	double bin_weight = sample_weights->GetBinContent(w_bin);
         	double tbin_val = plot->GetBinContent(x_bin, w_bin);
      //   	if (bin_weight != 0.)
      //   	{
         //		std::cout << "N EVENTS: " << plot->GetBinContent(x_bin, w_bin) << std::endl;
         //		std::cout << "bin weight: " << bin_weight << std::endl;

       //  	} 
         	if (tbin_val < 1) continue;

            bin_val += tbin_val * bin_weight ;
            dummy_hist->SetBinContent(1, tbin_val);
            dummy_hist->SetBinErrorOption(TH1::kPoisson);

            if (bin_val > 0.0){


           		// std::cout << "BIN ERROR UP: " << dummy_hist->GetBinErrorUp(1) << std::endl;
           		// std::cout << "BIN val     : " << dummy_hist->GetBinContent(1) << std::endl;
            }

            sum_up_weights += dummy_hist->GetBinErrorUp(1) * dummy_hist->GetBinErrorUp(1) * bin_weight * bin_weight;
            sum_lo_weights += dummy_hist->GetBinErrorLow(1) * dummy_hist->GetBinErrorLow(1) * bin_weight * bin_weight;

         }

         sum_up_weights = TMath::Sqrt(sum_up_weights);
         sum_lo_weights = TMath::Sqrt(sum_lo_weights);

         std::cout << sum_up_weights << std::endl;

         if (bin_val > 0.0) out_plot->SetBinContent(x_bin, bin_val);
         if (bin_val > 0.0) out_plot->SetBinError(x_bin, 0.50 * (sum_up_weights + sum_lo_weights));

      }


      double sig_int_val = MassPi0_PT_proj[pt_index]->Integral();
      double side_int_val = MassPi0_PT_proj_side[pt_index]->Integral();

      FgScale(MassPi0_PT_proj_side[pt_index]);
   // std::cout << "\n\n\n******************" << std::endl;
      FgScale(MassPi0_PT_proj[pt_index]);

   // for (int j = 1; j < MassPi0_PT_proj[pt_index]->GetNbinsX(); j++){
   //    std::cout << j << ": " << "err: " << MassPi0_PT_proj[pt_index]->GetBinError(j)  <<  std::endl;
   // }
   // std::cout << "******************\n\n\n" << std::endl;


      MassPi0_PT_proj[pt_index]->SetLineColor(kBlack);
      MassPi0_PT_proj_side[pt_index]->SetLineColor(kBlue);
      //MassPi0_PT_proj[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj[pt_index]->SetMarkerColorAlpha(kBlack, 0);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerColorAlpha(kBlue, 0);

      MassPi0_PT_proj[pt_index]->SetTitle("TwoProng_MassPi0");
      
      if (draw) MassPi0_PT_proj[pt_index]->Draw();
      if (draw) MassPi0_PT_proj_side[pt_index]->Draw("SAME");
      
      TString sig_int = TString(getSciNotation(sig_int_val));
      TString side_int = TString(getSciNotation(side_int_val));

      PT_legend[pt_index]->SetHeader(MassPi0_PT_proj[pt_index]->GetName() ,"C");
      PT_legend[pt_index]->AddEntry(MassPi0_PT_proj[pt_index], TString("signal, integral: ") + sig_int );
      PT_legend[pt_index]->AddEntry(MassPi0_PT_proj_side[pt_index], TString("sideband, integral: ") + side_int );
      if (draw) PT_legend[pt_index]->Draw("SAME");

      if (draw)PT_canvas_v[pt_index]->cd(2);
      TH1* ratio_plot = getRatios(MassPi0_PT_proj[pt_index], MassPi0_PT_proj_side[pt_index], draw, tag);
      if (ratio_plot != nullptr){
         ratio_plot->SetName(tag + TString("_MassPi0_PT_") + TString(name_strs.str()) + TString("_ratio"));
         ratio_plots.push_back(ratio_plot);
      } else {
         ratio_plots.push_back(nullptr);
      }
      if (ratio_plot == nullptr) std::cout << " boop " << std::endl;
      if (draw and ratio_plots[pt_index] != nullptr) ratio_plots[pt_index]->Draw();
      if (draw) horz_line->Draw("SAME");

      //all_canvas->cd();
      //if (pt_index == 0) MassPi0_PT_proj[pt_index]->Draw();
      //else MassPi0_PT_proj[pt_index]->Draw("SAME");
   }

   //all_canvas->cd();
   //PT_legend->Draw("SAME");

   return ratio_plots;

}


template <typename Data_Set>
std::vector<TH1*> GenerateAnalysis1(Data_Set data, bool draw = false){
   TH1::SetDefaultSumw2();
   std::vector<int> PT_BOUNDS = user_input::PT_BOUNDS1;
   TString tag = data._tag + TString("_r");

   TH2* fTwoProng_MassPi0_PT = data._PT_Mass_Pi0_pt;
   TH2* fTwoProng_MassPi0_PT_side = data._PT_Mass_Pi0_pt_sideband;

   std::vector<TLegend*> PT_legend;
   for (int pt_val : PT_BOUNDS){
      PT_legend.push_back(new TLegend(0.71, 0.77, .98, .98));
   }
   TLine* horz_line = new TLine(0, 1, 5, 1);
    

   //PT PROJECTIONS
   std::vector<TH1*> MassPi0_PT_proj, MassPi0_PT_proj_side;
   std::vector<TCanvas*> PT_canvas_v;
   std::vector<TH1*> ratio_plots;
   if (draw) TCanvas* all_canvas = new TCanvas("all_c", "Pt_bins");
   for (int pt_index = 0; pt_index < PT_BOUNDS.size(); pt_index++){

      std::ostringstream name_strs;
      name_strs << pt_index;

      MassPi0_PT_proj.push_back(fTwoProng_MassPi0_PT->ProjectionX(getName(TString("PT_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]), PT_BOUNDS[pt_index] / 5, PT_BOUNDS[pt_index + 1] / 5, "e"));
      
      MassPi0_PT_proj_side.push_back(fTwoProng_MassPi0_PT_side->ProjectionX(getName(TString("PT_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]) + TString("side"), PT_BOUNDS[pt_index] / 5, PT_BOUNDS[pt_index + 1] / 5, "e"));
      if (draw) PT_canvas_v.push_back(new TCanvas(getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]), getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]) + "side" ));
      if (draw) PT_canvas_v[pt_index]->Divide(1, 2);
      if (draw) PT_canvas_v[pt_index]->cd(1);

      double sig_int_val = MassPi0_PT_proj[pt_index]->Integral();
      double side_int_val = MassPi0_PT_proj_side[pt_index]->Integral();

      Scale(MassPi0_PT_proj_side[pt_index]);
   // std::cout << "\n\n\n******************" << std::endl;
      Scale(MassPi0_PT_proj[pt_index]);

   // for (int j = 1; j < MassPi0_PT_proj[pt_index]->GetNbinsX(); j++){
   //    std::cout << j << ": " << "err: " << MassPi0_PT_proj[pt_index]->GetBinError(j)  <<  std::endl;
   // }
   // std::cout << "******************\n\n\n" << std::endl;

      MassPi0_PT_proj[pt_index]->SetLineColor(kBlack);
      MassPi0_PT_proj_side[pt_index]->SetLineColor(kBlue);
      //MassPi0_PT_proj[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj[pt_index]->SetMarkerColorAlpha(kBlack, 0);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerColorAlpha(kBlue, 0);

      MassPi0_PT_proj[pt_index]->SetTitle("TwoProng_MassPi0");
      
      if (draw) MassPi0_PT_proj[pt_index]->Draw();
      if (draw) MassPi0_PT_proj_side[pt_index]->Draw("SAME");
      
      TString sig_int = TString(getSciNotation(sig_int_val));
      TString side_int = TString(getSciNotation(side_int_val));

      PT_legend[pt_index]->SetHeader(MassPi0_PT_proj[pt_index]->GetName() ,"C");
      PT_legend[pt_index]->AddEntry(MassPi0_PT_proj[pt_index], TString("signal, integral: ") + sig_int );
      PT_legend[pt_index]->AddEntry(MassPi0_PT_proj_side[pt_index], TString("sideband, integral: ") + side_int );
      if (draw) PT_legend[pt_index]->Draw("SAME");

      if (draw)PT_canvas_v[pt_index]->cd(2);
      TH1* ratio_plot = getRatios(MassPi0_PT_proj[pt_index], MassPi0_PT_proj_side[pt_index], draw, tag);
      if (ratio_plot != nullptr){
         ratio_plot->SetName(tag + TString("_MassPi0_PT_") + TString(name_strs.str()) + TString("_ratio"));
         ratio_plots.push_back(ratio_plot);
      } else {
         ratio_plots.push_back(nullptr);
      }
      if (ratio_plot == nullptr) std::cout << " boop " << std::endl;
      if (draw and ratio_plots[pt_index] != nullptr) ratio_plots[pt_index]->Draw();
      if (draw) horz_line->Draw("SAME");

      //all_canvas->cd();
      //if (pt_index == 0) MassPi0_PT_proj[pt_index]->Draw();
      //else MassPi0_PT_proj[pt_index]->Draw("SAME");
   }

   //all_canvas->cd();
   //PT_legend->Draw("SAME");



   return ratio_plots;

}

template <typename Data_Set>
std::vector<std::vector<TH1*>> GetMassDistributions(Data_Set data, bool draw = false){

   std::vector<int> PT_BOUNDS = user_input::PT_BOUNDS1;
   TString tag = data._tag;

//GETTING DATA FROM OUTPUT FILE **********************************************************************

TH2* PT_Mass_Pi0_01 = data._PT_Mass_Pi0_01;
TH2* PT_Mass_Pi0_02 = data._PT_Mass_Pi0_02;
TH2* PT_Mass_Pi0_03 = data._PT_Mass_Pi0_03;
TH2* PT_Mass_Pi0_04 = data._PT_Mass_Pi0_04;
TH2* PT_Mass_Pi0_05 = data._PT_Mass_Pi0_05;
TH2* PT_Mass_Pi0_06 = data._PT_Mass_Pi0_06;
TH2* PT_Mass_Pi0_07 = data._PT_Mass_Pi0_07;
TH2* PT_Mass_Pi0_08 = data._PT_Mass_Pi0_08;
TH2* PT_Mass_Pi0_09 = data._PT_Mass_Pi0_09;
TH2* PT_Mass_Pi0_10 = data._PT_Mass_Pi0_10;


TH2* PT_Mass_Pi0_01_side = data._PT_Mass_Pi0_01_side;
TH2* PT_Mass_Pi0_02_side = data._PT_Mass_Pi0_02_side;
TH2* PT_Mass_Pi0_03_side = data._PT_Mass_Pi0_03_side;
TH2* PT_Mass_Pi0_04_side = data._PT_Mass_Pi0_04_side;
TH2* PT_Mass_Pi0_05_side = data._PT_Mass_Pi0_05_side;
TH2* PT_Mass_Pi0_06_side = data._PT_Mass_Pi0_06_side;
TH2* PT_Mass_Pi0_07_side = data._PT_Mass_Pi0_07_side;
TH2* PT_Mass_Pi0_08_side = data._PT_Mass_Pi0_08_side;
TH2* PT_Mass_Pi0_09_side = data._PT_Mass_Pi0_09_side;
TH2* PT_Mass_Pi0_10_side = data._PT_Mass_Pi0_10_side;

std::vector<TH2*> ptmass_v, ptmass_v_side;

ptmass_v.push_back(PT_Mass_Pi0_01);
ptmass_v.push_back(PT_Mass_Pi0_02);
ptmass_v.push_back(PT_Mass_Pi0_03);
ptmass_v.push_back(PT_Mass_Pi0_04);
ptmass_v.push_back(PT_Mass_Pi0_05);
ptmass_v.push_back(PT_Mass_Pi0_06);
ptmass_v.push_back(PT_Mass_Pi0_07);
ptmass_v.push_back(PT_Mass_Pi0_08);
ptmass_v.push_back(PT_Mass_Pi0_09);
ptmass_v.push_back(PT_Mass_Pi0_10);

ptmass_v_side.push_back(PT_Mass_Pi0_01_side);
ptmass_v_side.push_back(PT_Mass_Pi0_02_side);
ptmass_v_side.push_back(PT_Mass_Pi0_03_side);
ptmass_v_side.push_back(PT_Mass_Pi0_04_side);
ptmass_v_side.push_back(PT_Mass_Pi0_05_side);
ptmass_v_side.push_back(PT_Mass_Pi0_06_side);
ptmass_v_side.push_back(PT_Mass_Pi0_07_side);
ptmass_v_side.push_back(PT_Mass_Pi0_08_side);
ptmass_v_side.push_back(PT_Mass_Pi0_09_side);
ptmass_v_side.push_back(PT_Mass_Pi0_10_side);

TH1* sample_weights = data._sample_weights;
//**********************************************************************************************************
TH1* dummy_hist = new TH1D("dummy_hist", "dummy", 1, 0., 1.);

   std::vector<TLegend*> PT_legend;
   for (int pt_val : PT_BOUNDS){
      PT_legend.push_back(new TLegend(0.71, 0.77, .98, .98));
   }
   TLine* horz_line = new TLine(0, 1, 5, 1);
    
   //PT PROJECTIONS
   std::vector<TH1*> MassPi0_PT_proj, MassPi0_PT_proj_side;
   std::vector<TCanvas*> PT_canvas_v;
   std::vector<TH1*> ratio_plots;


   if (draw) TCanvas* all_canvas = new TCanvas("all_c", "Pt_bins");
   for (int pt_index = 0; pt_index < PT_BOUNDS.size(); pt_index++){

       std::ostringstream name_strs;
       name_strs << pt_index;

      TString name = tag + TString("_MassPi0_PT_") + TString(name_strs.str());
      TString title = TString("MassPi0_PT_") + TString(getPTRange(pt_index, PT_BOUNDS1)) + TString("_w/weights");

      MassPi0_PT_proj.push_back(new TH1D(name, title, 30, 0.0, 5.0));
      MassPi0_PT_proj_side.push_back(new TH1D(name + TString("_side"), title, 30, 0.0, 5.0));
      if (draw) PT_canvas_v.push_back(new TCanvas(getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]), getName(TString("canvas_"), PT_BOUNDS[pt_index], PT_BOUNDS[pt_index + 1]) + "side" ));
      if (draw) PT_canvas_v[pt_index]->Divide(1, 2);
      if (draw) PT_canvas_v[pt_index]->cd(1);

      TH2* plot = ptmass_v_side[pt_index];  //////////////////////////DEFINE THIS
      TH1* out_plot = MassPi0_PT_proj_side[pt_index];

      for (int x_bin = 1; x_bin < PT_Mass_Pi0_01->GetNbinsX(); x_bin++){

         double bin_val = 0.;
         double sum_up_weights = 0.;
         double sum_lo_weights = 0.;

         for (int w_bin = 1; w_bin <= PT_Mass_Pi0_01->GetNbinsY(); w_bin++){
            double bin_weight = sample_weights->GetBinContent(w_bin);
            double tbin_val = plot->GetBinContent(x_bin, w_bin);
      //    if (bin_weight != 0.)
      //    {
         //    std::cout << "N EVENTS: " << plot->GetBinContent(x_bin, w_bin) << std::endl;
         //    std::cout << "bin weight: " << bin_weight << std::endl;

       //   } 
            if (tbin_val < 1) continue;

            bin_val += tbin_val * bin_weight ;
            dummy_hist->SetBinContent(1, tbin_val);
            dummy_hist->SetBinErrorOption(TH1::kPoisson);

            if (bin_val > 0.0){


               // std::cout << "BIN ERROR UP: " << dummy_hist->GetBinErrorUp(1) << std::endl;
               // std::cout << "BIN val     : " << dummy_hist->GetBinContent(1) << std::endl;
            }

            sum_up_weights += dummy_hist->GetBinErrorUp(1) * dummy_hist->GetBinErrorUp(1) * bin_weight * bin_weight;
            sum_lo_weights += dummy_hist->GetBinErrorLow(1) * dummy_hist->GetBinErrorLow(1) * bin_weight * bin_weight;

         }

         sum_up_weights = TMath::Sqrt(sum_up_weights);
         sum_lo_weights = TMath::Sqrt(sum_lo_weights);

         //std::cout << sum_up_weights << std::endl;

        out_plot->SetBinContent(x_bin, bin_val);
        out_plot->SetBinError(x_bin, 0.50 * (sum_up_weights + sum_lo_weights));

      }

      plot = ptmass_v[pt_index];
      out_plot = MassPi0_PT_proj[pt_index];

      for (int x_bin = 1; x_bin < PT_Mass_Pi0_01->GetNbinsX(); x_bin++){

         double bin_val = 0.;
         double sum_up_weights = 0.;
         double sum_lo_weights = 0.;

         for (int w_bin = 1; w_bin <= PT_Mass_Pi0_01->GetNbinsY(); w_bin++){
            double bin_weight = sample_weights->GetBinContent(w_bin);
            double tbin_val = plot->GetBinContent(x_bin, w_bin);
      //    if (bin_weight != 0.)
      //    {
         //    std::cout << "N EVENTS: " << plot->GetBinContent(x_bin, w_bin) << std::endl;
         //    std::cout << "bin weight: " << bin_weight << std::endl;

       //   } 
            if (tbin_val < 1) continue;

            bin_val += tbin_val * bin_weight ;
            dummy_hist->SetBinContent(1, tbin_val);
            dummy_hist->SetBinErrorOption(TH1::kPoisson);

            if (bin_val > 0.0){


               // std::cout << "BIN ERROR UP: " << dummy_hist->GetBinErrorUp(1) << std::endl;
               // std::cout << "BIN val     : " << dummy_hist->GetBinContent(1) << std::endl;
            }

            sum_up_weights += dummy_hist->GetBinErrorUp(1) * dummy_hist->GetBinErrorUp(1) * bin_weight * bin_weight;
            sum_lo_weights += dummy_hist->GetBinErrorLow(1) * dummy_hist->GetBinErrorLow(1) * bin_weight * bin_weight;

         }

         sum_up_weights = TMath::Sqrt(sum_up_weights);
         sum_lo_weights = TMath::Sqrt(sum_lo_weights);

         std::cout << sum_up_weights << std::endl;

         if (bin_val > 0.0) out_plot->SetBinContent(x_bin, bin_val);
         if (bin_val > 0.0) out_plot->SetBinError(x_bin, 0.50 * (sum_up_weights + sum_lo_weights));

      }


      double sig_int_val = MassPi0_PT_proj[pt_index]->Integral();
      double side_int_val = MassPi0_PT_proj_side[pt_index]->Integral();

      FgScale(MassPi0_PT_proj_side[pt_index]);
   // std::cout << "\n\n\n******************" << std::endl;
      FgScale(MassPi0_PT_proj[pt_index]);

   // for (int j = 1; j < MassPi0_PT_proj[pt_index]->GetNbinsX(); j++){
   //    std::cout << j << ": " << "err: " << MassPi0_PT_proj[pt_index]->GetBinError(j)  <<  std::endl;
   // }
   // std::cout << "******************\n\n\n" << std::endl;


      MassPi0_PT_proj[pt_index]->SetLineColor(kBlack);
      MassPi0_PT_proj_side[pt_index]->SetLineColor(kBlue);
      //MassPi0_PT_proj[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerStyle(kFullSquare);
      //MassPi0_PT_proj[pt_index]->SetMarkerColorAlpha(kBlack, 0);
      //MassPi0_PT_proj_side[pt_index]->SetMarkerColorAlpha(kBlue, 0);

      MassPi0_PT_proj[pt_index]->SetTitle("TwoProng_MassPi0");
      
      if (draw) MassPi0_PT_proj[pt_index]->Draw();
      if (draw) MassPi0_PT_proj_side[pt_index]->Draw("SAME");
      
      
   }

   //all_canvas->cd();
   //PT_legend->Draw("SAME");
   std::vector<std::vector<TH1*>> slice_plots = {MassPi0_PT_proj, MassPi0_PT_proj_side};
   return slice_plots;

}


}; //NAMESPACE analysis