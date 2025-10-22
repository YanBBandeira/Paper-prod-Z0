#include <fstream>
void Y_dlogQ2(){
//====================================================================
  TCanvas *Y_dlogQ2 = new TCanvas("Y_dlogQ2", "Y_dlogQ2",163,81,780,624);
  gStyle->SetOptStat(0);
    Y_dlogQ2->SetFillColor(0);
    Y_dlogQ2->SetBorderMode(0);    
    Y_dlogQ2->SetBorderSize(0);
    Y_dlogQ2->SetLeftMargin(0.12);
    Y_dlogQ2->SetRightMargin(0.03);
    Y_dlogQ2->SetTopMargin(0.18);
    Y_dlogQ2->SetBottomMargin(0.18);
    Y_dlogQ2->SetFrameFillColor(0);
    Y_dlogQ2->SetFrameLineWidth(2);
    Y_dlogQ2->SetFrameBorderMode(0);
    Y_dlogQ2->SetLogx(0);
    Y_dlogQ2->SetLogy();
    Y_dlogQ2->SetGridx(0);
    Y_dlogQ2->SetGridy(0);
    Y_dlogQ2->SetTickx();
    Y_dlogQ2->SetTicky();

   ifstream plik;
   float ds, ptsum, ds_sse, ptsum_sse;
//=====================================================================
// read data from files
//=====================================================================   
//--------------------------------------------
   plik.open("YAIE_logQ2.dat");
   int N1=0;
   float x1[50000], y1[50000];
   while(!plik.eof()){
    plik >> ptsum >> ds;
    x1[N1] = ptsum;
    y1[N1] = ds;	
    N1++;
   };
   plik.close();

   plik.open("YFIE_logQ2.dat");
   int N2=0;
   float x2[50000], y2[50000];
   while(!plik.eof()){
    plik >> ptsum >> ds;
    x2[N2] = ptsum;
    y2[N2] = ds;	
    N2++;
   };
   plik.close();

   plik.open("YLIE_logQ2.dat");
   int N3=0;
   float x3[50000], y3[50000];
   while(!plik.eof()){
    plik >> ptsum_sse >> ds_sse;
    x3[N3] = ptsum_sse;
    y3[N3] = ds_sse;	
    N3++;
   };
   plik.close();
  
   plik.open("YAEI_logQ2.dat");
   int N4=0;
   float x4[50000], y4[50000];
   while(!plik.eof()){
    plik >> ptsum_sse >> ds_sse;
    x4[N4] = ptsum_sse;
    y4[N4] = ds_sse;	
    N4++;
   };
   plik.close();

   plik.open("YFEI_logQ2.dat");
   int N5=0;
   float x5[50000], y5[50000];
   while(!plik.eof()){
    plik >> ptsum_sse >> ds_sse;
    x5[N5] = ptsum_sse;
    y5[N5] = ds_sse;	
    N5++;
   };
   plik.close();   

   plik.open("YLEI_logQ2.dat");
   int N6=0;
   float x6[50000], y6[50000];
   while(!plik.eof()){
    plik >> ptsum_sse >> ds_sse;
    x6[N6] = ptsum_sse;
    y6[N6] = ds_sse;	
    N6++;
   };
   plik.close();
  

//=====================================================================
// graphs parameters
//=====================================================================    
//------------------------------------

   TGraph * gr1 = new TGraph(N1,x1,y1);   
   gr1->SetName("gr1");
   gr1->SetTitle("");
   gr1->SetLineColor(207);
   gr1->SetLineWidth(4);
   gr1->SetLineStyle(1);
   gr1->Draw("ALS");
   
   TGraph * gr2 = new TGraph(N2,x2,y2);   
   gr2->SetName("gr2");
   gr2->SetTitle("");
   gr2->SetLineColor(209);
   gr2->SetLineWidth(4);
   gr2->SetLineStyle(9);
   gr2->Draw("LS");
   
   TGraph * gr3 = new TGraph(N3,x3,y3);   
   gr3->SetName("gr3");
   gr3->SetTitle("");
   gr3->SetLineColor(57);
   gr3->SetLineWidth(4);
   gr3->SetLineStyle(7);
   gr3->Draw("LS");
  /*    
   TGraph * gr4 = new TGraph(N4,x4,y4);   
   gr4->SetName("gr4");
   gr4->SetTitle("");
   gr4->SetLineColor(213);
   gr4->SetLineWidth(4);
   gr4->SetLineStyle(1);
   gr4->Draw("LS"); 

   TGraph * gr5 = new TGraph(N5,x5,y5);   
   gr5->SetName("gr5");
   gr5->SetTitle("");
   gr5->SetLineColor(227);
   gr5->SetLineWidth(4);
   gr5->SetLineStyle(9);
   gr5->Draw("LS");

   TGraph * gr6 = new TGraph(N6,x6,y6);   
   gr6->SetName("gr6");
   gr6->SetTitle("");
   gr6->SetLineColor(216);
   gr6->SetLineWidth(4);
   gr6->SetLineStyle(5);
   gr6->Draw("LS");
    
   TGraph * gr7 = new TGraph(N7,x7,y7);   
   gr7->SetName("gr7");
   gr7->SetTitle("");
   gr7->SetLineColor(kBlack);
   gr7->SetLineWidth(4);
   gr7->SetLineStyle(1);
   gr7->Draw("LS");
   
   TGraph * gr8 = new TGraph(N8,x8,y8);   
   gr8->SetName("gr8");
   gr8->SetTitle("");
   gr8->SetLineColor(213);
   gr8->SetLineWidth(4);
   gr8->SetLineStyle(1);   
   gr8->Draw("LS");    
*/
//=====================================================================
// axys description
//=====================================================================

    TAxis *X = gr1->GetXaxis();
    X->SetTitle("log_{10}(Q_{1}^{2}) (GeV^{2})");
    X->CenterTitle();
    X->SetTitleFont(42);
    X->SetLabelFont(42);
    X->SetTitleSize(0.055);
    X->SetTitleOffset(1.2);
    X->SetRangeUser(-3.5,4.);

    TAxis *Y = gr1->GetYaxis();
    Y->SetTitle("d#sigma/dlog_{10}(Q_{1}^{2}) (pb/GeV^{2})");
    Y->CenterTitle();
    Y->SetTitleFont(42);
    Y->SetLabelFont(42);
    Y->SetTitleSize(0.055);
    Y->SetTitleOffset(1.);
    Y->SetRangeUser(1e-5,2);
  
//=====================================================================
// writing on the graph
//=====================================================================

// Basics information

   TLatex *tex1 = new TLatex(-3.2,6.5e-1,"pp #rightarrow #mu^{+}#mu^{-} p(X) p(Y)");
    tex1->SetTextColor(1);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->Draw();
         
   TLatex *tex2 = new TLatex(-3.24,2.5e-1,"#sqrt{s}=13 TeV");
    tex2->SetTextColor(1);
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->Draw();

// Additional information
 
   TLatex *tex3 = new TLatex(2.05,7.5e-1,"p_{#perp 1} , p_{#perp 2} < 15 GeV");
    tex3->SetTextColor(kRed+3);
    tex3->SetTextFont(42);
    tex3->SetTextSize(0.036);
    tex3->Draw();

   TLatex *tex4 = new TLatex(2.24,2.3e-1,"-2.5 < Y_{ll} < 2.5");
    tex4->SetTextColor(kRed+3);
    tex4->SetTextFont(42);
    tex4->SetTextSize(0.036);
    tex4->Draw();    

   TLatex *tex5 = new TLatex(2.14,1e-1,"inelastic-elastic");
    tex5->SetTextColor(kRed+3);
    tex5->SetTextFont(42);
    tex5->SetTextSize(0.036);
    tex5->Draw();      
    
//=====================================================================
// legend
//=====================================================================
         
   TLegend *leg = new TLegend(0.12,0.82,0.97,0.98);
   leg->SetNColumns(3);
   leg->SetFillColor(0);
   leg->SetLineColor(1);
   leg->SetLineWidth(2);
   leg->SetTextSize(0.033);
   leg->AddEntry("gr1","ALLM F_{2}","l");
//   leg->AddEntry("gr4","ALLM elastic-inelastic","l");
   leg->AddEntry("gr2","FFJLM F_{2}","l");   
//   leg->AddEntry("gr5","FFJLM elastic-inelastic","l");
   leg->AddEntry("gr3","LUX-like F_{2}+F_{L}","l");
//   leg->AddEntry("gr6","LUX-like elastic-inelastic","l");
   leg->Draw();

//=====================================================================
// writing to file
//===================================================================== 

   Y_dlogQ2->Update(); 
   Y_dlogQ2->cd();
   Y->Draw();
   X->Draw();
   Y_dlogQ2->Print("ppl_dsig_dlogQ2Y.eps");
   Y_dlogQ2->Print("ppl_dsig_dlogQ2Y.pdf");
}

