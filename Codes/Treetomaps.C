
//====================================================================
//      This program generates 3D maps from file.root
//      calculations from the Lorentz vector
//====================================================================

#define max_events 3000000   //chagne!!!!!!
#define MX         0
#define MY         1
#define mll        2
#define t1         3
#define t2         4
#define logt1      5
#define logt2      6
#define pt1        7
#define pt2        8
#define logW12     9
#define logW22     10
#define logxBj1    11
#define logxBj2    12

#define max_dists  13        //change!!!!!!
int lwyk =         14;       //chagne!!!!!!

void
Treetomaps()
{
 TFile *f_pptoll[lwyk];
 TTree *t_pptoll[lwyk];
  
 TLorentzVector l1, l2, op1, op2, ip1, ip2; 
 TLorentzVector pJet1, pJet2;
 Double_t       weight, x1, x2;
 TLine *line;
 TString        tmp;
 Double_t xBj1, xBj2, t1val, t2val;
 Double_t mp1,mp2,MX2, MY2, W12, W1;


 const Int_t num_2d_distributions = lwyk;
 TCanvas *c2[num_2d_distributions];
 TH2D *h2[num_2d_distributions];
 gStyle->SetOptStat(0);
 gStyle->SetPadRightMargin(0.2);
 gStyle->SetPadLeftMargin(0.2);

//-------------------------------------------------------------------
//          number of points and range of axes for each variable
//-------------------------------------------------------------------
 f_pptoll[0] = new TFile("events.root");

 tmp.Form("MX_MY");   h2[0] = new TH2D("MX_MY",  "LUX-like F_{2}", 100, 0., 1000., 100, 0., 250.);
 tmp.Form("t1_t2");   h2[1] = new TH2D("t1_t2",  "LUX-like F_{2}", 100, 0., 1000., 100, 0., 0.016);
 tmp.Form("t1_MX");   h2[2] = new TH2D("t1_MX",  "LUX-like F_{2}", 100, 0., 25, 100, 0., 200.);
 tmp.Form("pt1_pt2"); h2[3] = new TH2D("pt1_pt2","LUX-like F_{2}", 100, 0.2, 10., 100, 0.2, 10.);
 tmp.Form("W12_t1");  h2[4] = new TH2D("W12_t1", "LUX-like F_{2}", 50, 0.0, 5., 50, 1., 4.); 
 tmp.Form("W1_t1");   h2[5] = new TH2D("W1_t1",  "LUX-like F_{2}", 50, 0.2, 5., 50, 1., 10.); 


 tmp.Form("logt1_MX");  h2[6] = new TH2D("logt1_MX",  "LUX-like F_{2}", 100, -5., 4., 100, 0., 400.);
 tmp.Form("logt1_MY");  h2[7] = new TH2D("logt1_MY",  "LUX-like F_{2}", 100, -6., 6., 100, 0., 1000.);
 tmp.Form("logt2_MX");  h2[8] = new TH2D("logt2_MX",  "LUX-like F_{2}", 100, -8., 0., 100, 0., 400.);
 tmp.Form("logt2_MY");  h2[9]=  new TH2D("logt2_MY",  "LUX-like F_{2}", 100, -6., 6., 100, 0., 1000.); 
 tmp.Form("logt1_Mll"); h2[10] = new TH2D("logt1_Mll","LUX-like F_{2}", 100, -3.0, 3.0, 100, 8., 50.);
 
 
 tmp.Form("logt1_logt2");   h2[11] = new TH2D("logt1_logt2",  "LUX-like F_{2}", 100, -4., 4., 100, -8., 0.);
 tmp.Form("logt1_logxBj1"); h2[12] = new TH2D("logt1_logxBj1","LUX-like F_{2}", 100, -4., 2., 100, -4.5, 0.);
 tmp.Form("logW12_logt1");  h2[13] = new TH2D("logW12_logt1", "LUX-like F_{2}", 100, 0.0, 5., 100, -4., 2.);


//-------------------------------------------------------------------
//           processing - Lorentz operator
//-------------------------------------------------------------------
  Double_t px[9], py[9], pz[9], E[9];
  Double_t xsect;
 
  cout << "Extracting data from .root file" << endl;

  t_pptoll[0] = (TTree*)(f_pptoll[0]->Get("events"));
  t_pptoll[0]->SetBranchAddress("xsect", &xsect);
  t_pptoll[0]->SetBranchAddress("px", px);
  t_pptoll[0]->SetBranchAddress("py", py);
  t_pptoll[0]->SetBranchAddress("pz", pz);
  t_pptoll[0]->SetBranchAddress("E",  E);

  cout << "Calculation cross setcion for each variables" << endl;
   
  for (Int_t i=0; i<max_events; i++) {
    t_pptoll[0]->GetEntry(i);
    weight = xsect/max_events;
    
    l1.SetXYZT (px[5], py[5], pz[5], E[5]);
    l2.SetXYZT (px[6], py[6], pz[6], E[6]);
    ip1.SetXYZT(px[0], py[0], pz[0], E[0]);
    ip2.SetXYZT(px[1], py[1], pz[1], E[1]);
    op1.SetXYZT(px[2], py[2], pz[2], E[2]);
    op2.SetXYZT(px[4], py[4], pz[4], E[4]);

//-------------------------------------------------------------------
//           definition of variables
//-------------------------------------------------------------------
   mp1 = ip1.M2();
   mp2 = ip2.M2();
   MX2 = op1.M2();                  //(MX1)^2
   MY2 = op2.M2();                  //(MX2)^2

   t1val= -(op1-ip1).M2();          //(q1)^2
   t2val= -(op2-ip2).M2();          //(q2)^2  
   
   xBj1 = t1val/(t1val + MX2- mp1);
   xBj2 = t2val/(t2val + MY2- mp2);
   
   W12 = (t1val/xBj1)-t1val+mp1;    //(W1)^2
   W1  = sqrt(W12);

//-------------------------------------------------------------------
//           calculation of the cross sections
//-------------------------------------------------------------------             
    h2[0]  ->Fill(op1.M(), op2.M(), weight);
    h2[1]  ->Fill(t1val,   t2val,   weight);
    h2[2]  ->Fill(t1val,   op1.M(), weight);
    h2[3]  ->Fill(l1.Pt(), l2.Pt(), weight);
    h2[4]  ->Fill(t1val,   W12,     weight);            
    h2[5]  ->Fill(t1val,   W1,      weight);  

    h2[6]  ->Fill(log10(t1val), op1.M(), weight);
    h2[7]  ->Fill(log10(t1val), op2.M(), weight);
    h2[8]  ->Fill(log10(t2val), op1.M(), weight);
    h2[9]  ->Fill(log10(t2val), op2.M(), weight);
    h2[10] ->Fill(log10(t1val), (l1+l2).M(), weight);

    h2[11] ->Fill(log10(t1val), log10(t2val), weight);
    h2[12] ->Fill(log10(t1val), log10(xBj1),  weight);
    h2[13] ->Fill(log10(W12),   log10(t1val), weight);
  }

//==================================================================
//           GENERATING GRAPHS
//==================================================================
 cout << "Generating graphs..." << endl;

//------------------------------------------------------------------
//           Axis names
//------------------------------------------------------------------
 TString xtit[num_2d_distributions], ytit[num_2d_distributions], ztit[num_2d_distributions]; 
  xtit[0] = "M_{X} (GeV)";       ytit[0] = "M_{Y} (GeV)";      
  ztit[0] = "d^{2}#sigma/dM_{X}dM_{Y} (pb/GeV^{2})";
  
  xtit[1] = "Q_{1}^{2} (GeV^{2})"; ytit[1] = "Q_{2}^{2} (GeV^{2})";
  ztit[1] = "d^{2}#sigma/dQ_{1}^{2}dQ_{2}^{2} (pb/GeV^{4})";
  
  xtit[2] = "Q_{1}^{2} (GeV^{2})"; ytit[2] = "M_{X} (GeV)";
  ztit[2] = "d^{2}#sigma/dQ_{1}^{2}dM_{X} (pb/GeV^{3})";
  
  xtit[3] = "p_{t1} (GeV)";      ytit[3] = "p_{t2} (GeV)";  
  ztit[3] = "d^{2}#sigma/dp_{t1}dp_{t2} (pb/GeV^{2})";
  
  xtit[4] = "Q_{1}^{2} (GeV^{2})"; ytit[4] = "W_{1}^{2} (GeV^{2})";     
  ztit[4] = "d^{2}#sigma/dQ_{1}^{2}dW_{1}^{2} (pb/GeV^{4})";
  
  xtit[5] = "Q_{1}^{2} (GeV^{2})"; ytit[5] = "W_{1} (GeV)"; 
  ztit[5] = "d^{2}#sigma/dQ_{1}^{2}dW_{1} (pb/GeV^{3})";


  xtit[6] = "log_{10}(Q_{1}^{2})"; ytit[6] = "M_{X} (GeV)";
  ztit[6] = "d^{2}#sigma/dlog_{10}(Q_{1}^{2})dM_{X} (pb/GeV)";
  
  xtit[7] = "log_{10}(Q_{1}^{2})"; ytit[7] = "M_{Y} (GeV)";
  ztit[7] = "d^{2}#sigma/dlog_{10}(Q_{1}^{2})dM_{Y} (pb/GeV)";
  
  xtit[8] = "log_{10}(Q_{2}^{2})"; ytit[8] = "M_{X} (GeV)";
  ztit[8] = "d^{2}#sigma/dlog_{10}(Q_{2}^{2})dM_{X} (pb/GeV)";
  
  xtit[9] = "log_{10}(Q_{2}^{2})"; ytit[9] = "M_{Y} (GeV)";
  ztit[9] = "d^{2}#sigma/dlog_{10}(Q_{2}^{2})dM_{Y} (pb/GeV)";
  
  xtit[10]= "log_{10}(Q_{1}^{2})"; ytit[10]= "M_{e^{+}e^{-}} (GeV)";  
  ztit[10] = "d^{2}#sigma/dlog_{10}(Q_{1}^{2})dM_{e^{+}e^{-}} (pb/GeV)";  


  xtit[11] = "log_{10}(Q_{1}^{2})"; ytit[11] = "log_{10}(Q_{2}^{2})";
  ztit[11] = "d^{2}#sigma/dlog_{10}(Q_{1}^{2}) dlog_{10}(Q_{2}^{2}) (pb)";  

  xtit[12] = "log_{10}(Q_{1}^{2})"; ytit[12] = "log_{10}(x_{Bj1})";
  ztit[12] = "d^{2}#sigma/dlog_{10}(Q_{1}^{2}) dlog_{10}(x_{Bj1}) (pb)";  

  xtit[13] = "log_{10}(W_{1}^{2})"; ytit[13] = "log_{10}(Q_{1}^{2})";
  ztit[13] = "d^{2}#sigma/dlog_{10}(W_{1}^{2}) dlog_{10}(Q_{1}^{2}) (pb)";  

  for (Int_t i=0; i<num_2d_distributions; i++) {
    c2[i] = new TCanvas(h2[i]->GetName(), "", 600, 600);
    h2[i]->Draw("colz");  
    h2[i]->GetXaxis()->SetTitle(xtit[i]);
    h2[i]->GetXaxis()->SetTitleOffset(1.2);
    h2[i]->GetYaxis()->SetTitle(ytit[i]);
    h2[i]->GetYaxis()->SetTitleOffset(2.);
    h2[i]->GetZaxis()->SetTitle(ztit[i]);
    h2[i]->GetZaxis()->SetTitleOffset(1.5); 
    h2[i]->SetMaximum(3000);

    c2[i]->SetLogz();
    c2[i]->SaveAs("L_IMR_"+TString(h2[i]->GetName())+".pdf");
//    c2[i]->SaveAs(TString(h2[i]->GetName())+"_colz.eps");

    h2[i]->Draw("lego2");
    gStyle->SetPalette();
    h2[i]->GetXaxis()->SetTitle(xtit[i]);
    h2[i]->GetXaxis()->SetTitleOffset(2.2);
    h2[i]->GetYaxis()->SetTitle(ytit[i]);
    h2[i]->GetYaxis()->SetTitleOffset(2.);
    h2[i]->GetZaxis()->SetTitle(ztit[i]);
    h2[i]->GetZaxis()->SetTitleOffset(1.85);    
    h2[i]->SetMaximum(3000.);  //max value of color scale
    h2[i]->SetMinimum(1.);    //min value of color scale
    h2[i]->SetLineWidth(0.);
//    h2[i]->SetLineColor(214);

    c2[i]->SetPhi(-120);
    c2[i]->SetLogz();
    c2[i]->SaveAs("L_IMR_"+TString(h2[i]->GetName())+"_lego.pdf");  
//    c2[i]->SaveAs(TString(h2[i]->GetName())+"_lego.eps");     
  }

}
