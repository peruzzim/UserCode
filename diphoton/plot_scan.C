{
  TFile *f = new TFile("outphoton_ggjets_bkg_rconescan.root");
  f->cd("scan_cone");

  scan_cone_EB_0p450->SetMarkerStyle(20);
  scan_cone_EB_0p450->SetMarkerColor(kRed);
  scan_cone_EB_0p450->Draw("P");
  scan_cone_EB_0p500->SetMarkerStyle(20);
  scan_cone_EB_0p500->SetMarkerColor(kMagenta);
  scan_cone_EB_0p500->Draw("Psame");
  scan_cone_EB_0p800->SetMarkerStyle(20);
  scan_cone_EB_0p800->SetMarkerColor(kBlue);
  scan_cone_EB_0p800->Draw("Psame");
  scan_cone_EB_0p1000->SetMarkerStyle(20);
  scan_cone_EB_0p1000->SetMarkerColor(kGreen);
  scan_cone_EB_0p1000->Draw("Psame");
  scan_cone_EB_0p1200->SetMarkerStyle(20);
  scan_cone_EB_0p1200->SetMarkerColor(kOrange);
  scan_cone_EB_0p1200->Draw("Psame");




}
