{
  for (int i=0; i<8; i++){
    TCanvas *c = new TCanvas();
    c->cd();
    (RooPlot*)(gDirectory->Get(Form("roovar_dR_obs_roodset_EBEB_invmass_b%d",i)))->Draw();
    c->SaveAs(Form("plots_dR_%d.pdf",i));
  }
}
