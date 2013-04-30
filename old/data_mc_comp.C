void data_mc_comp(int code, int cat){

const char* masspointsfile = "masspoints_ee_";


 TString pathmc("output_mc/");
 TString pathdata("output_data/");
    pathmc.Append(masspointsfile);
    pathdata.Append(masspointsfile);
    pathmc+=code;
    pathdata+=code;
    pathmc.Append("_cat");
    pathdata.Append("_cat");
    pathmc+=cat;
    pathdata+=cat;
    pathmc.Append(".txt");
    pathdata.Append(".txt");

    ifstream inmc;
    inmc.open(pathmc.Data());
    ifstream indata;
    indata.open(pathdata.Data());

    TH1F *histomc = new TH1F("histomc","histomc",100,80,100);
    TH1F *histodata = new TH1F("histodata","histodata",100,80,100);

    /*
    TH1F *histomc00 = new TH1F("histomc00","histomc00",100,80,100);
    TH1F *histodata00 = new TH1F("histodata00","histodata00",100,80,100);
    */

    while(1){
      float x;
      float w;
      inmc >> x;
      inmc >> w;
      if (!inmc.good()) break;
      histomc->Fill(x,w);
    }

   while(1){
      float x;
      float w;
      indata >> x;
      indata >> w;
      if (!indata.good()) break;
      histodata->Fill(x,w);
    }

   /*
   {
     ifstream in;
     in.open("output_mc/masspoints_ee_0_cat0.txt");
  while(1){
      float x;
      in >> x;
      if (!in.good()) break;
      histomc00->Fill(x);
    }
  in.close();
  in.open("output_data/masspoints_ee_0_cat0.txt");
   while(1){
      float x;
      in >> x;
      if (!in.good()) break;
      histodata00->Fill(x);
    }
   in.close();
   }
   */

   inmc.close();
   indata.close();

   histomc->SetLineColor(kRed);

   histomc->Scale(((float)histodata->Integral())/histomc->Integral());

   /*
   histomc00->SetLineColor(kRed);
   histomc00->Scale(((float)histodata00->GetEntries())/histomc00->GetEntries());

   histomc00->SetLineStyle(3);
   histodata00->SetLineStyle(3);
   histomc00->Draw();
   histodata00->Draw("same");
   */

   histomc->Draw();
   histodata->Draw("same");

}
