#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


void count_entries(const char* listname="list.txt", const char* treename="analyze/Analysis"){
  TFile *f;
  TTree *tree;
  Int_t count=0;
  Int_t filecount=0;

ifstream list;
list.open(listname);
char filename[300];

while (true){
list >> filename;
if (!list.good()) break;
cout << filename << endl;
 f=TFile::Open(filename,"read");
 filecount++;
 f->GetObject(treename,tree);
 count+=tree->GetEntries();
 f->Close();
}

list.close();

 cout << "Total entries: " << count << " from " << filecount << " files" <<  endl;

}
