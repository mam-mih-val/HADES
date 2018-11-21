void reader()
{
    auto fStyle = new TFile("style.root");
    auto style = (TStyle*) fStyle->Get("style");
    gROOT->ForceStyle();
    style->cd();
    fStyle->Close();
    auto f = new TFile("AuAu.root");
    auto t = (TTree*) f->Get("DataTree");
    auto ev = new DataTreeEvent;
    auto DTEvent = (TBranch*) t->GetBranch("DTEvent");
    DTEvent->SetAddress(&ev);

    TH1F* multy_MDC = new TH1F("number of tracks, MDC","number of tracks MDC & hits TOF+RPC",100,1,100);
    multy_MDC->SetLineColor(1);
    TH1F* multy_TOF = new TH1F("number of hits, TOF","number of hits, TOF+RPC",100,1,100);
    multy_TOF->SetLineColor(2);
    Long64_t n_events = (Long64_t) t->GetEntries();
    Int_t n_tracks = 0;
    Int_t n_hits = 0;
    for(int i=0;i<n_events;i++)
    {
        DTEvent->GetEntry(i);
        n_tracks = ev->GetNVertexTracks();
        n_hits = ev->GetNTOFHits();
        multy_MDC->Fill(n_tracks);
        multy_TOF->Fill(n_hits);
    }
    multy_MDC->Draw();
    multy_TOF->Draw("same");
}