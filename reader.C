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

    TH1F* multy_MDC = new TH1F("number of tracks, MDC","number of tracks MDC",100,1,100);
    multy_MDC->GetXaxis()->SetTitle("MDC tracks");

    TH1F* multy_TOF = new TH1F("number of hits, TOF","number of hits, TOF+RPC",100,1,100);
    multy_TOF->GetXaxis()->SetTitle("TOF hits");

    TH2F* q_vs_tr = new TH2F("charge vs tracks","charge in FW vs number of tracks",100,1,100,100,1,8000);
    q_vs_tr->GetXaxis()->SetTitle("MDC tracks");
    q_vs_tr->GetYaxis()->SetTitle("FW charge");

    TH2F* q_vs_h = new TH2F("charge vs hits","charge in FW vs number of hits",100,1,100,100,1,8000);
    q_vs_h->GetXaxis()->SetTitle("TOF hits");
    q_vs_h->GetYaxis()->SetTitle("FW charge");

    TH2F* vertex = new TH2F("vertex position","vertex position",20,-50,50,20,-150,50);
    vertex->GetXaxis()->SetTitle("x");
    vertex->GetYaxis()->SetTitle("z");

    TH1F* histo_pt = new TH1F("pt distribution","pt distribution",100,0.0,2.5);
    histo_pt->GetXaxis()->SetTitle("pt, [GeV/c]");

    Long64_t n_events = (Long64_t) t->GetEntries();
    Int_t n_tracks = 0;
    Int_t n_hits = 0;
    Double_t charge = 0;
    Double_t vert_pos[3];
    Double_t x1,x2;
    DataTreeTrack* track;
    Double_t pt=0;
    for(int i=0;i<n_events;i++)
    {
        DTEvent->GetEntry(i);
        n_tracks = ev->GetNVertexTracks();
        n_hits = ev->GetNTOFHits();
        charge = ev->GetPSDEnergy();

        multy_MDC->Fill(n_tracks);
        multy_TOF->Fill(n_hits);
        q_vs_tr->Fill(n_tracks,charge);
        q_vs_h->Fill(n_hits,charge);
        for(int j=0;j<3;j++)
        {
            vert_pos[j]=ev->GetVertexPositionComponent(j);
        }
        vertex->Fill(vert_pos[0],vert_pos[2]);
        for(int j=0;j<n_tracks;j++)
        {
            track = ev->GetVertexTrack(j);
            pt = track->GetPt();
            histo_pt->Fill(pt);
        }
    }

    histo_pt->Draw();
    TFile* w = new TFile("histo.root","recreate");
    w->cd();
    multy_MDC->Write();
    multy_TOF->Write();
    q_vs_tr->Write();
    q_vs_h->Write();
    vertex->Write();
    histo_pt->Write();
    style->Write();
    w->Close();
}