#include "reader.h"

reader::reader()
{
    fChain = new TChain("DataTree");
    fChain->Add("AuAu.root");
    ev = new DataTreeEvent;
    DTEvent = (TBranch*) fChain->GetBranch("DTEvent");
    DTEvent->SetAddress(&ev);
    MeanQx = 0; MeanQy = 0;
//    cout << "success reading" << endl;
}
reader::~reader()
{
    delete ev;
    delete fChain;
}
void reader::GetEvent(Int_t n)
{
    fChain->GetEvent(n);
}

void reader::Main_Loop()
{
    TProfile* pt_vs_v1 = new TProfile("y vs v_{1}","rapidity vs directed flow",8,-0.5,1.5);
    pt_vs_v1->GetXaxis()->SetTitle("rapidity, y");
    pt_vs_v1->GetYaxis()->SetTitle("v_{1}");
    this->GetQMean();
    Int_t N_ev = (Int_t) fChain->GetEntries();
    DataTreeTrack* tr;
    Float_t phi=0, pt=0, v1, y, w;
    Float_t* PsiEP = new Float_t[2];
    for(int i=0;i<N_ev;i++)
    {
        fChain->GetEntry(i);
        PsiEP = this->GetPsiEP();
        Int_t n_tracks = (Int_t) ev->GetNVertexTracks();
        for(int j=0;j<n_tracks;j++)
        {
            tr = ev->GetVertexTrack(j);
            phi = tr->GetPhi();
            pt = tr->GetPt();
            y = tr->GetRapidity();
            if(y>=0)
                w = 1;
            if(y<0)
                w = -1;
            v1 = this->Get_v(PsiEP,phi);
//            cout<< "v1=" << v1 << " phi=" << phi << " Psi=" << PsiEP[0] << endl;
            pt_vs_v1->Fill(y,v1);
        }
    }
    TCanvas* canv = new TCanvas("canv","y vs v_{1}",1500,700);
    canv->cd();
    pt_vs_v1->Draw();
    TLine* l = new TLine(-1,0,1,0);
    //l->Draw("same");
}

Float_t reader::Get_v( Float_t* PsiEP, Float_t phi, Int_t n=1  )
{
    Float_t MeanPsiEP = ( PsiEP[0] + PsiEP[1] ) / 2;
    Float_t v_n=0;
    v_n = cos( n*(phi-MeanPsiEP) )/cos( n*(PsiEP[0]-PsiEP[1]) );
    return v_n;
}

void reader::GetQMean() 
{
    Int_t N = fChain->GetEntries();
    Int_t n_psd_modules;
    DataTreePSDModule* psd_module;
    Float_t q, phi;
    Float_t Qx, Qy;
    Float_t EofAll;
    TH2F* histo_Q = new TH2F("Q-vectors distribution","Q-vectors distribution",100,-2.5,2.5,100,-2.5,2.5);
    histo_Q->GetXaxis()->SetTitle("Qx");
    histo_Q->GetYaxis()->SetTitle("Qy");
    for(Long64_t i=0; i<N; i++) 
    {
        Qx=0; Qy=0;
        fChain->GetEntry(i);
        EofAll = ev->GetPSDEnergy();
        n_psd_modules = ev->GetNPSDModules();
        for(int j=0;j<n_psd_modules;j++)
        {
            psd_module = ev->GetPSDModule(j);
            phi =        psd_module->GetPhi();
            q   =        psd_module->GetEnergy();
            Qx+=q*cos(phi); Qy+=q*sin(phi);
        }
        histo_Q->Fill(Qx,Qy);
    }
    MeanQx = histo_Q->GetMean(1);
    MeanQy = histo_Q->GetMean(2);
//    delete histo_Q;
}

Float_t* reader::GetPsiEP() 
{
    Int_t n_psd_modules;
    DataTreePSDModule* psd_module;
    Float_t q, phi;
    Float_t Qx[2], Qy[2];
    Float_t EofAll;

    for(int i=0;i<2;i++)
    {
        Qx[0] = 0;
        Qy[0] = 0;
    }
    EofAll = ev->GetPSDEnergy();
    n_psd_modules = ev->GetNPSDModules();
    for(int j=0;j<n_psd_modules;j++)
    {
        psd_module = ev->GetPSDModule(j);
        phi =        psd_module->GetPhi();
        q   =        psd_module->GetEnergy();
        Int_t p = gRandom->Uniform(0,1);
        Qx[p]+= q*cos(phi); 
        Qy[p]+= q*sin(phi);
    }
    Float_t* PsiEP = new Float_t[2];
    for(int j=0;j<2;j++)
    {
        Qx[j]-=MeanQx; 
        Qy[j]-=MeanQy;
        PsiEP[j] = atan2(Qy[j],Qx[j]);
//        cout << "Qx=" << Qx[j] << " Qy=" << Qy[j] << " PsiEP=" << PsiEP[j] << endl;
    }
    return PsiEP;
}
Bool_t  reader::evSelector()
{
    if (ev->GetNVertexTracks() < 5)
        return 0;
    if (ev->GetNTOFHits() < 5)
        return 0;
    return 1;
}