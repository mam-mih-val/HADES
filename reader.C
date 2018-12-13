#include "reader.h"

reader::reader()
{
    fChain = new TChain("DataTree");
    fChain->Add("AuAu.root");
    ev = new DataTreeEvent;
    DTEvent = (TBranch*) fChain->GetBranch("DTEvent");
    DTEvent->SetAddress(&ev);

    MeanQx = new TProfile("c vs Qx","centrality vs Qx",20,0,100);
    MeanQx->GetYaxis()->SetTitle("<Qx>");
    MeanQx->GetXaxis()->SetTitle("centrality");

    MeanQy = new TProfile("c vs Qy","centrality vs Qy",20,0,100);
    MeanQy->GetYaxis()->SetTitle("<Qy>");
    MeanQy->GetXaxis()->SetTitle("centrality");

    ResPsi = new TProfile("resolution of SE, r","Resolution of subevent method, recentred",20,0,100);
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

void reader::Q_histos()
{
    this->GetQMean();
    Float_t* Q = new Float_t[2];
    Int_t n_ev = fChain->GetEntries();
    TProfile* Qx_vs_c = new TProfile("Qx vs c, recentred", "Qx vs centrality after recentring",20,0,100);
    TProfile* Qy_vs_c = new TProfile("Qy vs c, recentred", "Qy vs centrality after recentring",20,0,100);
    TH1F* Qx_nr_histo = new TH1F("Qx nr","Qx not recented",100,-1,1);
    TH1F* Qx_r_histo = new TH1F("Qx r","Qx recented",100,-1,1);
    TH1F* Qy_nr_histo = new TH1F("Qy nr","Qy not recented",100,-1,1);
    TH1F* Qy_r_histo = new TH1F("Qy r","Qy recented",100,-1,1);
    Float_t cent;
    Float_t meanx, meany;
    Int_t bin;
    for(int i=0;i<n_ev;i++)
    {
        fChain->GetEntry(i);
        if(ev->GetPSDEnergy()==0)
            continue;
        cent = ev->GetCentrality();
        if(cent<0)
            continue;
        bin = cent/5+1;
        meanx = MeanQx->GetBinContent(bin);
        meany = MeanQy->GetBinContent(bin);
        Q = this->GetQ();
        Qx_nr_histo->Fill(Q[0]);
        Qx_r_histo->Fill(Q[0]-meanx);
        Qy_nr_histo->Fill(Q[1]);
        Qy_r_histo->Fill(Q[1]-meany);
        Qx_vs_c->Fill(cent,Q[0]-meanx);
        Qy_vs_c->Fill(cent,Q[1]-meany);
        this->ShowProgress(i,n_ev);
    }
    TFile* f = new TFile("histos.root","recreate");
    f->cd();
    MeanQx->Write();
    MeanQy->Write();
    Qx_vs_c->Write();
    Qy_vs_c->Write();
    Qx_nr_histo->Write();
    Qx_r_histo->Write();
    Qy_nr_histo->Write();
    Qy_r_histo->Write();
}

void reader::v_histos()
{
    this->GetQMean();
    this->GetResPsi();
    auto ResSubEv_nr = new TProfile("Res of SE, nr","Resolution of Subevent, not recentred",20,0,100);
    Int_t n_ev = fChain->GetEntries();
    Float_t* PsiEP_nr = new Float_t[2];
    Float_t  cent;

    for(int i=0; i<n_ev; i++)
    {
        fChain->GetEntry(i);
        if(ev->GetPSDEnergy() == 0)
            continue;
        cent=ev->GetCentrality();
        if (cent < 0)
            continue;
        PsiEP_nr = this->GetPsiEP_SE(0);
        Float_t res_nr = cos(PsiEP_nr[0]-PsiEP_nr[1]);
        ResSubEv_nr->Fill(cent,res_nr);
    }
    TFile* f = new TFile("histos_se.root","recreate");
    f->cd();
    ResSubEv_nr->Write();
    ResPsi->Write();
    f->Close();
}

Float_t reader::Get_v( Float_t PsiEP, Float_t phi, Int_t n=1  )
{
    Float_t cent = ev->GetCentrality();
    Int_t bin = (cent/5)+1;
    //cout << bin << endl;
    Float_t res = ResPsi->GetBinContent(bin);
    Float_t v_n=0;
    v_n = cos( n*(phi-PsiEP) ) / sqrt(res);
    if(v_n != v_n)
        cout << "PsiEP=" << PsiEP << " phi=" << phi << " res=" << res << " v1=" << v_n << "centrality =" << cent <<endl;
    return v_n;
}

void reader::GetFlows()
{
    this->GetQMean();
    this->GetResPsi();
    cout <<"building pt&y vs directed flow"<<endl;

    auto y_vs_v1_cent = new TProfile("y vs v1, 0-20","rapidity vs directed flow, centrality 0-20\%, protons",15,-1,2);
    auto pt_vs_v1_cent = new TProfile("pt vs v1, 0-20","transverse momentum vs directed flow, centrality 0-20\%, protons",10,0.2,2.0);
    
    auto y_vs_v1_mid = new TProfile("y vs v1, 20-30","rapidity vs directed flow, centrality 20-30\%, protons",15,-1,2);
    auto pt_vs_v1_mid = new TProfile("pt vs v1, 20-30","transverse momentum vs directed flow, centrality 20-30\%, protons",10,0.2,2.0);

    auto y_vs_v1_per = new TProfile("y vs v1, 30-50","rapidity vs directed flow, centrality 30-50\%, protons",15,-1,2);
    auto pt_vs_v1_per = new TProfile("pt vs v1, 30-50","transverse momentum vs directed flow, centrality 30-50\%, protons",10,0.2,2.0);

    Int_t n_ev = fChain->GetEntries();
    DataTreeTrack* tr;
    Float_t PsiEP;
    Float_t pt, y, w, phi;
    Float_t cent;
    Float_t v1;
    Float_t y0 = 0.5*log(1.23*197+156.743) - 0.5*log(1.23*197-156.743);
    for(int i=0;i<n_ev;i++)
    {
        fChain->GetEntry(i);
        if( ev->GetPSDEnergy()==0 )
            continue;
        Int_t n_tr = ev->GetNVertexTracks();
        PsiEP = this->GetPsiEP();
        if (n_tr < 5)
            continue;
        cent = ev->GetCentrality();
        if(cent < 0)
        {
            continue;
        }
        for(int j=0;j<n_tr;j++)
        {
            tr = ev->GetVertexTrack(j);
            if(tr->GetPdgId()!=14)
                continue;
            y = tr->GetRapidity();
            pt= tr->GetPt();
            y-=y0;
            phi = tr->GetPhi();
            v1 = this->Get_v(PsiEP,phi);
            if(y>=0)
                w=1;
            if(y<0)
                w=-1;
            if (cent < 20) 
            {
                y_vs_v1_cent->Fill(y,v1);
                pt_vs_v1_cent->Fill(pt,w*v1);
            }
            if( (cent>=20)&&(cent<30) )
            {
                y_vs_v1_mid->Fill(y,v1);
                pt_vs_v1_mid->Fill(pt,w*v1);
            }
            if( (cent>=30)&&(cent<50) )
            {
                y_vs_v1_per->Fill(y,v1);
                pt_vs_v1_per->Fill(pt,w*v1);
            }
        }
        this->ShowProgress(i,n_ev);
    }
    TFile* f = new TFile("histos_flow.root","recreate");
    f->cd();
    y_vs_v1_cent->Write();
    pt_vs_v1_cent->Write();
    y_vs_v1_mid->Write();
    pt_vs_v1_mid->Write();
    y_vs_v1_per->Write();
    pt_vs_v1_per->Write();
    f->Close();
}

void reader::GetQMean() 
{
    cout << "Filling <Q> as a function of centrality" << endl;
    Float_t* Q = new Float_t[2];
    Int_t N_ev = (Int_t) fChain->GetEntries();
    Float_t cent;
    for(int i=0;i<N_ev;i++)
    {
        fChain->GetEntry(i);
        if( ev->GetPSDEnergy() == 0 )
            continue;
        cent = ev->GetCentrality();
        if( cent < 0)
            continue;
        Q = this->GetQ();
        MeanQx->Fill(cent, Q[0]);
        MeanQy->Fill(cent, Q[1]);
        this->ShowProgress(i,N_ev);
    }
}

void reader::GetResPsi()
{
    cout << "Filling subevent resolution function of centrality" << endl;
    Int_t n_ev = fChain->GetEntries();
    Float_t* PsiEP = new Float_t[2];
    Float_t  cent;
    for(int i=0; i<n_ev; i++)
    {
        fChain->GetEntry(i);
        if(ev->GetPSDEnergy() == 0)
            continue;
        cent = ev->GetCentrality();
        PsiEP = this->GetPsiEP_SE();
        Float_t res = cos(PsiEP[0]-PsiEP[1]);
        ResPsi->Fill(cent,res);
        this->ShowProgress(i,n_ev);
    }
}

Float_t* reader::GetPsiEP_SE(Bool_t recentring=1) // random subevent method
{
    Float_t cent = ev->GetCentrality();
    Float_t meanx = 0, meany=0;
    Int_t bin = (cent/5)+1;;
    Float_t Qx[2], Qy[2];
    meanx = MeanQx->GetBinContent(bin);
    meany = MeanQy->GetBinContent(bin);
    for(int i=0;i<2;i++)
    {
        Qx[i]=0;
        Qy[i]=0;
    }
    Int_t n_psd_modules = ev->GetNPSDModules();
    DataTreePSDModule* psd_module;
    Float_t phi=0,q=0;
    Float_t EofAll = ev->GetPSDEnergy();
    for(int i=0;i<n_psd_modules;i++)
    {
        psd_module = ev->GetPSDModule(i);
        if(psd_module->GetId() < 0)
            continue;
        phi =        psd_module->GetPhi();
        q   =        psd_module->GetEnergy();
        Int_t p = gRandom->Uniform(0,2);
        Qx[p] += q*cos(phi);
        Qy[p] += q*sin(phi);
    }

    Float_t* PsiEP = new Float_t[2];
    
    for(int i=0;i<2;i++)
    {
        Qx[i]/=0.5*EofAll;
        Qy[i]/=0.5*EofAll;
        if(recentring)
        {
            Qx[i]-= meanx; 
            Qy[i]-= meany;
        }
        PsiEP[i] = atan2(Qy[i],Qx[i]);
//        cout << "Qx = " << Qx[i] << " Qy = "<< Qy[i] << " recentring = " << recentring << " PsiEP = " << PsiEP[i] << " E_all = " << EofAll << endl;
    }
    return PsiEP;
}

Float_t reader::GetPsiEP()
{
    Float_t* Q =(Float_t*) this->GetQ();
    Float_t cent = ev->GetCentrality();
    Int_t bin = (cent/5) + 1;
    Float_t meanx=MeanQx->GetBinContent(bin);
    Float_t meany=MeanQy->GetBinContent(bin);
    Q[0]-=meanx; Q[1]-=meany;
//    cout << "Qx=" << Q[0] << " Qy=" << Q[1] << endl;
    Float_t PsiEP=atan2(Q[1],Q[0]);
    return PsiEP;
}

Bool_t  reader::evSelector()
{
    if (ev->GetNVertexTracks() < 5)
        return 0;
    if (ev->GetNTOFHits() < 5)
        return 0;
    if (ev->GetPSDEnergy() == 0)
        return 0;
    return 1;
}

Float_t* reader::GetQ()
{
    Int_t n_psd_modules;
    DataTreePSDModule* psd_module;
    Float_t q, phi;
    Float_t Qx=0, Qy=0;
    Float_t EofAll=0;
    Float_t* Q = new Float_t[2];
    EofAll = ev->GetPSDEnergy();
    if (EofAll == 0)
    {
        Q[0]=1; Q[1] =1;
        return Q;
    }
    n_psd_modules = ev->GetNPSDModules();
    for(int j=0;j<n_psd_modules;j++)
    {
        psd_module = ev->GetPSDModule(j);
        if(psd_module->GetId() < 0)
            continue;
        phi =        psd_module->GetPhi();
        q   =        psd_module->GetEnergy();
        Qx+= q*cos(phi); 
        Qy+= q*sin(phi);
    }
    Qx/=EofAll; Qy/=EofAll;
    Q[0]=Qx; Q[1]=Qy;
    return Q;
}

void reader::ShowProgress(Int_t i, Int_t N)
{
    Float_t barsize=50;
    Int_t progress = barsize*i/N-1;
    if(i == 0)
    {
        cout << "[" << flush;
        for(int j=0;j<barsize-2;j++)
            cout << "." << flush;
        cout << "]" << flush;
    }
    if( (i>0)&&(i<N-2) )
    {
        for(int j=0;j<barsize;j++)
            cout<<"\r"<<flush;
        cout << "[" <<flush;
        for(int j=0;j<=progress;j++)
            cout<<"#"<<flush;
        for(int j=0;j<barsize-progress-2;j++)
            cout<<"."<<flush;
        cout<<"]"<<flush;
    }
    if( i>N-2 )
    {
        for(int j=0;j<barsize;j++)
            cout<<"\r"<<flush;
        cout << "[" << flush;
        for(int j=0;j<barsize-2;j++)
            cout<<"#"<<flush;
        cout<<"]"<<endl;
    }
}