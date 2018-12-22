#include "reader.h"
#include "HADES_constants.h"

reader::reader()
{
    fChain = new TChain("DataTree");
    //fChain->Add("~/pool/1/mam2mih/out.root");
    fChain->Add("AuAu.root");
    ev = new DataTreeEvent;
    DTEvent = (TBranch*) fChain->GetBranch("DTEvent");
    DTEvent->SetAddress(&ev);

    MeanQx = new TProfile("c vs Qx","centrality vs Qx",10,0,50);
    MeanQx->GetYaxis()->SetTitle("<Qx>");
    MeanQx->GetXaxis()->SetTitle("centrality");

    MeanQy = new TProfile("c vs Qy","centrality vs Qy",10,0,50);
    MeanQy->GetYaxis()->SetTitle("<Qy>");
    MeanQy->GetXaxis()->SetTitle("centrality");

    ResPsi = new TProfile("resolution of SE, r","Resolution of subevent method, recentred",10,0,50);

    this->GetQMean();
    this->GetResPsi();
    cout << "Initialization was succsess" << endl;
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
    cout << "Q-vector components distributions are building" << endl;
    Float_t* Q = new Float_t[2];
    Int_t n_ev = fChain->GetEntries();
    TProfile* Qx_vs_c = new TProfile("Qx vs c, recentred", "Qx vs centrality after recentring",10,0,50);
    TProfile* Qy_vs_c = new TProfile("Qy vs c, recentred", "Qy vs centrality after recentring",10,0,50);
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
        if ( ! this->evSelector() )
            continue;
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
    TFile* f = new TFile("histos_Q.root","recreate");
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
    auto ResSubEv_nr = new TProfile("Res of SE, nr","Resolution of Subevent, not recentred",10,0,50);
    Int_t n_ev = fChain->GetEntries();
    Float_t* PsiEP_nr = new Float_t[2];
    Float_t  cent;

    for(int i=0; i<n_ev; i++)
    {
        fChain->GetEntry(i);
        if( !this->evSelector() )
            continue;
        if(ev->GetPSDEnergy() == 0)
            continue;
        cent=ev->GetCentrality();
        if (cent < 0)
            continue;
        PsiEP_nr = this->GetPsiEP_SE(0);
        Float_t res_nr = 2*cos(PsiEP_nr[0]-PsiEP_nr[1]);
        ResSubEv_nr->Fill(cent,res_nr);
    }
    TFile* f = new TFile("histos_res.root","recreate");
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
    cout <<"Directed flow as function of rapidity and transverse momentum is building"<< endl;

    auto y_vs_v1_cent = new TProfile("y vs v1, 0-20","rapidity vs directed flow, centrality 0-20\%, protons",10,-1,1);
    auto pt_vs_v1_cent = new TProfile("pt vs v1, 0-20","transverse momentum vs directed flow, centrality 0-20\%, protons",10,0.,2.0);
    
    auto y_vs_v1_mid = new TProfile("y vs v1, 20-30","rapidity vs directed flow, centrality 20-30\%, protons",10,-1,1);
    auto pt_vs_v1_mid = new TProfile("pt vs v1, 20-30","transverse momentum vs directed flow, centrality 20-30\%, protons",10,0.,2.0);

    auto y_vs_v1_per = new TProfile("y vs v1, 30-50","rapidity vs directed flow, centrality 30-50\%, protons",10,-1,1);
    auto pt_vs_v1_per = new TProfile("pt vs v1, 30-50","transverse momentum vs directed flow, centrality 30-50\%, protons",10,0.,2.0);

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
        if ( ! this->evSelector() )
            continue;
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
            if ( ! this->trSelector(j) ) 
                continue;
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
    cout << "Mean Q-vectors componets as functions of centrality are building" << endl;
    Float_t* Q = new Float_t[2];
    Int_t N_ev = (Int_t) fChain->GetEntries();
    Float_t cent;
    for(int i=0;i<N_ev;i++)
    {
        fChain->GetEntry(i);
        if ( !this->evSelector() )
            continue;
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
    cout << "Subevent resolution as a function of centrality is building" << endl;
    Int_t n_ev = fChain->GetEntries();
    Float_t* PsiEP = new Float_t[2];
    Float_t  cent;
    for(int i=0; i<n_ev; i++)
    {
        fChain->GetEntry(i);
        if ( !this->evSelector() )
            continue;
        if(ev->GetPSDEnergy() == 0)
            continue;
        cent = ev->GetCentrality();
        PsiEP = this->GetPsiEP_SE();
        Float_t res = 2*cos(PsiEP[0]-PsiEP[1]);
        ResPsi->Fill(cent,res);
        this->ShowProgress(i,n_ev);
    }
}

Float_t* reader::GetPsiEP_SE(Bool_t recentring=1) // random subevent method
{
    Float_t cent = ev->GetCentrality();
    Float_t meanx = 0, meany=0;
    Int_t bin = (cent/5)+1;
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
    vector<DataTreePSDModule*> psd_v;
    Float_t phi=0,q=0;
    Float_t EofAll = ev->GetPSDEnergy();
    for(int i=0;i<n_psd_modules;i++)
    {
        psd_module = ev->GetPSDModule(i);
        if(psd_module->GetId() < 0)
            continue;
        psd_v.push_back(psd_module);
    }
    random_shuffle(psd_v.begin(),psd_v.end());
    n_psd_modules = psd_v.size();
    for(int i=0;i<n_psd_modules;i++)
    {
        phi =        psd_v[i]->GetPhi();
        q   =        psd_v[i]->GetEnergy();
        Int_t p =i%2;
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
    if (  ev->GetVertexPositionComponent(2) > 0 || ev->GetVertexPositionComponent(2) < -60 )
        return 0;
    Float_t Rx = ev->GetVertexPositionComponent(0), Ry = ev->GetVertexPositionComponent(1);
    if ( sqrt(Rx*Rx+Ry*Ry) > 3 )
        return 0;
    if ( ev->GetVertexQuality() < 0.5 || ev->GetVertexQuality() > 40 )
        return 0;    
    if ( !ev->GetTrigger(HADES_constants::kGoodVertexClust)->GetIsFired() ) 
        return 0;
    if ( ! ev->GetTrigger(HADES_constants::kGoodVertexCand)->GetIsFired() )
        return 0;
    if( !ev->GetTrigger(HADES_constants::kGoodSTART)->GetIsFired() )
        return 0;
    if( !ev->GetTrigger(HADES_constants::kNoPileUpSTART)->GetIsFired() )
        return 0;
    if( !ev->GetTrigger(HADES_constants::kGoodSTARTVETO)->GetIsFired() )
        return 0;
    if( !ev->GetTrigger(HADES_constants::kGoodSTARTMETA)->GetIsFired() )
        return 0;
    if( !ev->GetTrigger(HADES_constants::kNoVETO)->GetIsFired() )
        return 0;
    return 1;
}

Bool_t reader::trSelector(Int_t idx)
{
    DataTreeTrack* tr = ev->GetVertexTrack(idx);
    DataTreeTOFHit* hit = ev->GetTOFHit(idx);
    Float_t tof = hit->GetTime();
    Float_t len = hit->GetPathLength();
    //cout << len/tof/299.792458 << endl;
    if ( tr->GetDCAComponent(1) > 15 )
        return 0;
    if ( len/tof/299.792458 > 1 ) 
        return 0;
    if ( hit->GetPositionComponent(0) < -5 || hit->GetPositionComponent(0) > 5 )
        return 0;
    if ( hit->GetPositionComponent(1) < -5 || hit->GetPositionComponent(1) > 5 )
        return 0;
    if ( tr->GetChi2() > 100 )
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

void reader::MassQA()
{
    cout << "mass QA histos are building" << endl;
    TH2F* m2_vs_p_all = new TH2F("m^{2} vs p, all","squared mass vs momentum, all particles",100,0,3.5,100,0,2);
    TH2F* m2_vs_p_p = new TH2F("m^{2} vs p, p^{+}","squared mass vs momentum, protons",100,0,3.5,100,0,2);
    
    TH2F* beta_vs_p_all = new TH2F("#beta vs p, all","#beta vs momentum, all particles",100,0,3.5,100,0,1);
    TH2F* beta_vs_p_p = new TH2F("#beta vs p, p^{+}","#beta vs momentum, protons",100,0,3.5,100,0,1);
    
    TH2F* dEdx_vs_p_all = new TH2F("#frac{dE}{dx} vs p, all","Ionization loss vs momentum, all particles",100,0,3.5,100,0,20);
    TH2F* dEdx_vs_p_p = new TH2F("#frac{dE}{dx} vs p, p","Ionization loss vs momentum, protons",100,0,3.5,100,0,20);

    TH1F* p_histo = new TH1F("p","momentum distribution",100,0,3.5);
    Int_t n_ev = fChain->GetEntries();
    DataTreeTrack* tr;
    DataTreeTOFHit* hit;
    Float_t tof, len;
    for(int i=0;i<n_ev;i++)
    {
        fChain->GetEntry(i);
        if (!this->evSelector())
            continue;
        int n_tr = ev->GetNVertexTracks();
        for(int j=0;j<n_tr;j++)
        {
            if(!this->trSelector(j))
            {
                continue;
            }
            tr = ev->GetVertexTrack(j);
            hit = ev->GetTOFHit(j);
            TLorentzVector P = tr->GetMomentum();
            p_histo->Fill(P.P());
            int q = tr->GetCharge();
            tof = hit->GetTime(); len = hit->GetPathLength();
            float beta = len/tof/299.792458;
            float m2 = P.P()*P.P()*(1-beta*beta)/(beta*beta);
            float dEdx = tr->GetdEdx(2);
            //cout << dEdx << endl;
            m2_vs_p_all->Fill(P.P()*q,m2);
            dEdx_vs_p_all->Fill(q*P.P(),dEdx);
            beta_vs_p_all->Fill( P.P()*q,beta );
            if ( tr->GetPdgId() == 14 )
            {
                beta_vs_p_p->Fill( P.P()*q,beta );
                m2_vs_p_p->Fill(P.P()*q,m2);
                dEdx_vs_p_p->Fill(q*P.P(),dEdx);
            }
        }
        this->ShowProgress(i,n_ev);
    }
    TFile* f = new TFile("histos_mass.root","recreate");
    f->cd();
    m2_vs_p_all->Write();
    m2_vs_p_p->Write();
    beta_vs_p_all->Write();
    beta_vs_p_p->Write();
    dEdx_vs_p_all->Write();
    dEdx_vs_p_p->Write();
    f->Close();
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
    if( i>N-3 )
    {
        for(int j=0;j<barsize;j++)
            cout<<"\r"<<flush;
        cout << "[" << flush;
        for(int j=0;j<barsize-1;j++)
            cout<<"#"<<flush;
        cout<<"]"<<endl;
    }
}