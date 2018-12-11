class reader 
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    public:
    DataTreeEvent* ev;
    TProfile* MeanQx; 
    TProfile* MeanQy;
    TProfile* ResPsi;
    reader();
    ~reader();
    void    GetEvent(Int_t n);
    Bool_t  evSelector();
    Bool_t  Check_track( Int_t tr_num );
    Float_t Get_v(Float_t* PsiEP, Float_t phi, Int_t n=1);
    Float_t* GetPsiEP(Bool_t recentring=1);
    Float_t* GetQ();
    void    GetResPsi();
    void    Q_histos();
    void    v_histos();
    void    GetQMean(); 
    void    GetFlows();
    void    ShowProgress(Int_t i, Int_t N);
};