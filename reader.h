class reader 
{
    private:
    TChain* fChain;
    TBranch* DTEvent;
    public:
    DataTreeEvent* ev;
    Float_t MeanQx, MeanQy;
    reader();
    ~reader();
    void    GetEvent(Int_t n);
    Bool_t  evSelector();
    Bool_t  Check_track( Int_t tr_num );
    Float_t Get_v(Float_t* PsiEP, Float_t phi, Int_t n=1);
    Float_t* GetPsiEP( );
    void    Main_Loop();
    void    GetQMean(); 
};