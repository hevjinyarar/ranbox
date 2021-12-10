#define tree_cxx
#include "tree.h"
#define _USE_MATH_DEFINES
#include <math.h>
#define SGN(x)  ( ((x) < 0) ?  -1 : ( ((x) == 0 ) ? 0 : 1) )
#define ANG(x,y)  ( M_PI-fabs(fabs(y-x)-M_PI) )
#define VECSUMX(x,phix,y,phiy,z,phiz) (x*cos(phix)+y*cos(phiy)+z*cos(phiz))
#define VECSUMY(x,phix,y,phiy,z,phiz) (x*sin(phix)+y*sin(phiy)+z*sin(phiz))
#define VECSUMZ(x,thetax,y,thetay,z,thetaz) (x/tan(thetax)+y/tan(thetay)+z/tan(thetaz))
#define EN(m1,pt1,m2,pt2) (sqrt(pow(m1,2)+pow(pt1,2))+sqrt(pow(m2,2)+pow(pt2,2)))
#define MASS(en,px,py,pz) (sqrt(pow(en,2)-pow(px,2)-pow(py,2)-pow(pz,2)))
#define NORM(x,y,z) (sqrt(pow(x,2)+pow(y,2)+pow(z,2)))
#define THETA(x) (2*atan(exp(-x)))
#define ETA(x) (-log(tan(x/2)))
#include <TH2.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <vector>

struct par{int ch;float phi;float pt;float mass;float angle;float eta;float theta;float en;float p;};

void tree::Loop(){
   // TFile *out= new TFile("histosRight.root","RECREATE");
    TRandom2 *rand=new TRandom2();
    TH1F *h1_a1g = new TH1F("1_angle1g","Angle between mu1 and same charge pion, set 1, Good",70,-0.2,1.7);
    TH1F *h1_a2g = new TH1F("1_angle2g","Angle between mu1 and opp charge pion, set 1, Good",90,-1.7,1.7);
    TH1F *h1_a1b = new TH1F("1_angle1b","Angle between mu1 and same charge pion, set 1,Bad",50,-0.2,3.5);
    TH1F *h1_a2b = new TH1F("1_angle2b","Angle between mu1 and opp charge pion, set 1, Bad",70,-3.5,3.5);
    
    TH1F *h2_a1g = new TH1F("2_angle1g","Angle between mu1 and same charge pion, set 2, Good",40,-0.2,1.7);
    TH1F *h2_a2g = new TH1F("2_angle2g","Angle between mu1 and opp charge pion, set 2, Good",40,-1.7,1.7);
    TH1F *h2_a1b = new TH1F("2_angle1b","Angle between mu1 and same charge pion, set 2, Bad",40,-0.2,3.5);
    TH1F *h2_a2b = new TH1F("2_angle2b","Angle between mu1 and opp charge pion, set 2, Bad",40,-3.5,3.5);

    TH1F *h3_a1g = new TH1F("3_angle1g","Angle between mu1 and same charge pion, set 3, Good",30,-0.2,1.7);
    TH1F *h3_a2g = new TH1F("3_angle2g","Angle between mu1 and opp charge pion, set 3, Good",30,-1.7,1.7);
    TH1F *h3_a1b = new TH1F("3_angle1b","Angle between mu1 and same charge pion, set 3, Bad",40,-0.2,3.5);
    TH1F *h3_a2b = new TH1F("3_angle2b","Angle between mu1 and opp charge pion, set 3, Bad",40,-3.5,3.5);
    
    TH1F *h4_a1g = new TH1F("4_angle1g","Angle between mu1 and same charge pion, set 4, Good",30,-0.2,1.7);
    TH1F *h4_a2g = new TH1F("4_angle2g","Angle between mu1 and opp charge pion, set 4, Good",30,-1.7,1.7);
    TH1F *h4_a1b = new TH1F("4_angle1b","Angle between mu1 and same charge pion, set 4, Bad",40,-0.2,3.5);
    TH1F *h4_a2b = new TH1F("4_angle2b","Angle between mu1 and opp charge pion, set 4, Bad",40,-3.5,3.5);
    
    TH1F *hrhomass1_1 = new TH1F("hrhomass1_1","Rhomass1, set 1",70,0,4);
    TH1F *hrhomass1_2 = new TH1F("hrhomass1_2","Rhomass2, set 1",70,0,4);
    TH1F *whrhomass1_1 = new TH1F("whrhomass1_1","Rhomass",70,0,4);
    TH1F *hrhomass1_1c = new TH1F("hrhomass1_1c","Rhomass1, set 1",70,0,4);
    TH1F *hrhomass1_2c = new TH1F("hrhomass1_2c","Rhomass2, set 1",70,0,4);
    
    TH1F *hrhomass2_1 = new TH1F("hrhomass3_1","Rhomass1, set 2",25,0,4);
    TH1F *hrhomass2_2 = new TH1F("hrhomass3_2","Rhomass2, set 2",25,0,4);
    TH1F *whrhomass2_1 = new TH1F("whrhomass3_1","Rhomass",25,0,4);
    TH1F *hrhomass2_1c = new TH1F("hrhomass3_1c","Rhomass1,set 2",25,0,4);
    TH1F *hrhomass2_2c = new TH1F("hrhomass3_2c","Rhomass2,set 2",25,0,4);
    
    TH2F *hmet1_lg = new TH2F("hmet1_lg","Angle between mu1 and MET, EN<10 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet1_mg = new TH2F("hmet1_mg","Angle between mu1 and MET, 10<EN<20 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet1_hg = new TH2F("hmet1_hg","Angle between mu1 and MET, 20<EN<30 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet1_vhg = new TH2F("hmet1_vhg","Angle between mu1 and MET, EN>30 GeV",40,-0.2,3.5,40,0,10);
    
    TH2F *hmet2_lg = new TH2F("hmet2_lg","Angle between mu1 and MET, EN<10 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet2_mg = new TH2F("hmet2_mg","Angle between mu1 and MET, 10<EN<20 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet2_hg = new TH2F("hmet2_hg","Angle between mu1 and MET, 20<EN<30 GeV",40,-0.2,3.5,40,0,10);
    TH2F *hmet2_vhg = new TH2F("hmet2_vhg","Angle between mu1 and MET, EN>30 GeV",40,-0.2,3.5,40,0,10);
    
    TH1F *heta1_1g=new TH1F("heta1_1g","Eta mu & same ch pi, set 1, Good",50,-4,4);
    TH1F *heta1_1b=new TH1F("heta1_1b","Eta mu & same ch pi, set 1, Bad",50,-6,6);
    TH1F *heta1_2g=new TH1F("heta1_2g","Eta mu & opp ch pi, set 1, Good",70,-4,4);
    TH1F *heta1_2b=new TH1F("heta1_2b","Eta mu & opp ch pi, set 1, Bad",70,-6,6);
    
    TH1F *heta2_1g=new TH1F("heta2_1g","Eta mu & same ch pi, set 2, Good",40,-4,4);
    TH1F *heta2_1b=new TH1F("heta2_1b","Eta mu & same ch pi, set 2, Bad",40,-6,6);
    TH1F *heta2_2g=new TH1F("heta2_2g","Eta mu & opp ch pi, set 2, Good",50,-4,4);
    TH1F *heta2_2b=new TH1F("heta2_2b","Eta mu & opp ch pi, set 2, Bad",50,-6,6);
    
    TH1F *heta3_1g=new TH1F("heta3_1g","Eta mu & same ch pi, set 3, Good",30,-4,4);
    TH1F *heta3_1b=new TH1F("heta3_1b","Eta mu & same ch pi, set 3, Bad",50,-6,6);
    TH1F *heta3_2g=new TH1F("heta3_2g","Eta mu & opp ch pi, set 3, Good",30,-4,4);
    TH1F *heta3_2b=new TH1F("heta3_2b","Eta mu & opp ch pi, set 3, Bad",50,-6,6);
    
    TH1F *heta4_1g=new TH1F("heta4_1g","Eta mu & same ch pi, set 4, Good",30,-4,4);
    TH1F *heta4_1b=new TH1F("heta4_1b","Eta mu & same ch pi, set 4, Bad",50,-6,6);
    TH1F *heta4_2g=new TH1F("heta4_2g","Eta mu & opp ch pi, set 4, Good",30,-4,4);
    TH1F *heta4_2b=new TH1F("heta4_2b","Eta mu & opp ch pi, set 4, Bad",50,-6,6);
    
    TH1F *right1=new TH1F("right1","Right candidates in Bad subset",2,0,1);
    TH1F *right2=new TH1F("right2","Right candidates in Bad subset",2,0,1);
    TH1F *right3=new TH1F("right3","Right candidates in Bad subset",2,0,1);
    TH1F *right4=new TH1F("right4","Right candidates in Bad subset",2,0,1);
    
    TH1F *cuts1=new TH1F("cuts1","",6,0,6);
    TH1F *cuts2=new TH1F("cuts2","",6,0,6);
    TH1F *cuts3=new TH1F("cuts3","",6,0,6);
    TH1F *cuts4=new TH1F("cuts4","",6,0,6);
    
    TH1F *hrhomass1l = new TH1F("hrhomass1l","Rhomass, Pt<4,  set 1",50,0,2);
    TH1F *hrhomass1m = new TH1F("hrhomass1m","Rhomass, 4<Pt<7, set 1",50,0,2);
    TH1F *hrhomass1h = new TH1F("hrhomass1h","Rhomass, 7<Pt<14, set 1",50,0,2);
    TH1F *hrhomass1vh = new TH1F("hrhomass1vh","Rhomass, Pt>14, set 1",40,0,2);
    TH1F *hwrhomass1 = new TH1F("hwrhomass1","Rhomass",40,0,2);
    
    TH1F *hkmass1l = new TH1F("hkmass1l","Rhomass, Pt<4,  set 1",50,0,2);
    TH1F *hkmass1m = new TH1F("hkmass1m","Rhomass, 4<Pt<7, set 1",50,0,2);
    TH1F *hkmass1h = new TH1F("hkmass1h","Rhomass, 7<Pt<14, set 1",50,0,2);
    TH1F *hkmass1vh = new TH1F("hkmass1vh","Rhomass, Pt>14, set 1",50,0,2);
    TH1F *hwkmass1 = new TH1F("hwkmass1","Rhomass",40,0,2);
    
    TH1F *hkmass2l = new TH1F("hkmass2l","Rhomass, Pt<4,  set 1",25,0,2);
    TH1F *hkmass2m = new TH1F("hkmass2m","Rhomass, 4<Pt<7, set 1",25,0,2);
    TH1F *hkmass2h = new TH1F("hkmass2h","Rhomass, 7<Pt<14, set 1",25,0,2);
    TH1F *hkmass2vh = new TH1F("hkmass2vh","Rhomass, Pt>14, set 1",25,0,2);
    TH1F *hwkmass2 = new TH1F("hwkmass2","Rhomass",25,0,2);
    
    TH1F *hrhomass2l = new TH1F("hrhomass2l","Rhomass, Pt<4, set 2",40,0,2);
    TH1F *hrhomass2m = new TH1F("hrhomass2m","Rhomass, 4<Pt<7, set 2",40,0,2);
    TH1F *hrhomass2h = new TH1F("hrhomass2h","Rhomass, 7<Pt<14, set 2",40,0,2);
    TH1F *hrhomass2vh = new TH1F("hrhomass2vh","Rhomass, Pt>14, set 2",40,0,2);
    TH1F *hwrhomass2 = new TH1F("hwrhomass2","Rhomass",40,0,2);
    
    TH1F *hkkmass1_1 = new TH1F("hkkmass1_1","Rhomass1,set 1",70,0,50);
    TH1F *hkkmass1_2 = new TH1F("hkkmass1_2","Rhomass2,set 1",70,0,4);
    TH1F *hkkmass2_1 = new TH1F("hkkmass2_1","Rhomass2,set 1",70,0,4);
    TH1F *hkkmass2_2 = new TH1F("hkkmass2_2","Rhomass2,set 2",70,0,4);
    
    TH1F *hmetrespi = new TH1F("hmetrespi","Angle between Met and Res Pt, set 1 ",100,-3.5,3.5);
    TH1F *hmetrate=new TH1F("Met Rate","",2,0,1);
    
    TH1F *hkkmasstot1 = new TH1F("hkkmasstot1","Rhomass1,set 1",500,0,4);
    TH1F *hkkmasstot2 = new TH1F("hkkmasstot2","Rhomass2,set 1",500,0,4);
    /*
    TH1F *hmet1_l = new TH1F("anglemet1_l","Angle between mu1 and MET, EN<10 GeV",40,-0.2,3.5);
    TH1F *hmet1_m = new TH1F("anglemet1_m","Angle between mu1 and MET, 10<EN<20 GeV",40,-0.2,3.5);
    TH1F *hmet1_h = new TH1F("anglemet1_h","Angle between mu1 and MET, 20<EN<30 GeV",40,-0.2,3.5);
    TH1F *hmet1_vh = new TH1F("anglemet1_vh","Angle between mu1 and MET, EN>30 GeV",40,-0.2,3.5);
   
    TH1F *hmet3_l = new TH1F("anglemet3_l","Angle between mu1 and MET, EN<10 GeV",20,-0.2,3.5);
    TH1F *hmet3_m = new TH1F("anglemet3_m","Angle between mu1 and MET, 10<EN<20 GeV",20,-0.2,3.5);
    TH1F *hmet3_h = new TH1F("anglemet3_h","Angle between mu1 and MET, 20<EN<30 GeV",20,-0.2,3.5);
    TH1F *hmet3_vh = new TH1F("anglemet3_vh","Angle between mu1 and MET, EN>30 GeV",20,-0.2,3.5);
   
    TH1F *hwhmet_l = new TH1F("wanglemet_l","Angle between mu1 and MET, EN<10 GeV",40,-0.2,3.5);
    TH1F *hwhmet_m = new TH1F("wanglemet_m","Angle between mu1 and MET, 10<EN<20 GeV",40,-0.2,3.5);
    TH1F *hwhmet_h = new TH1F("wanglemet_h","Angle between mu1 and MET, 20<EN<30 GeV",40,-0.2,3.5);
    TH1F *hwhmet_vh = new TH1F("wanglemet_vh","Angle between mu1 and MET, EN>30 GeV",40,-0.2,3.5);
   
    TH1F *whmet_l = new TH1F("hwanglemet_l","Angle between mu1 and MET, EN<10 GeV",40,-0.2,3.5);
    TH1F *whmet_m = new TH1F("hwanglemet_m","Angle between mu1 and MET, 10<EN<20 GeV",40,-0.2,3.5);
    TH1F *whmet_h = new TH1F("hwanglemet_h","Angle between mu1 and MET, 20<EN<30 GeV",40,-0.2,3.5);
    TH1F *whmet_vh = new TH1F("hwanglemet_vh","Angle between mu1 and MET, EN>30 GeV",40,-0.2,3.5);
    
    TH1F *hrig1 = new TH1F("rig1","right 1 set",2,0,1);
    TH1F *hrig2 = new TH1F("rig2","right 2 set",2,0,1);
     */
    //  THF2 *******************************************************************************************

   
    
   
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<< "nentries "<<nentries<<endl;
    par MU;
    par PI1;
    par PI2;
    par PI3;
    par RESPI;
    par RHO1;
    par RHO2;
    par MET;
    float sumet=0;
    int a=0;
    int b=0;
    int c=0;
    int w=0;
    int pw=0;
    int numw=0;
    int numhw=0;
    int numc=0;
    int nset1=0;
    int nset2=0;
    int nset3=0;
    int nset4=0;
    int totc=0;
    int totw=0;
    
    int fang1=0;
    int fang2=0;
    int fang3=0;
    int feta1=0;
    int feta2=0;
    int feta3=0;
    
    int tang1=0;
    int tang2=0;
    int tang3=0;
    int teta1=0;
    int teta2=0;
    int teta3=0;
    
    int oang1=0;
    int oang2=0;
    int oang3=0;
    int oeta1=0;
    int oeta2=0;
    int oeta3=0;
    
    int trang1=0;
    int trang2=0;
    int trang3=0;
    int treta1=0;
    int treta2=0;
    int treta3=0;
    float rhomass1=0;
    float rhomass2=0;
    float wrhomass=0;
    float kmass1=0;
    float kmass2=0;
    float wkmass=0;
    int n=0;
    int r1=0;
    int r2=0;
    int r3=0;
    int r4=0;
    int rigtot=0;
    const float maxphi=1.5;
    const float maxeta=2;
    float ptpix=0;
    float ptpiy=0;
    float piresphi=0;
    float rep2g[2]={0,0};
    float rep2b[2]={0,0};
    float rep4g[3]={0,0,0};
    float rep4b[3]={0,0,0};
    float rann=0;
    int allpisame=0;
    float metrespiphi=0;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        b_BsTauTau_mu1_eta->GetEntry(jentry);
        b_BsTauTau_mu1_pt->GetEntry(jentry);
        b_BsTauTau_mu1_q->GetEntry(jentry);
        b_BsTauTau_mu1_phi->GetEntry(jentry);
        
        b_BsTauTau_tau_pi1_eta->GetEntry(jentry);
        b_BsTauTau_tau_pi1_charge->GetEntry(jentry);
        b_BsTauTau_tau_pi1_phi->GetEntry(jentry);
        b_BsTauTau_tau_pi1_pt->GetEntry(jentry);
        b_BsTauTau_tau_pi1_mass->GetEntry(jentry);
        
        b_BsTauTau_tau_pi2_eta->GetEntry(jentry);
        b_BsTauTau_tau_pi2_charge->GetEntry(jentry);
        b_BsTauTau_tau_pi2_phi->GetEntry(jentry);
        b_BsTauTau_tau_pi2_pt->GetEntry(jentry);
        b_BsTauTau_tau_pi2_mass->GetEntry(jentry);
        
        b_BsTauTau_tau_pi3_eta->GetEntry(jentry);
        b_BsTauTau_tau_pi3_charge->GetEntry(jentry);
        b_BsTauTau_tau_pi3_phi->GetEntry(jentry);
        b_BsTauTau_tau_pi3_pt->GetEntry(jentry);
        b_BsTauTau_tau_pi3_mass->GetEntry(jentry);
        
        b_BsTauTau_B_phi->GetEntry(jentry);
        
        b_BsTauTau_tau_rhomass1->GetEntry(jentry);
        b_BsTauTau_tau_rhomass2->GetEntry(jentry);
        b_BsTauTau_tau_pt->GetEntry(jentry);
        b_BsTauTau_tau_phi->GetEntry(jentry);
        b_BsTauTau_tau_mass->GetEntry(jentry);
        
        b_MET_puppi_phi->GetEntry(jentry);
        b_MET_puppi_et->GetEntry(jentry);
        b_MET_sumEt->GetEntry(jentry);
        
        b_BsTauTau_tau_isRight->GetEntry(jentry);
        b_BsTauTau_tau_isRight1->GetEntry(jentry);
        b_BsTauTau_tau_isRight2->GetEntry(jentry);
        b_BsTauTau_tau_isRight3->GetEntry(jentry);
        c=-1;
        w=-1;
        a=0;
        b=0;
        pw=0;

        rep2g[0]=0;
        rep2g[1]=0;
        rep2b[0]=0;
        rep2b[1]=0;
        
        rep4g[0]=0;
        rep4g[1]=0;
        rep4b[0]=0;
        rep4b[1]=0;
        

        if(BsTauTau_tau_pi1_phi->size()>1){
            for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
            PI1.ch=BsTauTau_tau_pi1_charge->at(s);
            PI2.ch=BsTauTau_tau_pi2_charge->at(s);
            PI3.ch=BsTauTau_tau_pi3_charge->at(s);
            MU.ch=BsTauTau_mu1_q->at(0);
                if((MU.ch+PI1.ch+PI2.ch+PI3.ch)!=0) pw++;
            }
            if( pw==BsTauTau_tau_pi1_phi->size() ) w=1;
            else c=1;
        }
        else{
            PI1.ch=BsTauTau_tau_pi1_charge->at(0);
            PI2.ch=BsTauTau_tau_pi2_charge->at(0);
            PI3.ch=BsTauTau_tau_pi3_charge->at(0);
            MU.ch=BsTauTau_mu1_q->at(0);
            if((MU.ch+PI1.ch+PI2.ch+PI3.ch)==0) c++;
            else if((MU.ch+PI1.ch+PI2.ch+PI3.ch)!=0) w++;
        }
 
        if(c==0){     //  SET 1   *********************************************************************
          for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
              PI1.ch=BsTauTau_tau_pi1_charge->at(s);
              PI2.ch=BsTauTau_tau_pi2_charge->at(s);
              PI3.ch=BsTauTau_tau_pi3_charge->at(s);
              MU.ch=BsTauTau_mu1_q->at(0);
              MU.phi=BsTauTau_mu1_phi->at(0);
              MU.eta=BsTauTau_mu1_eta->at(0);
                if (PI1.ch==MU.ch) {                         // PION 1
                    PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                    PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
                 }
                else if(PI2.ch==MU.ch){                        // PION 2
                    PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                    PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
                }
                else if(PI3.ch==MU.ch){                         // PION 3
                    PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                  
                    PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
              }
              PI1.angle=ANG(MU.phi,PI1.phi);
              PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
              PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
              PI1.theta=THETA(PI1.eta);
              PI2.theta=THETA(PI2.eta);
              PI3.theta=THETA(PI3.eta);
              PI1.p=PI1.pt/sin(PI1.theta);
              PI2.p=PI2.pt/sin(PI2.theta);
              PI3.p=PI3.pt/sin(PI3.theta);
              PI1.eta=MU.eta-PI1.eta;
              PI2.eta=MU.eta-PI2.eta;
              PI3.eta=MU.eta-PI3.eta;
              
                 
              if(PI1.angle<maxphi){ oang1++;
                  if(fabs(PI2.angle)<maxphi){oang2++;
                      if(fabs(PI3.angle)<maxphi){oang3++;
                          if(fabs(PI1.eta)<maxeta){oeta1++;
                              if(fabs(PI2.eta)<maxeta){oeta2++;
                                  if(fabs(PI3.eta)<maxeta){oeta3++;
                                  }}}}}}
            if(PI1.angle<maxphi && fabs(PI2.angle)<maxphi && fabs(PI3.angle)<maxphi && fabs(PI1.eta)<maxeta && fabs(PI2.eta)<maxeta && fabs(PI3.eta)<maxeta ){
                  h1_a1g->Fill(PI1.angle);
                  h1_a2g->Fill(PI2.angle);
                  h1_a2g->Fill(PI3.angle);
                  heta1_1g->Fill(PI1.eta);
                  heta1_2g->Fill(PI2.eta);
                  heta1_2g->Fill(PI2.eta);
                  if(BsTauTau_tau_isRight->at(s)==true){
                      right1->Fill(0.7);
                      r1++;
                  }
                  else right1->Fill(0.3);
                  
                  ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                  ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                  piresphi=atan2(ptpiy,ptpix);
                  RESPI.angle=ANG(MU.phi,piresphi);
                  RESPI.pt=NORM(ptpix,ptpiy,0);
                  rhomass1=MASS(EN(PI1.mass,PI1.p,PI2.mass,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                  rhomass2=MASS(EN(PI1.mass,PI1.p,PI3.mass,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                  wrhomass=MASS(EN(PI3.mass,PI3.p,PI2.mass,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                  hwrhomass1->Fill(wrhomass);

                  kmass1=MASS(EN(0.498,PI1.p,0.498,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                  kmass2=MASS(EN(0.498,PI1.p,0.498,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                  wkmass=MASS(EN(0.498,PI3.p,0.498,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                  hwkmass1->Fill(wkmass);
                 hkkmasstot1->Fill(kmass1);
                 hkkmasstot1->Fill(kmass2);
                  if(RESPI.pt<4){
                      hrhomass1l->Fill(rhomass1);
                      hrhomass1l->Fill(rhomass2);
                      hkmass1l->Fill(kmass1);
                      hkmass1l->Fill(kmass1);
                  }
                      else if(RESPI.pt<7){
                          hrhomass1m->Fill(rhomass1);
                          hrhomass1m->Fill(rhomass2);
                          hkmass1m->Fill(kmass1);
                          hkmass1m->Fill(kmass2);
                      }
                          else if(RESPI.pt<14){
                              hrhomass1h->Fill(rhomass1);
                              hrhomass1h->Fill(rhomass2);
                              hkmass1h->Fill(kmass1);
                              hkmass1h->Fill(kmass2);
                          }
                              else{
                                  hrhomass1vh->Fill(rhomass1);
                                  hrhomass1vh->Fill(rhomass2);
                                  hkmass1vh->Fill(kmass1);
                                  hkmass1vh->Fill(kmass2);
                              }
                 
                  MET.angle=ANG(MU.phi,MET.phi);
                  if(MET.en<10) hmet1_lg->Fill(MET.angle,MET.en/sqrt(sumet));
                  else if(MET.en<20) hmet1_mg->Fill(MET.angle,MET.en/sqrt(sumet));
                  else if(MET.en<30) hmet1_hg->Fill(MET.angle,MET.en/sqrt(sumet));
                  else hmet1_vhg->Fill(MET.angle,MET.en/sqrt(sumet));
                
               }
              
              else{
                  if(PI1.angle>maxphi) cuts1->Fill(0.5);
                  if(fabs(PI2.angle)>maxphi) cuts1->Fill(1.5);
                  if(fabs(PI3.angle)>maxphi) cuts1->Fill(2.5);
                  if(fabs(PI1.eta)>maxeta) cuts1->Fill(3.5);
                  if(fabs(PI2.eta)>maxeta) cuts1->Fill(4.5);
                  if(fabs(PI3.eta)>maxeta) cuts1->Fill(5.5);
                  
              h1_a1b->Fill(PI1.angle);
              h1_a2b->Fill(PI2.angle);
              h1_a2b->Fill(PI3.angle);
              heta1_1b->Fill(PI1.eta);
              heta1_2b->Fill(PI2.eta);
              heta1_2b->Fill(PI2.eta);
              }
//----------------------------------------------------------------------------------------------------
              MET.angle=ANG(MU.phi,MET.phi)*SGN(MET.phi-MU.phi)*SGN(PI1.phi-MU.phi);
              if( -M_PI/4<MET.angle<M_PI/4 || MET.angle<-3*M_PI/4 || MET.angle>3*M_PI/4 ){
                  if( (-M_PI/4<PI1.angle<M_PI/4 || PI1.angle<-3*M_PI/4 || PI1.angle>3*M_PI/4 ) && (-M_PI/4<PI2.angle<M_PI/4 || PI2.angle<-3*M_PI/4 || PI2.angle>3*M_PI/4 ) &&
                     (-M_PI/4<PI3.angle<M_PI/4 || PI3.angle<-3*M_PI/4 || PI3.angle>3*M_PI/4 ) ){
                      ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                      ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                      piresphi=atan2(ptpiy,ptpix);
                      ptpix=VECSUMX(ptpix,piresphi,MU.pt,0,0,0);
                      piresphi=atan2(ptpiy,ptpix);
                      RESPI.angle=ANG(MU.phi,piresphi);
                      RESPI.pt=NORM(ptpix,ptpiy,0);
                      metrespiphi=ANG(RESPI.angle,MET.angle)*SGN(RESPI.angle-MET.angle)*SGN(PI1.phi-MU.phi);
                      hmetrespi->Fill(metrespiphi);
                      if(-M_PI/4<metrespiphi<M_PI/4 || metrespiphi<-3*M_PI/4 || metrespiphi>3*M_PI/4 ){
                          hmetrate->Fill(0.3);
                      }
                      else{
                          hmetrate->Fill(0.6);
                      }
                  }
              }
                  
//----------------------------------------------------------------------------------------------------
              
                  nset1++;
       //   }
         }
        }
        else if(c==1){      //  SET 2   **************************************************************
          for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
              PI1.ch=BsTauTau_tau_pi1_charge->at(s);
              PI2.ch=BsTauTau_tau_pi2_charge->at(s);
              PI3.ch=BsTauTau_tau_pi3_charge->at(s);
              MU.ch=BsTauTau_mu1_q->at(0);
              MU.phi=BsTauTau_mu1_phi->at(0);
              MU.eta=BsTauTau_mu1_eta->at(0);
              if (PI1.ch==MU.ch) {                         // PION 1
                  PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                  PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                  PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                  PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                  PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                  PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                
                  PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                  PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                  PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                  PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                  PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                  PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                  PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                  PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                  PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  MET.en=MET_puppi_et->at(0);
                  MET.phi=MET_puppi_phi->at(0);
                  sumet=MET_sumEt->at(0);
        
                 }
              else if(PI2.ch==MU.ch){                         // PION 2
                  PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                  PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                  PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                  PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                  PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                  PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                 
                  PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                  PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                  PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                  PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                  PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                  PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                  PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                  PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                  PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  MET.en=MET_puppi_et->at(0);
                 MET.phi=MET_puppi_phi->at(0);
                 sumet=MET_sumEt->at(0);
            
              }
              else if(PI3.ch==MU.ch){                         // PION 3
                  PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                  PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                  PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                  PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                  PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                  PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                
                  PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                  PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                  PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                  PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                  PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                  PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                  PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                  PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                  PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                  MET.en=MET_puppi_et->at(0);
                  MET.phi=MET_puppi_phi->at(0);
                  sumet=MET_sumEt->at(0);
                }
              
              PI1.angle=ANG(MU.phi,PI1.phi);
              PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
              PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
              PI1.theta=THETA(PI1.eta);
              PI2.theta=THETA(PI2.eta);
              PI3.theta=THETA(PI3.eta);
              PI1.p=PI1.pt/sin(PI1.theta);
              PI2.p=PI2.pt/sin(PI2.theta);
              PI3.p=PI3.pt/sin(PI3.theta);
              PI1.eta=MU.eta-PI1.eta;
              PI2.eta=MU.eta-PI2.eta;
              PI3.eta=MU.eta-PI3.eta;
     //         cout<<PI1.angle<<" "<<PI2.angle<<" "<<PI3.angle<<" "<<PI1.eta<<" "<<PI2.eta<<" "<<PI3.eta<<" 1"<<endl;
              if(PI1.angle<maxphi && fabs(PI2.angle)<maxphi && fabs(PI3.angle)<maxphi && fabs(PI1.eta)<maxeta && fabs(PI2.eta)<maxeta && fabs(PI3.eta)<maxeta ){
                  if(PI1.pt>rep2g[0]){
                      rep2g[0]=PI1.pt;
                      rep2g[1]=(float)s;
                  }
              }
              else{
                  a++;
                  if(PI1.pt>rep2b[0]){
                      rep2b[0]=PI1.pt;
                      rep2b[1]=(float)s;
                  }
              }
          }
            if(a==BsTauTau_tau_pi1_eta->size()){
                int s=(int)rep2b[1];
                PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                MU.ch=BsTauTau_mu1_q->at(0);
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                if (PI1.ch==MU.ch) {                         // PION 1
                    PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                   
                    PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
          
                   }
                else if(PI2.ch==MU.ch){                         // PION 2
                    PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                 
                    PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
              
                }
                else if(PI3.ch==MU.ch){                         // PION 3
                    PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                    PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                    PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                    PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                    PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                    PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                  
                    PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                    PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                    PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                    PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                    PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                    PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                    PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                    PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                    PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    MET.en=MET_puppi_et->at(0);
                    MET.phi=MET_puppi_phi->at(0);
                    sumet=MET_sumEt->at(0);
                  }
                
                nset2++;
                PI1.angle=ANG(MU.phi,PI1.phi);
                PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
                PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
                PI1.theta=THETA(PI1.eta);
                PI2.theta=THETA(PI2.eta);
                PI3.theta=THETA(PI3.eta);
                PI1.p=PI1.pt/sin(PI1.theta);
                PI2.p=PI2.pt/sin(PI2.theta);
                PI3.p=PI3.pt/sin(PI3.theta);
                PI1.eta=MU.eta-PI1.eta;
                PI2.eta=MU.eta-PI2.eta;
                PI3.eta=MU.eta-PI3.eta;
                
            h2_a1b->Fill(PI1.angle);
            h2_a2b->Fill(PI2.angle);
            h2_a2b->Fill(PI3.angle);
            heta2_1b->Fill(PI1.eta);
            heta2_2b->Fill(PI2.eta);
            heta2_2b->Fill(PI3.eta);
                if(PI1.angle>maxphi) cuts2->Fill(0.5);
                if(fabs(PI2.angle)>maxphi) cuts2->Fill(1.5);
                if(fabs(PI3.angle)>maxphi) cuts2->Fill(2.5);
                if(fabs(PI1.eta)>maxeta) cuts2->Fill(3.5);
                if(fabs(PI2.eta)>maxeta) cuts2->Fill(4.5);
                if(fabs(PI3.eta)>maxeta) cuts2->Fill(5.5);
                    
            }
            else{
          //  else if(BsTauTau_tau_isRight->at((int)rep2g[1])==true){
            int s=(int)rep2g[1];
            PI1.ch=BsTauTau_tau_pi1_charge->at(s);
            PI2.ch=BsTauTau_tau_pi2_charge->at(s);
            PI3.ch=BsTauTau_tau_pi3_charge->at(s);
            MU.ch=BsTauTau_mu1_q->at(0);
            MU.phi=BsTauTau_mu1_phi->at(0);
            MU.eta=BsTauTau_mu1_eta->at(0);
            if (PI1.ch==MU.ch) {                         // PION 1
                PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                PI3.phi=BsTauTau_tau_pi3_phi->at(s);
              
                PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                MET.en=MET_puppi_et->at(0);
                MET.phi=MET_puppi_phi->at(0);
                sumet=MET_sumEt->at(0);
      
               }
            else if(PI2.ch==MU.ch){                         // PION 2
                PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                PI3.phi=BsTauTau_tau_pi3_phi->at(s);
               
                PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                MET.en=MET_puppi_et->at(0);
                MET.phi=MET_puppi_phi->at(0);
                sumet=MET_sumEt->at(0);
          
            }
            else if(PI3.ch==MU.ch){                         // PION 3
                PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                PI3.phi=BsTauTau_tau_pi1_phi->at(s);
               
                PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                MET.en=MET_puppi_et->at(0);
                MET.phi=MET_puppi_phi->at(0);
                sumet=MET_sumEt->at(0);
              }
            nset2++;
                PI1.angle=ANG(MU.phi,PI1.phi);
                PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
                PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
                PI1.theta=THETA(PI1.eta);
                PI2.theta=THETA(PI2.eta);
                PI3.theta=THETA(PI3.eta);
                PI1.p=PI1.pt/sin(PI1.theta);
                PI2.p=PI2.pt/sin(PI2.theta);
                PI3.p=PI3.pt/sin(PI3.theta);
                PI1.eta=MU.eta-PI1.eta;
                PI2.eta=MU.eta-PI2.eta;
                PI3.eta=MU.eta-PI3.eta;
            
            if(PI1.angle<maxphi){ tang1++;
                if(fabs(PI2.angle)<maxphi){tang2++;
                    if(fabs(PI3.angle)<maxphi){tang3++;
                        if(fabs(PI1.eta)<maxeta){teta1++;
                            if(fabs(PI2.eta)<maxeta){teta2++;
                                if(fabs(PI3.eta)<maxeta){teta3++;
                                }}}}}}
                
                ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                piresphi=atan2(ptpiy,ptpix);
                RESPI.angle=ANG(MU.phi,piresphi);
                RESPI.pt=NORM(ptpix,ptpiy,0);
                rhomass1=MASS(EN(PI1.mass,PI1.p,PI2.mass,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                rhomass2=MASS(EN(PI1.mass,PI1.p,PI3.mass,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                wrhomass=MASS(EN(PI3.mass,PI3.p,PI2.mass,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                hwrhomass2->Fill(wrhomass);

                kmass1=MASS(EN(0.498,PI1.p,0.498,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                kmass2=MASS(EN(0.498,PI1.p,0.498,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                wkmass=MASS(EN(0.498,PI3.p,0.498,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                hwkmass2->Fill(wkmass);
                hkkmasstot2->Fill(kmass1);
                hkkmasstot2->Fill(kmass2);
                
                if(RESPI.pt<4){
                    hrhomass2l->Fill(rhomass1);
                    hrhomass2l->Fill(rhomass2);
                    hkmass2l->Fill(kmass1);
                    hkmass2l->Fill(kmass1);
                }
                    else if(RESPI.pt<7){
                        hrhomass2m->Fill(rhomass1);
                        hrhomass2m->Fill(rhomass2);
                        hkmass2m->Fill(kmass1);
                        hkmass2m->Fill(kmass2);
                    }
                        else if(RESPI.pt<14){
                            hrhomass2h->Fill(rhomass1);
                            hrhomass2h->Fill(rhomass2);
                            hkmass2h->Fill(kmass1);
                            hkmass2h->Fill(kmass2);
                        }
                            else{
                                hrhomass2vh->Fill(rhomass1);
                                hrhomass2vh->Fill(rhomass2);
                                hkmass2vh->Fill(kmass1);
                                hkmass2vh->Fill(kmass2);
                            }
            MET.angle=ANG(MU.phi,MET.phi);
            if(MET.en<10) hmet2_lg->Fill(MET.angle,MET.en/sqrt(sumet));
            else if(MET.en<20) hmet2_mg->Fill(MET.angle,MET.en/sqrt(sumet));
            else if(MET.en<30) hmet2_hg->Fill(MET.angle,MET.en/sqrt(sumet));
            else hmet2_vhg->Fill(MET.angle,MET.en/sqrt(sumet));
                
            h2_a1g->Fill(PI1.angle);
            h2_a2g->Fill(PI2.angle);
            h2_a2g->Fill(PI3.angle);
            heta2_1g->Fill(PI1.eta);
            heta2_2g->Fill(PI2.eta);
            heta2_2g->Fill(PI3.eta);
          }
        }
        
        if(w==0){              //  SET 3  *************************************************************
          for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
              PI1.ch=BsTauTau_tau_pi1_charge->at(s);
              PI2.ch=BsTauTau_tau_pi2_charge->at(s);
              PI3.ch=BsTauTau_tau_pi3_charge->at(s);
              MU.ch=BsTauTau_mu1_q->at(0);
              MU.phi=BsTauTau_mu1_phi->at(0);
              MU.eta=BsTauTau_mu1_eta->at(0);
              if (MU.ch==PI1.ch==PI2.ch){
                  float ran=rand->Rndm();
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  }
              }
              else if(MU.ch==PI1.ch==PI3.ch){
                  float ran=rand->Rndm();
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                  }
              }
              else if(MU.ch==PI2.ch==PI3.ch){
                  float ran=rand->Rndm();
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                  }
              }
              PI1.angle=ANG(MU.phi,PI1.phi);
              PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
              PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
              PI1.theta=THETA(PI1.eta);
              PI2.theta=THETA(PI2.eta);
              PI3.theta=THETA(PI3.eta);
              PI1.p=PI1.pt/sin(PI1.theta);
              PI2.p=PI2.pt/sin(PI2.theta);
              PI3.p=PI3.pt/sin(PI3.theta);
              PI1.eta=MU.eta-PI1.eta;
              PI2.eta=MU.eta-PI2.eta;
              PI3.eta=MU.eta-PI3.eta;
              
              if(PI1.angle<maxphi){ trang1++;
                  if(fabs(PI2.angle)<maxphi){trang2++;
                      if(fabs(PI3.angle)<maxphi){trang3++;
                          if(fabs(PI1.eta)<maxeta){treta1++;
                              if(fabs(PI2.eta)<maxeta){treta2++;
                                  if(fabs(PI3.eta)<maxeta){treta3++;
                                  }}}}}}
              if(PI1.angle<maxphi && fabs(PI2.angle)<maxphi && fabs(PI3.angle)<maxphi && fabs(PI1.eta)<maxeta && fabs(PI2.eta)<maxeta && fabs(PI3.eta)<maxeta ){
            
                      h3_a1g->Fill(PI1.angle);
                      h3_a2g->Fill(PI2.angle);
                      h3_a2g->Fill(PI3.angle);
                      heta3_1g->Fill(PI1.eta);
                      heta3_2g->Fill(PI2.eta);
                      heta3_2g->Fill(PI2.eta);
                  }
              else{
                  if(PI1.angle>maxphi) cuts3->Fill(0.5);
                  if(fabs(PI2.angle)>maxphi) cuts3->Fill(1.5);
                  if(fabs(PI3.angle)>maxphi) cuts3->Fill(2.5);
                  if(fabs(PI1.eta)>maxeta) cuts3->Fill(3.5);
                  if(fabs(PI2.eta)>maxeta) cuts3->Fill(4.5);
                  if(fabs(PI3.eta)>maxeta) cuts3->Fill(5.5);

                  h3_a1b->Fill(PI1.angle);
                  h3_a2b->Fill(PI2.angle);
                  h3_a2b->Fill(PI3.angle);
                  heta3_1b->Fill(PI1.eta);
                  heta3_2b->Fill(PI2.eta);
                  heta3_2b->Fill(PI2.eta);
                  if(BsTauTau_tau_isRight->at(s)==true){
                      right3->Fill(0.7);
                      r3++;
                  }
                      else right3->Fill(0.3);
                  }
              }
            nset3++;
        }
        
 //  SET 4 *************************************************************************************************
        
        else if(w==1){
          for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
              PI1.ch=BsTauTau_tau_pi1_charge->at(s);
              PI2.ch=BsTauTau_tau_pi2_charge->at(s);
              PI3.ch=BsTauTau_tau_pi3_charge->at(s);
              MU.ch=BsTauTau_mu1_q->at(0);
              MU.phi=BsTauTau_mu1_phi->at(0);
              MU.eta=BsTauTau_mu1_eta->at(0);
              if (MU.ch==PI1.ch==PI2.ch){
                  float ran=rand->Rndm();
                  rann=ran;
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                  }
              }
              else if(MU.ch==PI1.ch==PI3.ch){
                  float ran=rand->Rndm();
                  rann=ran;
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                  }
              }
              else if(MU.ch==PI2.ch==PI3.ch){
                  float ran=rand->Rndm();
                  rann=ran;
                  if(ran>0.5){
                      PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                  }
                  else{
                      PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                      PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                      PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                      PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                      PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                      PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                      PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                      PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                      PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                      PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                      PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                      PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                      PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                      PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                      PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                  }
              }
              else if(PI2.ch==PI3.ch==PI1.ch){
                  break;
              }
              PI1.angle=ANG(MU.phi,PI1.phi);
              PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
              PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
              PI1.theta=THETA(PI1.eta);
              PI2.theta=THETA(PI2.eta);
              PI3.theta=THETA(PI3.eta);
              PI1.p=PI1.pt/sin(PI1.theta);
              PI2.p=PI2.pt/sin(PI2.theta);
              PI3.p=PI3.pt/sin(PI3.theta);
              PI1.eta=MU.eta-PI1.eta;
              PI2.eta=MU.eta-PI2.eta;
              PI3.eta=MU.eta-PI3.eta;
              
              if(PI1.angle<maxphi && fabs(PI2.angle)<maxphi && fabs(PI3.angle)<maxphi && fabs(PI1.eta)<maxeta && fabs(PI2.eta)<maxeta && fabs(PI3.eta)<maxeta ){
              if(PI1.pt>rep4g[0]){
                  rep4g[0]=PI1.pt;
                  rep4g[1]=(float)s;
                  rep4g[2]=rann;
                 }
              }
              else{
                b++;
                 if(PI1.pt>rep4b[0]){
                  rep4b[0]=PI1.pt;
                  rep4b[1]=(float)s;
                  rep4b[2]=rann;
                 }
              }
          }
            if(PI2.ch==PI3.ch==PI1.ch){
                allpisame++;
                continue;
            }
            if(b==BsTauTau_tau_pi1_eta->size()){
                int s=(int)rep4b[1];
                PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                MU.ch=BsTauTau_mu1_q->at(0);
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                if (MU.ch==PI1.ch==PI2.ch){
                    
                    if(rep4b[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    }
                }
                else if(MU.ch==PI1.ch==PI3.ch){
                    
                    if(rep4b[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                    }
                }
                else if(MU.ch==PI2.ch==PI3.ch){
                    
                    if(rep4b[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    }
                }
                PI1.angle=ANG(MU.phi,PI1.phi);
                PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
                PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
                PI1.theta=THETA(PI1.eta);
                PI2.theta=THETA(PI2.eta);
                PI3.theta=THETA(PI3.eta);
                PI1.p=PI1.pt/sin(PI1.theta);
                PI2.p=PI2.pt/sin(PI2.theta);
                PI3.p=PI3.pt/sin(PI3.theta);
                PI1.eta=MU.eta-PI1.eta;
                PI2.eta=MU.eta-PI2.eta;
                PI3.eta=MU.eta-PI3.eta;
                
          //      if(PI1.angle<maxphi && fabs(PI2.angle)<maxphi && fabs(PI3.angle)<maxphi && fabs(PI1.eta)<maxeta && fabs(PI2.eta)<maxeta && fabs(PI3.eta)<maxeta ){
         //           cout<<PI1.angle<<" "<<PI2.angle<<" "<<PI3.angle<<" "<<PI1.eta<<" "<<PI2.eta<<" "<<PI3.eta<<" BAD  "<<endl;
          //      }
                
                h4_a1b->Fill(PI1.angle);
                h4_a2b->Fill(PI2.angle);
                h4_a2b->Fill(PI3.angle);
                heta4_1b->Fill(PI1.eta);
                heta4_2b->Fill(PI2.eta);
                heta4_2b->Fill(PI3.eta);
                    if(PI1.angle>maxphi) cuts4->Fill(0.5);
                    if(fabs(PI2.angle)>maxphi) cuts4->Fill(1.5);
                    if(fabs(PI3.angle)>maxphi) cuts4->Fill(2.5);
                    if(fabs(PI1.eta)>maxeta) cuts4->Fill(3.5);
                    if(fabs(PI2.eta)>maxeta) cuts4->Fill(4.5);
                    if(fabs(PI3.eta)>maxeta) cuts4->Fill(5.5);
                nset4++;
                        
            }
            else{
              int s=(int)rep4g[1];
                PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                MU.ch=BsTauTau_mu1_q->at(0);
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                if (MU.ch==PI1.ch==PI2.ch){
                    
                    if(rep4g[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                    }
                }
                else if(MU.ch==PI1.ch==PI3.ch){
                    
                    if(rep4g[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi2_eta->at(s);
                    }
                }
                else if(MU.ch==PI2.ch==PI3.ch){
                    
                    if(rep4g[2]>0.5){
                        PI1.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    }
                    else{
                        PI1.ch=BsTauTau_tau_pi3_charge->at(s);
                        PI2.ch=BsTauTau_tau_pi2_charge->at(s);
                        PI3.ch=BsTauTau_tau_pi1_charge->at(s);
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                        PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                        PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                        PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                        PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                        PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                        PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                        PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                        PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                    }
                }
                PI1.angle=ANG(MU.phi,PI1.phi);
                PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
                PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
                PI1.theta=THETA(PI1.eta);
                PI2.theta=THETA(PI2.eta);
                PI3.theta=THETA(PI3.eta);
                PI1.p=PI1.pt/sin(PI1.theta);
                PI2.p=PI2.pt/sin(PI2.theta);
                PI3.p=PI3.pt/sin(PI3.theta);
                PI1.eta=MU.eta-PI1.eta;
                PI2.eta=MU.eta-PI2.eta;
                PI3.eta=MU.eta-PI3.eta;
                
                if(PI1.angle<maxphi){ fang1++;
                    if(fabs(PI2.angle)<maxphi){fang2++;
                        if(fabs(PI3.angle)<maxphi){fang3++;
                            if(fabs(PI1.eta)<maxeta){feta1++;
                                if(fabs(PI2.eta)<maxeta){feta2++;
                                    if(fabs(PI3.eta)<maxeta){feta3++;
                                    }}}}}}
     
             
                  //   else if(set2_g[14]==1){
                      h4_a1g->Fill(PI1.angle);
                      h4_a2g->Fill(PI2.angle);
                      h4_a2g->Fill(PI3.angle);
                      heta4_1g->Fill(PI1.eta);
                      heta4_2g->Fill(PI2.eta);
                      heta4_2g->Fill(PI3.eta);
                nset4++;
            }
         }
      }
   
    // CHECKING NUMBERS OF EVENT AND RIGHT CANDIDATES  *****************************************************
    
    int uc=0;
    int mc=0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        b_BsTauTau_tau_pi1_phi->GetEntry(jentry);
        for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
            if(BsTauTau_tau_pi1_phi->size()==1) uc++;
            if(BsTauTau_tau_pi1_phi->size()>1){
            mc++;
            break;
            }
        }
    }
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        b_BsTauTau_tau_isRight->GetEntry(jentry);
        for (int s=0; s<BsTauTau_tau_isRight->size(); s++){
            if(BsTauTau_tau_isRight->at(s)==true) rigtot++;
        }
    }
    cout<<endl;
    cout<<"set 1: "<<nset1<<endl;
    cout<<"set 2: "<<nset2<<endl;
    cout<<"set 3: "<<nset3<<endl;
    cout<<"set 4: "<<nset4<<endl;
    cout<<"Event with 3 pions with same charge "<<allpisame<<endl;
    cout<<"total event sum = "<<nset1+nset2+nset3+nset4<<endl;
  //  cout<<nset1+nset2<<" = "<<totc<<endl;
  //  cout<<nset3+nset4<<" = "<<totw<<endl;
    cout<<nset1+nset3<<" = "<<uc<<endl;
    cout<<nset2+nset4<<" = "<<mc<<endl;
    
    cout<<"% right in 1 set on tot right, good "<<(float)r1/rigtot*100<<"% on set 1 "<<(float)r1/nset1*100<<endl;
    cout<<"% right in 2 set, bad "<<(float)r2/rigtot*100<<endl;
    cout<<"% right in 3 set, bad "<<(float)r3/rigtot*100<<endl;
    cout<<"% right in 4 set, bad "<<(float)r4/rigtot*100<<endl;
               
    cout<<"remaning after each cut, set1 "<<(float)oang1/nset1*100<<" "<<(float)oang2/nset1*100<<" "<<(float)oang3/nset1*100<<" "<<(float)oeta1/nset1*100<<" "<<(float)oeta2/nset1*100<<" "<<(float)oeta3/nset1*100<<endl;
               
    cout<<"remaning after each cut, set2 "<<(float)tang1/nset2*100<<" "<<(float)tang2/nset2*100<<" "<<(float)tang3/nset2*100<<" "<<(float)teta1/nset2*100<<" "<<(float)teta2/nset2*100<<" "<<(float)teta3/nset2*100<<endl;
               
    cout<<"remaning after each cut, set3 "<<(float)trang1/nset3*100<<" "<<(float)trang2/nset3*100<<" "<<(float)trang3/nset3*100<<" "<<(float)treta1/nset3*100<<" "<<(float)treta2/nset3*100<<" "<<(float)treta3/nset3*100<<endl;
               
    cout<<"remaning after each cut, set4 "<<(float)fang1/nset4*100<<" "<<(float)fang2/nset4*100<<" "<<(float)fang3/nset4*100<<" "<<(float)feta1/nset4*100<<" "<<(float)feta2/nset4*100<<" "<<(float)feta3/nset4*100<<endl;
    
    cout<<"average failed criteria , set 1: "<<(float)(cuts1->GetEntries())/heta1_1b->GetEntries()<<endl;
    cout<<"average failed criteria , set 2: "<<(float)(cuts2->GetEntries())/heta2_1b->GetEntries()<<endl;
    cout<<"average failed criteria , set 3: "<<(float)(cuts3->GetEntries())/heta3_1b->GetEntries()<<endl;
    cout<<"average failed criteria , set 4: "<<(float)(cuts4->GetEntries())/heta4_1b->GetEntries()<<endl;
    
  // out->Write();
  // out->Close();
    cout<<"Done!"<<endl;
    cout<<"List of variables: eta | angle | rho | cuts | right | met "<<endl;
    cout<<"Insert variable name: "<<endl;
    
   // Draw results
    TString var;
    while(1){
        cin>>var;
        if (var=="q") break;
        else if(var=="eta"){
            TCanvas *ce1g= new TCanvas("ce1g","Eta1, All sets, Good",500,100,400,300);
            ce1g->cd();
            TPad *peg1= new TPad("peg1","eta1, set 1, good",0,0.50,0.50,1);
            peg1->Draw();
            peg1->cd();
            heta1_1g->Draw();
            peg1->Update();
            ce1g->cd();
            TPad *peg2= new TPad("peg2","eta1, set 2, good",0.50,0.50,1,1);
            peg2->Draw();
            peg2->cd();
            heta2_1g->Draw();
            peg2->Update();
            ce1g->cd();
            TPad *peg3= new TPad("peg3","eta1, set 3, good",0,0,0.50,0.50);
            peg3->Draw();
            peg3->cd();
            heta3_1g->Draw();
            peg3->Update();
            ce1g->cd();
            TPad *peg4= new TPad("peg4","eta1, set 4, good",0.50,0,1,0.5);
            peg4->Draw();
            peg4->cd();
            heta4_1g->Draw();
            peg4->Update();
            
            TCanvas *ce2g= new TCanvas("ce2g","Eta2, All sets, Good",500,100,400,300);
            ce2g->cd();
            TPad *peg5= new TPad("peg5","eta2, set 1, good",0,0.50,0.50,1);
            peg5->Draw();
            peg5->cd();
            heta1_2g->Draw();
            peg5->Update();
            ce2g->cd();
            TPad *peg6= new TPad("peg6","eta2, set 2, good",0.50,0.50,1,1);
            peg6->Draw();
            peg6->cd();
            heta2_2g->Draw();
            peg6->Update();
            ce2g->cd();
            TPad *peg7= new TPad("peg7","eta2, set 3, good",0,0,0.50,0.50);
            peg7->Draw();
            peg7->cd();
            heta3_2g->Draw();
            peg7->Update();
            ce2g->cd();
            TPad *peg8= new TPad("peg8","eta2, set 4, good",0.50,0,1,0.5);
            peg8->Draw();
            peg8->cd();
            heta4_2g->Draw();
            peg8->Update();
            
            TCanvas *ce1b= new TCanvas("ce1b","Eta1, All sets, Bad",500,100,400,300);
            ce1b->cd();
            TPad *peb1= new TPad("peb1","eta1, set 1, bad",0,0.50,0.50,1);
            peb1->Draw();
            peb1->cd();
            heta1_1b->SetLineColor(kRed);
            heta1_1b->Draw();
            peb1->Update();
            ce1b->cd();
            TPad *peb2= new TPad("peb2","eta1, set 2, bad",0.50,0.50,1,1);
            peb2->Draw();
            peb2->cd();
            heta2_1b->SetLineColor(kRed);
            heta2_1b->Draw();
            peb2->Update();
            ce1b->cd();
            TPad *peb3= new TPad("peb3","eta1, set 3, bad",0,0,0.50,0.50);
            peb3->Draw();
            peb3->cd();
            heta3_1b->SetLineColor(kRed);
            heta3_1b->Draw();
            peb3->Update();
            ce1b->cd();
            TPad *peb4= new TPad("peb4","eta1, set 4, bad",0.50,0,1,0.5);
            peb4->Draw();
            peb4->cd();
            heta4_1b->SetLineColor(kRed);
            heta4_1b->Draw();
            peb4->Update();
            
            TCanvas *ce2b= new TCanvas("ce2b","Eta2, All sets, Bad",500,100,400,300);
            ce2b->cd();
            TPad *peb5= new TPad("peb5","eta1, set 1, bad",0,0.50,0.50,1);
            peb5->Draw();
            peb5->cd();
            heta1_2b->SetLineColor(kRed);
            heta1_2b->Draw();
            peb5->Update();
            ce2b->cd();
            TPad *peb6= new TPad("peb6","eta1, set 2, bad",0.50,0.50,1,1);
            peb6->Draw();
            peb6->cd();
            heta2_2b->SetLineColor(kRed);
            heta2_2b->Draw();
            peb6->Update();
            ce2b->cd();
            TPad *peb7= new TPad("peb7","eta1, set 3, bad",0,0,0.50,0.50);
            peb7->Draw();
            peb7->cd();
            heta3_2b->SetLineColor(kRed);
            heta3_2b->Draw();
            peb7->Update();
            ce2b->cd();
            TPad *peb8= new TPad("peb8","eta1, set 4, bad",0.50,0,1,0.5);
            peb8->Draw();
            peb8->cd();
            heta4_2b->SetLineColor(kRed);
            heta4_2b->Draw();
            peb8->Update();
        }
        else if(var=="met"){
            TCanvas *cmet1g= new TCanvas("cmet1g","angle between met and mu, set 1, good",500,100,400,300);
            cmet1g->cd();
            TPad *pmetl= new TPad("pmetl","met, set 1, good",0,0.50,0.50,1);
            pmetl->Draw();
            pmetl->cd();
            hmet1_lg->Draw();
            hmet1_lg->SetMarkerStyle(7);
            pmetl->Update();
            cmet1g->cd();
            TPad *pmetm= new TPad("pmetm","met, set 1, good",0.50,0.50,1,1);
            pmetm->Draw();
            pmetm->cd();
            hmet1_mg->Draw();
            hmet1_mg->SetMarkerStyle(7);
            pmetm->Update();
            cmet1g->cd();
            TPad *pmeth= new TPad("pmeth","met, set 1, good",0,0,0.50,0.50);
            pmeth->Draw();
            pmeth->cd();
            hmet1_hg->Draw();
            hmet1_hg->SetMarkerStyle(7);
            pmeth->Update();
            cmet1g->cd();
            TPad *pmetvh= new TPad("pmetvh","met, set 1, good",0.50,0,1,0.5);
            pmetvh->Draw();
            pmetvh->cd();
            hmet1_vhg->Draw();
            hmet1_vhg->SetMarkerStyle(7);
            pmetvh->Update();
            
            TCanvas *cmet2g= new TCanvas("cmet2g","angle between met and mu, set 2, good",500,100,400,300);
            cmet2g->cd();
            TPad *pmetl2= new TPad("pmetl2","met, set 2, good",0,0.50,0.50,1);
            pmetl2->Draw();
            pmetl2->cd();
            hmet2_lg->Draw();
            hmet2_lg->SetMarkerStyle(7);
            pmetl2->Update();
            cmet2g->cd();
            TPad *pmetm2= new TPad("pmetm2","met, set 2, good",0.50,0.50,1,1);
            pmetm2->Draw();
            pmetm2->cd();
            hmet2_mg->Draw();
            hmet2_mg->SetMarkerStyle(7);
            pmetm2->Update();
            cmet2g->cd();
            TPad *pmeth2= new TPad("pmeth2","met, set 2, good",0,0,0.50,0.50);
            pmeth2->Draw();
            pmeth2->cd();
            hmet2_hg->Draw();
            hmet2_hg->SetMarkerStyle(7);
            pmeth2->Update();
            cmet2g->cd();
            TPad *pmetvh2= new TPad("pmetvh2","met, set 2, good",0.50,0,1,0.5);
            pmetvh2->Draw();
            pmetvh2->cd();
            hmet2_vhg->Draw();
            hmet2_vhg->SetMarkerStyle(7);
            pmetvh2->Update();
            
            TCanvas *cmetrespi= new TCanvas("cmetrespi","angle between met and mu, set 2, good",500,100,400,300);
            hmetrespi->Draw();
            cmetrespi->cd();
            
            TCanvas *cmetrate= new TCanvas("cmetrate","Rate Met",500,100,400,300);
            hmetrate->Draw();
            cmetrate->cd();
        }
        else if(var=="rho"){
            /*
            TCanvas *crhog= new TCanvas("crhog","",500,100,400,300);
            crhog->cd();
            TPad *prg= new TPad("prg","rhomass1, set 1, good",0,0.50,0.50,1);
            prg->Draw();
            prg->cd();
            hrhomass1_1->Draw();
            hrhomass1_1c->SetLineColor(kGreen);
            hrhomass1_1c->Draw("same");
            whrhomass1_1->SetLineColor(kRed);
            whrhomass1_1->Draw("same");
         //  hkkmass1_1->Draw("same");
         //   hkkmass1_1->SetLineColor(kBlack);
            prg->Update();
            crhog->cd();
            TPad *prg2= new TPad("prg2","rhomass2, set 1, good",0.50,0.50,1,1);
            prg2->Draw();
            prg2->cd();
            hrhomass1_2->Draw();
            hrhomass1_2c->SetLineColor(kGreen);
            hrhomass1_2c->Draw("same");
            whrhomass1_1->SetLineColor(kRed);
            whrhomass1_1->Draw("same");
         //   hkkmass1_2->Draw("same");
         //   hkkmass1_2->SetLineColor(kBlack);
            prg2->Update();
            crhog->cd();
            TPad *prg3= new TPad("prg3","rhomass1, set 2, good",0,0,0.50,0.50);
            prg3->Draw();
            prg3->cd();
            hrhomass2_1->Draw();
            hrhomass2_1c->SetLineColor(kGreen);
            hrhomass2_1c->Draw("same");
            whrhomass2_1->SetLineColor(kRed);
            whrhomass2_1->Draw("same");
         //   hkkmass2_1->Draw("same");
        //    hkkmass2_1->SetLineColor(kBlack);
            prg3->Update();
            crhog->cd();
            TPad *prg4= new TPad("prg4","rhomass2, set 2, good",0.50,0,1,0.5);
            prg4->Draw();
            prg4->cd();
            hrhomass2_2->Draw();
            hrhomass2_2c->SetLineColor(kGreen);
            hrhomass2_2c->Draw("same");
            whrhomass2_1->SetLineColor(kRed);
            whrhomass2_1->Draw("same");
         //   hkkmass2_2->Draw("same");
         //   hkkmass2_2->SetLineColor(kBlack);
            prg4->Update();
            TCanvas *ckmass1= new TCanvas("ckmass1","",500,100,400,300);
            hkkmass1_1->Draw();
            hwkmass->Draw("same");
            hwkmass->SetLineColor(kRed);
            ckmass1->cd();
            TCanvas *ckmass2= new TCanvas("ckmass2","",500,100,400,300);
            hkkmass1_2->Draw();
            hwkmass->Draw("same");
            hwkmass->SetLineColor(kRed);
            */
            TH1F *hwrhomass11=(TH1F*)hwrhomass1->Clone();
            TH1F *hwrhomass111=(TH1F*)hwrhomass1->Clone();
            TH1F *hwrhomass1111=(TH1F*)hwrhomass1->Clone();
            TH1F *hwkmass11=(TH1F*)hwkmass1->Clone();
            TH1F *hwkmass111=(TH1F*)hwkmass1->Clone();
            TH1F *hwkmass1111=(TH1F*)hwkmass1->Clone();
            TCanvas *crhog= new TCanvas("crhog","Rho and K, set 1",500,100,400,300);
            crhog->cd();
            TPad *prg= new TPad("prg","rhomass1, set 1, good",0,0.50,0.50,1);
            prg->Draw();
            prg->cd();
            hrhomass1l->Draw();
            hkmass1l->Draw("same");
            hkmass1l->SetLineColor(kBlack);
            hwkmass1->Scale((double)hkmass1l->GetEntries()/(double)hwkmass1->GetEntries());
            hwkmass1->SetLineColor(kMagenta);
            hwkmass1->Draw("Histsame");
            hwrhomass1->SetLineColor(kRed);
            hwrhomass1->Scale((double)hrhomass1l->GetEntries()/(double)hwrhomass1->GetEntries());
            hwrhomass1->Draw("Histsame");
            prg->Update();
            crhog->cd();
            TPad *prg2= new TPad("prg2","rhomass2, set 1, good",0.50,0.50,1,1);
            prg2->Draw();
            prg2->cd();
            hrhomass1m->Draw();
            hkmass1m->Draw("same");
            hkmass1m->SetLineColor(kBlack);
            hwkmass11->Scale((double)hkmass1m->GetEntries()/(double)hwkmass11->GetEntries());
            hwkmass11->SetLineColor(kMagenta);
            hwkmass11->Draw("Histsame");
            hwrhomass11->SetLineColor(kRed);
            hwrhomass11->Scale((double)hrhomass1m->GetEntries()/(double)hwrhomass11->GetEntries());
            hwrhomass11->Draw("Histsame");
            prg2->Update();
            crhog->cd();
            TPad *prg3= new TPad("prg3","rhomass1, set 2, good",0,0,0.50,0.50);
            prg3->Draw();
            prg3->cd();
            hrhomass1h->Draw();
            hkmass1h->Draw("same");
            hkmass1h->SetLineColor(kBlack);
            hwkmass111->Scale((double)hkmass1h->GetEntries()/(double)hwkmass111->GetEntries());
            hwkmass111->SetLineColor(kMagenta);
            hwkmass111->Draw("Histsame");
            hwrhomass111->SetLineColor(kRed);
            hwrhomass111->Scale((double)hrhomass1h->GetEntries()/(double)hwrhomass111->GetEntries());
            hwrhomass111->Draw("Histsame");
            prg3->Update();
            crhog->cd();
            TPad *prg4= new TPad("prg4","rhomass2, set 2, good",0.50,0,1,0.5);
            prg4->Draw();
            prg4->cd();
            hrhomass1vh->Draw();
            hkmass1vh->Draw("same");
            hkmass1vh->SetLineColor(kBlack);
            hwkmass1111->Scale((double)hkmass1vh->GetEntries()/(double)hwkmass1111->GetEntries());
            hwkmass1111->SetLineColor(kMagenta);
            hwkmass1111->Draw("Histsame");
            hwrhomass1111->SetLineColor(kRed);
            hwrhomass1111->Scale((double)hrhomass1vh->GetEntries()/(double)hwrhomass1111->GetEntries());
            hwrhomass1111->Draw("Histsame");
            prg4->Update();
            crhog->cd();
            
            
            
            
            
            
            TH1F *hwrhomass22=(TH1F*)hwrhomass2->Clone();
            TH1F *hwrhomass222=(TH1F*)hwrhomass2->Clone();
            TH1F *hwrhomass2222=(TH1F*)hwrhomass2->Clone();
            TH1F *hwkmass22=(TH1F*)hwkmass2->Clone();
            TH1F *hwkmass222=(TH1F*)hwkmass2->Clone();
            TH1F *hwkmass2222=(TH1F*)hwkmass2->Clone();
            TCanvas *crhog2= new TCanvas("crhog2","Rho and K, set 2",500,100,400,300);
            crhog2->cd();
            TPad *prg5= new TPad("prg5","rhomass1, set 1, good",0,0.50,0.50,1);
            prg5->Draw();
            prg5->cd();
            hrhomass2l->Draw();
            hkmass2l->Draw("same");
            hkmass2l->SetLineColor(kBlack);
            hwkmass2->Scale((double)hkmass2l->GetEntries()/(double)hwkmass2->GetEntries());
            hwkmass2->Draw("Histsame");
            hwkmass2->SetLineColor(kMagenta);
            hwrhomass2->SetLineColor(kRed);
            hwrhomass2->Scale((double)hrhomass2l->GetEntries()/(double)hwrhomass2->GetEntries());
            hwrhomass2->Draw("Histsame");
            prg5->Update();
            crhog2->cd();
            TPad *prg6= new TPad("prg6","rhomass2, set 1, good",0.50,0.50,1,1);
            prg6->Draw();
            prg6->cd();
            hrhomass2m->Draw();
            hkmass2m->Draw("same");
            hkmass2m->SetLineColor(kBlack);
            hwkmass22->Scale((double)hkmass2m->GetEntries()/(double)hwkmass22->GetEntries());
            hwkmass22->Draw("Histsame");
            hwkmass22->SetLineColor(kMagenta);
            hwrhomass22->SetLineColor(kRed);
            hwrhomass22->Scale((double)hrhomass2m->GetEntries()/(double)hwrhomass22->GetEntries());
            hwrhomass22->Draw("Histsame");
            prg6->Update();
            crhog2->cd();
            TPad *prg7= new TPad("prg7","rhomass1, set 2, good",0,0,0.50,0.50);
            prg7->Draw();
            prg7->cd();
            hrhomass2h->Draw();
            hkmass2h->Draw("same");
            hkmass2h->SetLineColor(kBlack);
            hwkmass222->Scale((double)hkmass2h->GetEntries()/(double)hwkmass222->GetEntries());
            hwkmass222->Draw("Histsame");
            hwkmass222->SetLineColor(kMagenta);
            hwrhomass222->SetLineColor(kRed);
            hwrhomass222->Scale((double)hrhomass2h->GetEntries()/(double)hwrhomass222->GetEntries());
            hwrhomass222->Draw("Histsame");
            prg7->Update();
            crhog2->cd();
            TPad *prg8= new TPad("prg8","rhomass2, set 2, good",0.50,0,1,0.5);
            prg8->Draw();
            prg8->cd();
            hrhomass2vh->Draw();
            hkmass2vh->Draw("same");
            hkmass2vh->SetLineColor(kBlack);
            hwkmass2222->Scale((double)hkmass2vh->GetEntries()/(double)hwkmass2222->GetEntries());
            hwkmass2222->Draw("Histsame");
            hwkmass2222->SetLineColor(kMagenta);
            hwrhomass2222->SetLineColor(kRed);
            hwrhomass2222->Scale((double)hrhomass2vh->GetEntries()/(double)hwrhomass2222->GetEntries());
            hwrhomass2222->Draw("Histsame");
            prg8->Update();
            crhog2->cd();
            
            TH1F *hwkmassx=(TH1F*)hwkmass1->Clone();
            TH1F *hwkmassy=(TH1F*)hwkmass2->Clone();
            TCanvas *ckmass1= new TCanvas("ckmass1","",500,100,400,300);
            hkkmasstot1->Draw();
            hwkmassx->Draw("Histsame");
            hwkmassx->SetLineColor(kMagenta);
            ckmass1->cd();
            TCanvas *ckmass2= new TCanvas("ckmass2","",500,100,400,300);
            hkkmasstot2->Draw();
            hwkmassy->Draw("Histsame");
            hwkmassy->SetLineColor(kMagenta);
            
        }
        else if(var=="angle"){
            TCanvas *ca1g= new TCanvas("ca1g","Angle1, All sets, Good",500,100,400,300);
            ca1g->cd();
            TPad *pag1= new TPad("pag1","eta1, set 1, good",0,0.50,0.50,1);
            pag1->Draw();
            pag1->cd();
            h1_a1g->Draw();
            pag1->Update();
            ca1g->cd();
            TPad *pag2= new TPad("pag2","eta1, set 2, good",0.50,0.50,1,1);
            pag2->Draw();
            pag2->cd();
            h2_a1g->Draw();
            pag2->Update();
            ca1g->cd();
            TPad *pag3= new TPad("pag3","eta1, set 3, good",0,0,0.50,0.50);
            pag3->Draw();
            pag3->cd();
            h3_a1g->Draw();
            pag3->Update();
            ca1g->cd();
            TPad *pag4= new TPad("pag4","eta1, set 4, good",0.50,0,1,0.5);
            pag4->Draw();
            pag4->cd();
            h4_a1g->Draw();
            pag4->Update();
            
            TCanvas *ca2g= new TCanvas("ca2g","Angle2, All sets, Good",500,100,400,300);
            ca2g->cd();
            TPad *pag5= new TPad("pag5","eta2, set 1, good",0,0.50,0.50,1);
            pag5->Draw();
            pag5->cd();
            h1_a2g->Draw();
            pag5->Update();
            ca2g->cd();
            TPad *pag6= new TPad("pag6","eta2, set 2, good",0.50,0.50,1,1);
            pag6->Draw();
            pag6->cd();
            h2_a2g->Draw();
            pag6->Update();
            ca2g->cd();
            TPad *pag7= new TPad("pag7","eta2, set 3, good",0,0,0.50,0.50);
            pag7->Draw();
            pag7->cd();
            h3_a2g->Draw();
            pag7->Update();
            ca2g->cd();
            TPad *pag8= new TPad("pag8","eta2, set 4, good",0.50,0,1,0.5);
            pag8->Draw();
            pag8->cd();
            h4_a2g->Draw();
            pag8->Update();
            
            TCanvas *ca1b= new TCanvas("ca1b","Angle1, All sets, Bad",500,100,400,300);
            ca1b->cd();
            TPad *pab1= new TPad("pab1","eta1, set 1, bad",0,0.50,0.50,1);
            pab1->Draw();
            pab1->cd();
            h1_a1b->SetLineColor(kRed);
            h1_a1b->Draw();
            pab1->Update();
            ca1b->cd();
            TPad *pab2= new TPad("pab2","eta1, set 2, bad",0.50,0.50,1,1);
            pab2->Draw();
            pab2->cd();
            h2_a1b->SetLineColor(kRed);
            h2_a1b->Draw();
            pab2->Update();
            ca1b->cd();
            TPad *pab3= new TPad("pab3","eta1, set 3, bad",0,0,0.50,0.50);
            pab3->Draw();
            pab3->cd();
            h3_a1b->SetLineColor(kRed);
            h3_a1b->Draw();
            pab3->Update();
            ca1b->cd();
            TPad *pab4= new TPad("pab4","eta1, set 4, bad",0.50,0,1,0.5);
            pab4->Draw();
            pab4->cd();
            h4_a1b->SetLineColor(kRed);
            h4_a1b->Draw();
            pab4->Update();
            
            TCanvas *ca2b= new TCanvas("ca2b","Angle2, All sets, Bad",500,100,400,300);
            ca2b->cd();
            TPad *pab5= new TPad("pab5","angle2, set 1, bad",0,0.50,0.50,1);
            pab5->Draw();
            pab5->cd();
            h1_a2b->SetLineColor(kRed);
            h1_a2b->Draw();
            pab5->Update();
            ca2b->cd();
            TPad *pab6= new TPad("pab6","angle2, set 2, bad",0.50,0.50,1,1);
            pab6->Draw();
            pab6->cd();
            h2_a2b->SetLineColor(kRed);
            h2_a2b->Draw();
            pab6->Update();
            ca2b->cd();
            TPad *pab7= new TPad("pab7","angle2, set 3, bad",0,0,0.50,0.50);
            pab7->Draw();
            pab7->cd();
            h3_a2b->SetLineColor(kRed);
            h3_a2b->Draw();
            pab7->Update();
            ca2b->cd();
            TPad *pab8= new TPad("pab8","angle2, set 4, bad",0.50,0,1,0.5);
            pab8->Draw();
            pab8->cd();
            h4_a2b->SetLineColor(kRed);
            h4_a2b->Draw();
            pab8->Update();
        }
        else if(var=="right"){
            TCanvas *cright1= new TCanvas("cright1","",500,100,400,300);
            right1->Draw();
            TCanvas *cright2= new TCanvas("cright2","",500,100,400,300);
            right2->Draw();
            TCanvas *cright3= new TCanvas("cright3","",500,100,400,300);
            right3->Draw();
            TCanvas *cright4= new TCanvas("cright4","",500,100,400,300);
            right4->Draw();
        }
        else if(var=="cuts"){
            TCanvas *ccuts1= new TCanvas("ccuts1","",500,100,400,300);
            cuts1->Draw();
            TCanvas *ccuts2= new TCanvas("ccuts2","",500,100,400,300);
            cuts2->Draw();
            TCanvas *ccuts3= new TCanvas("ccuts3","",500,100,400,300);
            cuts3->Draw();
            TCanvas *ccuts4= new TCanvas("ccuts4","",500,100,400,300);
            cuts4->Draw();
        }
    }

}
