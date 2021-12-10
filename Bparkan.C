#define bpark_cxx
#include "bpark.h"
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

struct par{int ch;float phi;float pt;float mass;float angle;float eta;float theta;float en;float p;};

void bpark::Loop(){
    // TFile *out= new TFile("BparkhistosRight.root","RECREATE");
     TRandom2 *rand=new TRandom2();
     TH1F *h1_a1 = new TH1F("1_angle1","Angle between mu1 and same charge pion, set 1",300,-0.2,3.5);
     TH1F *h1_a2 = new TH1F("1_angle2","Angle between mu1 and diff charge pion, set 1",300,-3.5,3.5);

     TH1F *h3_a1 = new TH1F("3_angle1","Angle between mu1 and same charge pion, set 3",300,-0.2,3.5);
     TH1F *h3_a2 = new TH1F("3_angle2","Angle between mu1 and diff charge pion, set 3",300,-3.5,3.5);
     
     TH1F *hb1 = new TH1F("angleMB1","Angle between mu1 and B, set 1",300,-0.2,3.5);
     TH1F *hres1 = new TH1F("angleMR1","Angle between mu1 and res pi momentum, set 1",300,-0.2,3.5);
     TH1F *hresB1 = new TH1F("angleBR1","Angle between B and res pi momentum, set 1",300,-0.2,3.5);
     
     TH1F *hb3 = new TH1F("angleMB3","Angle between mu1 and B, set 3",300,-0.2,3.5);
     TH1F *hres3 = new TH1F("angleMR3","Angle between mu1 and res pi momentum, set 3",300,-0.2,3.5);
     TH1F *hresB3 = new TH1F("angleBR3","Angle between B and res pi momentum, set 3",300,-0.2,3.5);
     
     TH1F *hwh1 = new TH1F("hwangle1","Angle between mu1 and same charge pion",300,-0.2,3.5);
     TH1F *hwh2 = new TH1F("hwangle2","Angle between mu1 and diff charge pion",300,-3.5,3.5);
     
     TH1F *wh1 = new TH1F("wangle1","Angle between mu1 and same charge pion",300,-0.2,3.5);
     TH1F *wh2 = new TH1F("wangle2","Angle between mu1 and diff charge pion2",300,-3.5,3.5);
     
     TH1F *hwhb = new TH1F("hwangleMB","Angle between mu1 and B",300,-0.2,3.5);
     TH1F *hwhres = new TH1F("hwangleMR","Angle between mu1 and res pi momentum",300,-0.2,3.5);
     TH1F *hwhresB = new TH1F("hwangleBR","Angle between B and res pi momentum",300,-0.2,3.5);
     
     TH1F *whb = new TH1F("wangleMB","Angle between mu1 and B",300,-0.2,3.5);
     TH1F *whres = new TH1F("wangleMR","Angle between mu1 and res pi momentum",300,-0.2,3.5);
     TH1F *whresB = new TH1F("wangleBR","Angle between B and res pi momentum",300,-0.2,3.5);
     
     TH1F *hrate = new TH1F("rate","Rate: 1' cand vs 1' cand after checking all vs 2' and others cand ",2,-1,1);
     
     TH1F *hrhomass1_1 = new TH1F("hrhomass1_1","Rhomass1, set 1",300,-0.2,2);
     TH1F *hrhomass1_2 = new TH1F("hrhomass1_2","Rhomass2, set 1",300,-0.2,2);
     TH1F *whrhomass1_1 = new TH1F("whrhomass1_1","Rhomass",300,-0.2,2);
     TH1F *hrhomass1_1c = new TH1F("hrhomass1_1c","Rhomass1, set 1",300,-0.2,2);
     TH1F *hrhomass1_2c = new TH1F("hrhomass1_2c","Rhomass2, set 1",300,-0.2,2);
     
     TH1F *hrhomass3_1 = new TH1F("hrhomass3_1","Rhomass1, set 3",200,-0.2,2);
     TH1F *hrhomass3_2 = new TH1F("hrhomass3_2","Rhomass2, set 3",200,-0.2,2);
     TH1F *whrhomass3_1 = new TH1F("whrhomass3_1","Rhomass",200,-0.2,2);
     TH1F *hrhomass3_1c = new TH1F("hrhomass3_1c","Rhomass1,set 3",200,-0.2,2);
     TH1F *hrhomass3_2c = new TH1F("hrhomass3_2c","Rhomass2,set 3",200,-0.2,2);
     
     TH1F *heta1=new TH1F("heta1","Eta,set 1",300,-4,4);
     TH1F *heta3=new TH1F("heta3","Eta,set 3",300,-4,4);
     TH1F *wheta=new TH1F("wheta","",300,-4,4);
     TH1F *hwheta=new TH1F("hwheta","",300,-4,4);
     
     TH1F *hmet1_l = new TH1F("anglemeth1_l","Angle between mu1 and MET",90,-0.2,3.5);
     TH1F *hmet1_m = new TH1F("anglemet1_m","Angle between mu1 and MET",90,-0.2,3.5);
     TH1F *hmet1_h = new TH1F("anglemet1_h","Angle between mu1 and MET",90,-0.2,3.5);
     TH1F *hmet1_vh = new TH1F("anglemet1_vh","Angle between mu1 and MET",90,-0.2,3.5);
    
     TH1F *hmet3_l = new TH1F("anglemeth3_l","Angle between mu1 and MET",70,-0.2,3.5);
     TH1F *hmet3_m = new TH1F("anglemet3_m","Angle between mu1 and MET",70,-0.2,3.5);
     TH1F *hmet3_h = new TH1F("anglemet3_h","Angle between mu1 and MET",70,-0.2,3.5);
     TH1F *hmet3_vh = new TH1F("anglemet3_vh","Angle between mu1 and MET",70,-0.2,3.5);
    
     TH1F *hwhmet_l = new TH1F("wanglemeth_l","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *hwhmet_m = new TH1F("wanglemet_m","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *hwhmet_h = new TH1F("wanglemet_h","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *hwhmet_vh = new TH1F("wanglemet_vh","Angle between mu1 and MET",100,-0.2,3.5);
    
     TH1F *whmet_l = new TH1F("hwanglemeth_l","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *whmet_m = new TH1F("hwanglemet_m","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *whmet_h = new TH1F("hwanglemet_h","Angle between mu1 and MET",100,-0.2,3.5);
     TH1F *whmet_vh = new TH1F("hwanglemet_vh","Angle between mu1 and MET",100,-0.2,3.5);
   
     TH1F *hrig1 = new TH1F("rig1","right 1 set",2,0,1);
     TH1F *hrig2 = new TH1F("rig2","right 2 set",2,0,1);
     //  THF2 *******************************************************************************************
     
     TH2F *hangle1_1 = new TH2F("A1_1vsMpt","Angle1_1 vs mu1 pt,set 1",300,0,60,300,-0.2,3.5);
     TH2F *hangle1_2 = new TH2F("A1_2vsMpt","Angle1_2 vs mu1 pt,set 1",300,0,60,300,-3.5,3.5);
     
     TH2F *hangle3_1 = new TH2F("A3_1vsMpt","Angle3_1 vs mu1 pt,set 3",300,0,60,300,-0.2,3.5);
     TH2F *hangle3_2 = new TH2F("A3_2vsMpt","Angle3_2 vs mu1 pt,set 3",300,0,60,300,-3.5,3.5);
     
     TH2F *hangle1_1vs2 = new TH2F("1_A1vsA2","Angle1 vs Angle2,set 1",300,-0.2,3.5,300,-3.5,3.5);
     TH2F *hangle3_1vs2 = new TH2F("3_A1vsA2","Angle1 vs Angle2,set 3",300,-0.2,3.5,300,-3.5,3.5);
     
     TH2F *hwhangle1 = new TH2F("hwA1vsMpt","",300,0,60,300,-0.2,3.5);
     TH2F *hwhangle2 = new TH2F("hwA2vsMpt","",300,0,60,300,-3.5,3.5);
     TH2F *hwhangle2d = new TH2F("hwA1vsA2","",300,-0.2,3.5,300,-3.5,3.5);
     
     TH2F *whangle1 = new TH2F("wA1vsMpt","",300,0,60,300,-0.2,3.5);
     TH2F *whangle2 = new TH2F("wA2vsMpt","",300,0,60,300,-3.5,3.5);
     TH2F *whangle2d = new TH2F("wA1vsA2","",300,-0.2,3.5,300,-3.5,3.5);
     
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<< "nentries "<<nentries<<endl;
    par MU;
    par PI1;
    par PI2;
    par PI3;
    par B;
    par RESPI;
    par TAU;
    par RHO1;
    par RHO2;
    par MET;
    float rate=0;
    double ptpix=0;
    double ptpiy=0;
    float angleRB=0;
    float anglepi12=0;
    float anglepi13=0;
    double piresphi=0;
    float rhomass1=0;
    float rhomass2=0;
    float wrhomass=0;
    float theta3pi=0;
    float etamu3pi=0;
    int j=0;
    int numw=0;
    int numhw=0;
    int numc=0;
    float a=0;
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
        
        b_MET_phi->GetEntry(jentry);
        b_MET_et->GetEntry(jentry);
        
        b_BsTauTau_tau_isRight->GetEntry(jentry);
        j=0;
        a=0;
        for (int s=0; s<BsTauTau_tau_pi1_phi->size(); s++){
            PI1.ch=BsTauTau_tau_pi1_charge->at(s);
            PI2.ch=BsTauTau_tau_pi2_charge->at(s);
            PI3.ch=BsTauTau_tau_pi3_charge->at(s);
            MU.ch=BsTauTau_mu1_q->at(0);
            
          if ((MU.ch+PI1.ch+PI2.ch+PI3.ch)==0){
            if (PI1.ch==MU.ch) {                         // PION 1
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                B.phi=BsTauTau_B_phi->at(s);
                PI1.mass=BsTauTau_tau_pi1_mass->at(s);
                PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                PI1.pt=BsTauTau_tau_pi1_pt->at(s);
                PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                PI1.eta=BsTauTau_tau_pi1_eta->at(s);
                PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                MET.en=MET_et->at(0);
                MET.phi=MET_phi->at(0);
               }
            else if(PI2.ch==MU.ch){                         // PION 2
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                B.phi=BsTauTau_B_phi->at(s);
                PI1.mass=BsTauTau_tau_pi2_mass->at(s);
                PI2.mass=BsTauTau_tau_pi1_mass->at(s);
                PI3.mass=BsTauTau_tau_pi3_mass->at(s);
                PI1.pt=BsTauTau_tau_pi2_pt->at(s);
                PI2.pt=BsTauTau_tau_pi1_pt->at(s);
                PI3.pt=BsTauTau_tau_pi3_pt->at(s);
                PI1.eta=BsTauTau_tau_pi2_eta->at(s);
                PI2.eta=BsTauTau_tau_pi1_eta->at(s);
                PI3.eta=BsTauTau_tau_pi3_eta->at(s);
                MET.en=MET_et->at(0);
                MET.phi=MET_phi->at(0);
            }
            else if(PI3.ch==MU.ch){                         // PION 3
                MU.phi=BsTauTau_mu1_phi->at(0);
                MU.eta=BsTauTau_mu1_eta->at(0);
                PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                B.phi=BsTauTau_B_phi->at(s);
                PI1.mass=BsTauTau_tau_pi3_mass->at(s);
                PI2.mass=BsTauTau_tau_pi2_mass->at(s);
                PI3.mass=BsTauTau_tau_pi1_mass->at(s);
                PI1.pt=BsTauTau_tau_pi3_pt->at(s);
                PI2.pt=BsTauTau_tau_pi2_pt->at(s);
                PI3.pt=BsTauTau_tau_pi1_pt->at(s);
                PI1.eta=BsTauTau_tau_pi3_eta->at(s);
                PI2.eta=BsTauTau_tau_pi2_eta->at(s);
                PI3.eta=BsTauTau_tau_pi1_eta->at(s);
                MET.en=MET_et->at(0);
                MET.phi=MET_phi->at(0);
            }
            PI1.angle=ANG(MU.phi,PI1.phi);
            PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
            PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
            B.angle=ANG(B.phi,MU.phi);
            MU.theta=THETA(MU.eta);
            PI1.theta=THETA(PI1.eta);
            PI2.theta=THETA(PI2.eta);
            PI3.theta=THETA(PI3.eta);
            PI1.p=PI1.pt/sin(PI1.theta);
            PI2.p=PI2.pt/sin(PI2.theta);
            PI3.p=PI3.pt/sin(PI3.theta);
            
              //    PIONS WITH CORRECT CHARGES BUT WRONG ANGLE  (&& BsTauTau_tau_pi1_phi->size()>1 ?)
              if (PI1.angle>1.5){
                  
                  ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                  ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                  piresphi=atan2(ptpiy,ptpix);
                  
                  RESPI.angle=ANG(MU.phi,piresphi);
                  angleRB=ANG(B.phi,piresphi);
                  
                  anglepi12=ANG(PI1.phi,PI2.phi);
                  anglepi13=ANG(PI1.phi,PI3.phi);
                  
                  MET.angle=ANG(MU.phi,MET.phi);
                  if(MET.en<10) hwhmet_l->Fill(MET.angle);
                  else if(MET.en<20) hwhmet_m->Fill(MET.angle);
                  else if(MET.en<30) hwhmet_h->Fill(MET.angle);
                  else hwhmet_vh->Fill(MET.angle);
                  
                  theta3pi=acos(VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)/NORM(VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)) );
                  etamu3pi=MU.eta-ETA(theta3pi);
                  
                  cout<<etamu3pi<<endl;
                  hwh1->Fill(PI1.angle);
                  hwhangle1->Fill(BsTauTau_mu1_pt->at(0),PI1.angle);
                  hwh2->Fill(PI2.angle);
                  hwhangle2->Fill(BsTauTau_mu1_pt->at(0),PI2.angle);
                  hwh2->Fill(PI3.angle);
                  hwhangle2->Fill(BsTauTau_mu1_pt->at(0),PI3.angle);
                  hwhb->Fill(B.angle);
                  hwhres->Fill(RESPI.angle);
                  hwhresB->Fill(angleRB);
                  hwhangle2d->Fill(PI1.angle,PI2.angle);
                  hwhangle2d->Fill(PI1.angle,PI3.angle);
                  hwheta->Fill(etamu3pi);
                  numhw++;
                  continue;
            }
            else if(PI1.angle<1.5){
              //  RIGHT EVENTS ( 3 PIONS WITH OPPOSITE CHARGE OF MU AND ANGLE1 < 1.5 )
            if(j==0){
                if(s==0){
                    rate=-0.67;
                    ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                    ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                    
                    piresphi=atan2(ptpiy,ptpix);
                    RESPI.angle=ANG(MU.phi,piresphi);
                    angleRB=ANG(B.phi,piresphi);
  // **************************  ANGLE BETWEEN PIONS  ****************************************
                    anglepi12=ANG(PI1.phi,PI2.phi);
                    anglepi13=ANG(PI1.phi,PI3.phi);
  // **************************      RHOMASS       ****************************************
                    RHO1.mass=BsTauTau_tau_rhomass1->at(s);
                    RHO2.mass=BsTauTau_tau_rhomass2->at(s);
                    rhomass1=MASS(EN(PI1.mass,PI1.p,PI2.mass,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                    rhomass2=MASS(EN(PI1.mass,PI1.p,PI3.mass,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                    wrhomass=MASS(EN(PI3.mass,PI3.p,PI2.mass,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                    if(PI2.pt<PI3.pt){
                    hrhomass1_1c->Fill(rhomass1);
                    hrhomass1_2c->Fill(rhomass2);
                    }
                    else {
                    hrhomass1_1c->Fill(rhomass2);
                    hrhomass1_2c->Fill(rhomass1);
                    }
                    hrhomass1_1->Fill(RHO1.mass);
                    hrhomass1_2->Fill(RHO2.mass);
                    whrhomass1_1->Fill(wrhomass);
    // **************************   MET    ****************************************
                    MET.angle=ANG(MU.phi,MET.phi);
                    if(MET.en<10) hmet1_l->Fill(MET.angle);
                    else if(MET.en<20) hmet1_m->Fill(MET.angle);
                    else if(MET.en<30) hmet1_h->Fill(MET.angle);
                    else hmet1_vh->Fill(MET.angle);
                    
    // **************************    ETA   ****************************************
                    theta3pi=acos(VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)/NORM(VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)) );
                    etamu3pi=MU.eta-ETA(theta3pi);
                    
                    h1_a1->Fill(PI1.angle);
                    h1_a2->Fill(PI2.angle);
                    h1_a2->Fill(PI3.angle);
                    hangle1_1->Fill(BsTauTau_mu1_pt->at(0),PI1.angle);
                    hangle1_2->Fill(BsTauTau_mu1_pt->at(0),PI2.angle);
                    hangle1_2->Fill(BsTauTau_mu1_pt->at(0),PI3.angle);
                    hangle1_1vs2->Fill(PI1.angle,PI2.angle);
                    hangle1_1vs2->Fill(PI1.angle,PI3.angle);
                    hrate->Fill(rate);
                    hb1->Fill(B.angle);
                    hres1->Fill(RESPI.angle);
                    hresB1->Fill(angleRB);
                    heta1->Fill(etamu3pi);
                   // if(fabs(etamu3pi)<1.5){
                    if(BsTauTau_tau_isRight->at(s)==true) hrig1->Fill(0.67);
                    else hrig1->Fill(0.33);
                    numc++;
                    j=1;
                  //  }
                }
                else {
                    rate=0.67;
                    ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                    ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                      
                    piresphi=atan2(ptpiy,ptpix);
                    RESPI.angle=ANG(MU.phi,piresphi);
                    angleRB=ANG(B.phi,piresphi);
                      
    // **************************  ANGLE BETWEEN PIONS  ****************************************
                    anglepi12=ANG(PI1.phi,PI2.phi);
                    anglepi13=ANG(PI1.phi,PI3.phi);
    // **************************      RHOMASS       ****************************************
                    RHO1.mass=BsTauTau_tau_rhomass1->at(s);
                    RHO2.mass=BsTauTau_tau_rhomass2->at(s);
                    rhomass1=MASS(EN(PI1.mass,PI1.p,PI2.mass,PI2.p),VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,0,1));
                    rhomass2=MASS(EN(PI1.mass,PI1.p,PI3.mass,PI3.p),VECSUMX(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMY(PI1.pt,PI1.phi,PI3.pt,PI3.phi,0,0),VECSUMZ(PI1.pt,PI1.theta,PI3.pt,PI3.theta,0,1));
                    wrhomass=MASS(EN(PI3.mass,PI3.p,PI2.mass,PI2.p),VECSUMX(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMY(PI3.pt,PI3.phi,PI2.pt,PI2.phi,0,0),VECSUMZ(PI3.pt,PI3.theta,PI2.pt,PI2.theta,0,1));
                    if(PI2.pt<PI3.pt){
                    hrhomass3_1c->Fill(rhomass1);
                    hrhomass3_2c->Fill(rhomass2);
                    }
                    else {
                    hrhomass3_1c->Fill(rhomass2);
                    hrhomass3_2c->Fill(rhomass1);
                    }
                    hrhomass3_1->Fill(RHO1.mass);
                    hrhomass3_2->Fill(RHO2.mass);
                    whrhomass3_1->Fill(wrhomass);
   // **************************   MET    ****************************************
                    MET.angle=ANG(MU.phi,MET.phi);
                    if(MET.en<10) hmet3_l->Fill(MET.angle);
                    else if(MET.en<20) hmet3_m->Fill(MET.angle);
                    else if(MET.en<30) hmet3_h->Fill(MET.angle);
                    else hmet3_vh->Fill(MET.angle);
    // **************************    ETA   ****************************************
                    theta3pi=acos(VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)/NORM(VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)) );
                    etamu3pi=MU.eta-ETA(theta3pi);
                    
                    h3_a1->Fill(PI1.angle);
                    h3_a2->Fill(PI2.angle);
                    h3_a2->Fill(PI3.angle);
                    hrate->Fill(rate);
                    hangle3_1->Fill(BsTauTau_mu1_pt->at(0),PI1.angle);
                    hangle3_2->Fill(BsTauTau_mu1_pt->at(0),PI2.angle);
                    hangle3_2->Fill(BsTauTau_mu1_pt->at(0),PI3.angle);
                    hangle3_1vs2->Fill(PI1.angle,PI2.angle);
                    hangle3_1vs2->Fill(PI1.angle,PI3.angle);
                    hb3->Fill(B.angle);
                    hres3->Fill(RESPI.angle);
                    hresB3->Fill(angleRB);
                    heta3->Fill(etamu3pi);
                    if(BsTauTau_tau_isRight->at(s)==true) hrig2->Fill(0.67);
                    else hrig2->Fill(0.33);
                    j=1;
                    numc++;
                }
              }
            }
          }
            
            // WRONG TRIPLETS ( UNCORRECT CHARGES )
        
           //else if((PI1.ch+PI2.ch+PI3.ch)==1 || (PI1.ch+PI2.ch+PI3.ch)==-1){
            else {
                MU.phi=BsTauTau_mu1_phi->at(0);
                if (MU.ch==PI1.ch==PI2.ch){
                    float ran=rand->Rndm();
                    if(ran>0.5){
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                    }
                    else{
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi3_phi->at(s);
                    }
                }
                else if(MU.ch==PI1.ch==PI3.ch){
                    float ran=rand->Rndm();
                    if(ran>0.5){
                        PI1.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                    }
                    else{
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi1_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi2_phi->at(s);
                    }
                }
                else if(MU.ch==PI2.ch==PI3.ch){
                    float ran=rand->Rndm();
                    if(ran>0.5){
                        PI1.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                    }
                    else{
                        PI1.phi=BsTauTau_tau_pi3_phi->at(s);
                        PI2.phi=BsTauTau_tau_pi2_phi->at(s);
                        PI3.phi=BsTauTau_tau_pi1_phi->at(s);
                    }
                }
                PI1.angle=ANG(MU.phi,PI1.phi);
                PI2.angle=ANG(MU.phi,PI2.phi)*SGN(PI1.phi-MU.phi)*SGN(PI2.phi-MU.phi);
                PI3.angle=ANG(MU.phi,PI3.phi)*SGN(PI1.phi-MU.phi)*SGN(PI3.phi-MU.phi);
                ptpix=VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                ptpiy=VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi);
                
                piresphi=atan2(ptpiy,ptpix);
                RESPI.angle=ANG(MU.phi,piresphi);
                angleRB=ANG(B.phi,piresphi);
                
                MET.angle=ANG(MU.phi,MET.phi);
                if(MET.en<10) whmet_l->Fill(MET.angle);
                else if(MET.en<20) whmet_m->Fill(MET.angle);
                else if(MET.en<30) whmet_h->Fill(MET.angle);
                else whmet_vh->Fill(MET.angle);
                
                theta3pi=acos(VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)/NORM(VECSUMX(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMY(PI1.pt,PI1.phi,PI2.pt,PI2.phi,PI3.pt,PI3.phi),VECSUMZ(PI1.pt,PI1.theta,PI2.pt,PI2.theta,PI3.pt,PI3.theta)) );
                etamu3pi=MU.eta-ETA(theta3pi);
                
             //   MET.angle=ANG(MU.phi,MET.phi);
              //  whmet->Fill(MET.angle);
                
                wh1->Fill(PI1.angle);
                whangle1->Fill(BsTauTau_mu1_pt->at(0),PI1.angle);
                wh2->Fill(PI2.angle);
                whangle2->Fill(BsTauTau_mu1_pt->at(0),PI2.angle);
                wh2->Fill(PI3.angle);
                whangle2->Fill(BsTauTau_mu1_pt->at(0),PI3.angle);
                whb->Fill(B.angle);
                whres->Fill(RESPI.angle);
                whresB->Fill(angleRB);
                whangle2d->Fill(PI1.angle,PI2.angle);
                whangle2d->Fill(PI1.angle,PI3.angle);
                wheta->Fill(etamu3pi);
                numw++;
            }
        }
     }
    
  // out->Write();
  // out->Close();

    cout<<"Done!"<<"Correct candidates: "<<numc<<endl;
    cout<<"Wrong candidates: "<<numw<<endl;
    cout<<"Half Wrong candidates: "<<numhw<<endl;
    cout<<"List of variables: eta | met | rho | angle | respi"<<endl;
    cout<<"Insert variable name: "<<endl;
   // Draw results
    TString var;
    while(1){
        cin>>var;
        if (var=="q") break;
        else if(var=="eta"){
            TCanvas *ceta1= new TCanvas("ceta1","",500,100,400,300);
            wheta->Scale((float)numc/numw);
            hwheta->Scale((float)numc/numhw);
            heta1->Draw();
            wheta->SetLineColor(kRed);
            hwheta->SetLineColor(kGreen);
            wheta->Draw("Histsame");
            hwheta->Draw("Histsame");
            heta1->GetYaxis()->SetRangeUser(0.,200.);
            ceta1->cd();
            TCanvas *ceta3= new TCanvas("ceta3","",500,100,400,300);
            heta3->Draw();
            ceta3->cd();
            TCanvas *cr= new TCanvas("cr","",500,100,400,300);
            hrate->Draw();
            cr->cd();
        }
        else if(var=="rho"){
            TCanvas *crho1_1= new TCanvas("crho1_1","",500,100,400,300);
            hrhomass1_1->Draw();
            hrhomass1_1->SetLineColor(kBlue);
            hrhomass1_1c->Draw("same");
            hrhomass1_1c->SetLineColor(kRed);
            whrhomass1_1->Draw("same");
            whrhomass1_1->SetLineColor(kGreen);
            crho1_1->cd();
            TCanvas *crho1_2= new TCanvas("crho1_2","",500,100,400,300);
            hrhomass1_2c->Draw();
            hrhomass1_2c->SetLineColor(kRed);
            whrhomass1_1->Draw("same");
            whrhomass1_1->SetLineColor(kGreen);
            hrhomass1_2->Draw("same");
            hrhomass1_2->SetLineColor(kBlue);
            crho1_2->cd();
            TCanvas *crho3_1= new TCanvas("crho3_1","",500,100,400,300);
            hrhomass3_1->Draw();
            hrhomass3_1->SetLineColor(kBlue);
            hrhomass3_1c->Draw("same");
            hrhomass3_1c->SetLineColor(kRed);
            whrhomass3_1->Draw("same");
            whrhomass3_1->SetLineColor(kGreen);
            crho3_1->cd();
            TCanvas *crho3_2= new TCanvas("crho3_2","",500,100,400,300);
            hrhomass3_2c->Draw();
            hrhomass3_2c->GetYaxis()->SetRangeUser(0.,2500);
            hrhomass3_2c->SetLineColor(kRed);
            whrhomass3_1->Draw("same");
            whrhomass3_1->SetLineColor(kGreen);
            hrhomass3_2->Draw("same");
            hrhomass3_2->SetLineColor(kBlue);
            crho3_2->cd();
            TCanvas *cr= new TCanvas("cr","",500,100,400,300);
            hrate->Draw();
            cr->cd();
        }
        else if(var=="met"){
            TCanvas *cmet1= new TCanvas("cmet1","",500,100,400,300);
            cmet1->cd();
            TPad *p1= new TPad("pad1","met low",0.01,0.51,0.49,0.99);
            p1->Draw();
            p1->cd();
            hmet1_l->Draw();
       //     whmet_l->SetLineColor(kRed);
        //    whmet_l->Draw("same");
        //    hwhmet_l->SetLineColor(kGreen);
        //    hwhmet_l->Draw("same");
            p1->Update();
            cmet1->cd();
            TPad *p2= new TPad("pad2","met med",0.51,0.51,0.99,0.99);
            p2->Draw();
            p2->cd();
            hmet1_m->Draw();
        //    whmet_m->SetLineColor(kRed);
       //     whmet_m->Draw("same");
       //     hwhmet_m->SetLineColor(kGreen);
       //     hwhmet_m->Draw("same");
            p2->Update();
            cmet1->cd();
            TPad *p3= new TPad("pad3","met high",0.01,0.01,0.49,0.49);
            p3->Draw();
            p3->cd();
            hmet1_h->Draw();
       //     whmet_h->SetLineColor(kRed);
       //     whmet_h->Draw("same");
       //     hwhmet_h->SetLineColor(kGreen);
       //     hwhmet_h->Draw("same");
            p3->Update();
            cmet1->cd();
            TPad *p4= new TPad("pad4","met very high",0.51,0.01,0.99,0.49);
            p4->Draw();
            p4->cd();
            hmet1_vh->Draw();
       //     whmet_vh->SetLineColor(kRed);
       //     whmet_vh->Draw("same");
       //     hwhmet_vh->SetLineColor(kGreen);
       //     hwhmet_vh->Draw("same");
            p4->Update();
            TCanvas *cmet2= new TCanvas("cmet2","",500,100,400,300);
            cmet2->cd();
            TPad *pp1= new TPad("ppad1","met low",0.01,0.51,0.49,0.99);
            pp1->Draw();
            pp1->cd();
            hmet3_l->Draw();
            pp1->Update();
            cmet2->cd();
            TPad *pp2= new TPad("ppad2","met med",0.51,0.51,0.99,0.99);
            pp2->Draw();
            pp2->cd();
            hmet3_m->Draw();
            pp2->Update();
            cmet2->cd();
            TPad *pp3= new TPad("ppad3","met high",0.01,0.01,0.49,0.49);
            pp3->Draw();
            pp3->cd();
            hmet3_h->Draw();
            pp3->Update();
            cmet2->cd();
            TPad *pp4= new TPad("ppad4","met very high",0.51,0.01,0.99,0.49);
            pp4->Draw();
            pp4->cd();
            hmet3_vh->Draw();
            pp4->Update();
        }
        else if(var=="angle"){
            TCanvas *c1_1= new TCanvas("c1_1","",500,100,400,300);
            h1_a1->Draw();
            wh1->Scale((float)numc/numw);
            hwh1->Scale((float)numc/numhw);
            wh1->SetLineColor(kRed);
            wh1->Draw("Histsame");
            hwh1->SetLineColor(kGreen);
            hwh1->Draw("Histsame");
            c1_1->cd();
            
            TCanvas *c3_1= new TCanvas("c3_1","",500,100,400,300);
            h3_a1->Draw();
         //   wh1->SetLineColor(kRed);
         //   wh1->Draw("same");
         //   hwh1->SetLineColor(kGreen);
         //   hwh1->Draw("same");
            c3_1->cd();
         
            TCanvas *c1_2= new TCanvas("c1_2","",500,100,400,300);
            h1_a2->Draw();
            wh2->SetLineColor(kRed);
            wh2->Draw("same");
            hwh2->SetLineColor(kGreen);
            hwh2->Draw("same");
            c1_2->cd();
            
            TCanvas *c3_2= new TCanvas("c3_2","",500,100,400,300);
            h3_a2->Draw();
          //  wh2->SetLineColor(kRed);
         //   wh2->Draw("same");
         //   hwh2->SetLineColor(kGreen);
         //   hwh2->Draw("same");
            c3_2->cd();
           
            TCanvas *c2d1_1= new TCanvas("c2d1_1","",500,100,400,300);
            hangle1_1->Draw();
            whangle1->Scale((float)numc/numw);
            hwhangle1->Scale((float)numc/numhw);
            hangle1_1->SetMarkerStyle(7);
            whangle1->SetMarkerStyle(7);
            whangle1->SetMarkerColor(kRed);
            whangle1->Draw("same");
            hwhangle1->SetMarkerStyle(7);
            hwhangle1->SetMarkerColor(kGreen);
            hwhangle1->Draw("same");
            c2d1_1->cd();
            
            TCanvas *c2d1_2= new TCanvas("c2d1_2","",500,100,400,300);
            hangle1_2->Draw();
            whangle2->Scale((float)numc/numw);
            hwhangle2->Scale((float)numc/numhw);
            hangle1_2->SetMarkerStyle(7);
            whangle2->SetMarkerStyle(7);
            whangle2->SetMarkerColor(kRed);
            whangle2->Draw("same");
            hwhangle2->SetMarkerStyle(7);
            hwhangle2->SetMarkerColor(kGreen);
            hwhangle2->Draw("same");
            c2d1_2->cd();

            TCanvas *c2d3_1= new TCanvas("c2d3_1","",500,100,400,300);
            hangle3_1->Draw();
            hangle3_1->SetMarkerStyle(7);
          //  whangle1->SetMarkerStyle(7);
          //  whangle1->SetMarkerColor(kRed);
          //  whangle1->Draw("same");
           // hwhangle1->SetMarkerStyle(7);
           // hwhangle1->SetMarkerColor(kGreen);
          //  hwhangle1->Draw("same");
            c2d3_1->cd();
           
            TCanvas *c2d3_2= new TCanvas("c2d3_2","",500,100,400,300);
            hangle3_2->Draw();
            hangle3_2->SetMarkerStyle(7);
          //  whangle2->SetMarkerStyle(7);
          //  whangle2->SetMarkerColor(kRed);
          //  whangle2->Draw("same");
          //  hwhangle2->SetMarkerStyle(7);
           // hwhangle2->SetMarkerColor(kGreen);
          //  hwhangle2->Draw("same");
            c2d3_2->cd();
            
            TCanvas *cr= new TCanvas("cr","",500,100,400,300);
            hrate->Draw();
            cr->cd();
    // ******************************************  A1vsA2 *******************************
            
            TCanvas *c1_1vs2= new TCanvas("c1_1vs2","",500,100,400,300);
            hangle1_1vs2->Draw();
            whangle2d->Scale((float)numc/numw);
            hwhangle2d->Scale((float)numc/numhw);
            hangle1_1vs2->SetMarkerStyle(7);
            whangle2d->SetMarkerStyle(7);
            whangle2d->SetMarkerColor(kRed);
            whangle2d->Draw("same");
            hwhangle2d->SetMarkerStyle(7);
            hwhangle2d->SetMarkerColor(kGreen);
            hwhangle2d->Draw("same");
            c1_1vs2->cd();
            
            TCanvas *c3_1vs2= new TCanvas("c3_1vs2","",500,100,400,300);
            hangle3_1vs2->Draw();
            hangle3_1vs2->SetMarkerStyle(7);
         //   whangle2d->SetMarkerStyle(7);
          //  whangle2d->SetMarkerColor(kRed);
          //  whangle2d->Draw("same");
          //  hwhangle2d->SetMarkerStyle(7);
         //   hwhangle2d->SetMarkerColor(kGreen);
          //  hwhangle2d->Draw("same");
            c3_1vs2->cd();
        }
 // ************************************  ANGLE B; ANGLE RB ; ANGLE RM ***********************************
        else if(var=="respi"){
            TCanvas *cb1= new TCanvas("cb1","",500,100,400,300);
            hb1->Draw();
            whb->Scale((float)numc/numw);
            hwhb->Scale((float)numc/numhw);
            whb->SetLineColor(kRed);
            whb->Draw("Histsame");
            hwhb->SetLineColor(kGreen);
            hwhb->Draw("Histsame");
            cb1->cd();
         
            TCanvas *crm1= new TCanvas("crm1","",500,100,400,300);
            hres1->Draw();
            whres->Scale((float)numc/numw);
            hwhres->Scale((float)numc/numhw);
            whres->SetLineColor(kRed);
            whres->Draw("Histsame");
            hwhres->SetLineColor(kGreen);
            hwhres->Draw("Histsame");
            crm1->cd();
         
            TCanvas *crb1= new TCanvas("crb1","",500,100,400,300);
            hresB1->Draw();
            whresB->Scale((float)numc/numw);
            hwhresB->Scale((float)numc/numhw);
            whresB->SetLineColor(kRed);
            whresB->Draw("Histsame");
            hwhresB->SetLineColor(kGreen);
            hwhresB->Draw("Histsame");
            crb1->cd();
            
            TCanvas *cb3= new TCanvas("cb3","",500,100,400,300);
            hb3->Draw();
           // whb->SetLineColor(kRed);
           // whb->Draw("same");
           // hwhb->SetLineColor(kGreen);
          //  hwhb->Draw("same");
            cb3->cd();
         
            TCanvas *crm3= new TCanvas("crm3","",500,100,400,300);
            hres3->Draw();
           // whres->SetLineColor(kRed);
           // whres->Draw("same");
           // hwhres->SetLineColor(kGreen);
           // hwhres->Draw("same");
            crm3->cd();
         
            TCanvas *crb3= new TCanvas("crb3","",500,100,400,300);
            hresB3->Draw();
          //  whresB->SetLineColor(kRed);
         //   whresB->Draw("same");
         //   hwhresB->SetLineColor(kGreen);
           // hwhresB->Draw("same");
            crb3->cd();
        }
        else if(var=="right"){
            TCanvas *cright1= new TCanvas("cright1","",500,100,400,300);
            hrig1->Draw();
     
            cright1->cd();
            TCanvas *cright2= new TCanvas("cright2","",500,100,400,300);
            hrig2->Draw();
           
            cright2->cd();
        }
    }
}
