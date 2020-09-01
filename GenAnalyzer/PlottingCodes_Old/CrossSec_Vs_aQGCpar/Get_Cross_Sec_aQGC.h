//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 24 13:53:02 2017 by ROOT version 6.06/01
// from TTree otree/GenParticles Basic Info
// found on file: /eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root
//////////////////////////////////////////////////////////

#ifndef Get_Cross_Sec_aQGC_h
#define Get_Cross_Sec_aQGC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class Get_Cross_Sec_aQGC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxngen_WJet1 = 1;
   const Int_t kMaxngen_WJet2 = 1;
   const Int_t kMaxngen_VBFjet1 = 1;
   const Int_t kMaxngen_VBFjet2 = 1;

   // Declaration of leaf types
   Int_t           isMuMinus;
   Double_t        LHELeptPt;
   Double_t        LHELeptEta;
   Double_t        LHELeptPhi;
   Double_t        LHELeptM;
   Double_t        LHELeptE;
   Double_t        LHENuPt;
   Double_t        LHENuEta;
   Double_t        LHENuPhi;
   Double_t        LHENuM;
   Double_t        LHENuE;
   Double_t        LHE_DeltaM_Wqrk0_pt;
   Double_t        LHE_DeltaM_Wqrk0_eta;
   Double_t        LHE_DeltaM_Wqrk0_phi;
   Double_t        LHE_DeltaM_Wqrk0_M;
   Double_t        LHE_DeltaM_Wqrk0_E;
   Double_t        LHE_DeltaM_Wqrk0_Mt;
   Double_t        LHE_DeltaM_Wqrk1_pt;
   Double_t        LHE_DeltaM_Wqrk1_eta;
   Double_t        LHE_DeltaM_Wqrk1_phi;
   Double_t        LHE_DeltaM_Wqrk1_M;
   Double_t        LHE_DeltaM_Wqrk1_E;
   Double_t        LHE_DeltaM_Wqrk1_Mt;
   Double_t        LHE_DeltaM_Iqrk0_pt;
   Double_t        LHE_DeltaM_Iqrk0_eta;
   Double_t        LHE_DeltaM_Iqrk0_phi;
   Double_t        LHE_DeltaM_Iqrk0_E;
   Double_t        LHE_DeltaM_Iqrk0_M;
   Double_t        LHE_DeltaM_Iqrk0_Mt;
   Double_t        LHE_DeltaM_Iqrk1_pt;
   Double_t        LHE_DeltaM_Iqrk1_eta;
   Double_t        LHE_DeltaM_Iqrk1_phi;
   Double_t        LHE_DeltaM_Iqrk1_E;
   Double_t        LHE_DeltaM_Iqrk1_M;
   Double_t        LHE_DeltaM_Iqrk1_Mt;
   Double_t        LHE_DeltaM_mWW;
   Double_t        LHE_DeltaM_mtWW;
   Double_t        LHE_DeltaM_mWLep;
   Double_t        LHE_DeltaM_mtWLep;
   Double_t        LHE_DeltaM_mWHad;
   Double_t        LHE_DeltaM_mtWHad;
   Double_t        LHE_DeltaM_costheta1;
   Double_t        LHE_DeltaM_costheta2;
   Double_t        LHE_DeltaM_phi;
   Double_t        LHE_DeltaM_costhetastar;
   Double_t        LHE_DeltaM_phi1;
   Double_t        LHE_DeltaM_dEtajj;
   Double_t        LHE_DeltaM_dPhijj;
   Double_t        LHE_DeltaM_mjj;
   Double_t        LHE_DeltaM_VBSCentrality;
   Double_t        LHE_MothInfo_Wqrk0_pt;
   Double_t        LHE_MothInfo_Wqrk0_eta;
   Double_t        LHE_MothInfo_Wqrk0_phi;
   Double_t        LHE_MothInfo_Wqrk0_M;
   Double_t        LHE_MothInfo_Wqrk0_E;
   Double_t        LHE_MothInfo_Wqrk0_Mt;
   Double_t        LHE_MothInfo_Wqrk1_pt;
   Double_t        LHE_MothInfo_Wqrk1_eta;
   Double_t        LHE_MothInfo_Wqrk1_phi;
   Double_t        LHE_MothInfo_Wqrk1_M;
   Double_t        LHE_MothInfo_Wqrk1_E;
   Double_t        LHE_MothInfo_Wqrk1_Mt;
   Double_t        LHE_MothInfo_Iqrk0_pt;
   Double_t        LHE_MothInfo_Iqrk0_eta;
   Double_t        LHE_MothInfo_Iqrk0_phi;
   Double_t        LHE_MothInfo_Iqrk0_E;
   Double_t        LHE_MothInfo_Iqrk0_M;
   Double_t        LHE_MothInfo_Iqrk0_Mt;
   Double_t        LHE_MothInfo_Iqrk1_pt;
   Double_t        LHE_MothInfo_Iqrk1_eta;
   Double_t        LHE_MothInfo_Iqrk1_phi;
   Double_t        LHE_MothInfo_Iqrk1_E;
   Double_t        LHE_MothInfo_Iqrk1_M;
   Double_t        LHE_MothInfo_Iqrk1_Mt;
   Double_t        LHE_MothInfo_mWW;
   Double_t        LHE_MothInfo_mtWW;
   Double_t        LHE_MothInfo_mWLep;
   Double_t        LHE_MothInfo_mtWLep;
   Double_t        LHE_MothInfo_mWHad;
   Double_t        LHE_MothInfo_mtWHad;
   Double_t        LHE_MothInfo_costheta1;
   Double_t        LHE_MothInfo_costheta2;
   Double_t        LHE_MothInfo_phi;
   Double_t        LHE_MothInfo_costhetastar;
   Double_t        LHE_MothInfo_phi1;
   Double_t        LHE_MothInfo_dEtajj;
   Double_t        LHE_MothInfo_dPhijj;
   Double_t        LHE_MothInfo_mjj;
   Double_t        LHE_MothInfo_VBSCentrality;
   vector<string>  *LHEWeightIDs;
   vector<double>  *LHEWeights;
   Int_t           ngen_Lept;
   Double_t        gen_LeptPt;
   Double_t        gen_LeptEta;
   Double_t        gen_LeptPhi;
   Double_t        gen_LeptM;
   Int_t           gen_LeptStatus;
   Int_t           gen_LeptId;
   Double_t        gen_LeptMother;
   Int_t           gen_LeptGrandMother;
   Int_t           gen_NuPdgId;
   Int_t           ngen_Nu;
   Double_t        gen_NuPt;
   Double_t        gen_NuEta;
   Double_t        gen_NuPhi;
   Double_t        gen_NuM;
   Double_t        gen_NuQ;
   Int_t           gen_Nustatus;
   Int_t           gen_NuMother;
   Int_t           gen_NuGrandMother;
   Int_t           gen_WJet1_PdgId;
   Int_t           ngen_WJet1_;
   Double_t        gen_WJet1_Pt;
   Double_t        gen_WJet1_Eta;
   Double_t        gen_WJet1_Phi;
   Double_t        gen_WJet1_M;
   Double_t        gen_WJet1_E;
   Double_t        gen_WJet1_Q;
   Int_t           gen_WJet1_status;
   Int_t           gen_WJet1_Mother;
   Int_t           gen_WJet1_GrandMother;
   Int_t           gen_WJet2_PdgId;
   Int_t           ngen_WJet2_;
   Double_t        gen_WJet2_Pt;
   Double_t        gen_WJet2_Eta;
   Double_t        gen_WJet2_Phi;
   Double_t        gen_WJet2_M;
   Double_t        gen_WJet2_E;
   Double_t        gen_WJet2_Q;
   Int_t           gen_WJet2_status;
   Int_t           gen_WJet2_Mother;
   Int_t           gen_WJet2_GrandMother;
   Int_t           gen_VBFjet1_PdgId;
   Int_t           ngen_VBFjet1_;
   Double_t        gen_VBFjet1_Pt;
   Double_t        gen_VBFjet1_Eta;
   Double_t        gen_VBFjet1_Phi;
   Double_t        gen_VBFjet1_M;
   Double_t        gen_VBFjet1_Q;
   Int_t           gen_VBFjet1_status;
   Int_t           gen_VBFjet1_Mother;
   Int_t           gen_VBFjet1_GrandMother;
   Int_t           gen_VBFjet2_PdgId;
   Int_t           ngen_VBFjet2_;
   Double_t        gen_VBFjet2_Pt;
   Double_t        gen_VBFjet2_Eta;
   Double_t        gen_VBFjet2_Phi;
   Double_t        gen_VBFjet2_M;
   Double_t        gen_VBFjet2_Q;
   Int_t           gen_VBFjet2_status;
   Int_t           gen_VBFjet2_Mother;
   Int_t           gen_VBFjet2_GrandMother;
   Double_t        gen_VBFjet1jet2_Pt;
   Double_t        gen_VBFjet1jet2_Eta;
   Double_t        gen_VBFjet1jet2_Phi;
   Double_t        gen_VBFjet1jet2_M;
   Double_t        gen_vbfjet_deltaR;
   Double_t        gen_WHad_Pt;
   Double_t        gen_WHad_M;
   Double_t        gen_WHad_Mt;
   Double_t        gen_WHad_deltaeta;
   Double_t        gen_WHad_deltaphi;
   Double_t        gen_WHad_deltar;
   Double_t        gen_deltaR_LepWHad;
   Double_t        gen_deltaphi_NuWHad;
   Double_t        gen_deltaphi_WlepWHad;
   Double_t        gen_mWW;
   Double_t        gen_mtWW;
   Double_t        gen_mWLep;
   Double_t        gen_mtWLep;
   Double_t        gen_mWHad;
   Double_t        gen_mtWHad;
   Double_t        gen_costheta1;
   Double_t        gen_costheta2;
   Double_t        gen_phi;
   Double_t        gen_costhetastar;
   Double_t        gen_phi1;
   Double_t        gen_dEtajj;
   Double_t        gen_dPhijj;
   Double_t        gen_mjj;
   Double_t        gen_VBSCentrality;
   vector<double>  *genQuarkStatus;
   Int_t           ngenJet;
   Int_t           nVBFJet;

   // List of branches
   TBranch        *b_isMuMinus;   //!
   TBranch        *b_LHELeptPt;   //!
   TBranch        *b_LHELeptEta;   //!
   TBranch        *b_LHELeptPhi;   //!
   TBranch        *b_LHELeptM;   //!
   TBranch        *b_LHELeptE;   //!
   TBranch        *b_LHENuPt;   //!
   TBranch        *b_LHENuEta;   //!
   TBranch        *b_LHENuPhi;   //!
   TBranch        *b_LHENuM;   //!
   TBranch        *b_LHENuE;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_pt;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_eta;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_phi;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_M;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_E;   //!
   TBranch        *b_LHE_DeltaM_Wqrk0_Mt;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_pt;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_eta;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_phi;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_M;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_E;   //!
   TBranch        *b_LHE_DeltaM_Wqrk1_Mt;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_pt;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_eta;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_phi;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_E;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_M;   //!
   TBranch        *b_LHE_DeltaM_Iqrk0_Mt;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_pt;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_eta;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_phi;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_E;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_M;   //!
   TBranch        *b_LHE_DeltaM_Iqrk1_Mt;   //!
   TBranch        *b_LHE_DeltaM_mWW;   //!
   TBranch        *b_LHE_DeltaM_mtWW;   //!
   TBranch        *b_LHE_DeltaM_mWLep;   //!
   TBranch        *b_LHE_DeltaM_mtWLep;   //!
   TBranch        *b_LHE_DeltaM_mWHad;   //!
   TBranch        *b_LHE_DeltaM_mtWHad;   //!
   TBranch        *b_LHE_DeltaM_costheta1;   //!
   TBranch        *b_LHE_DeltaM_costheta2;   //!
   TBranch        *b_LHE_DeltaM_phi;   //!
   TBranch        *b_LHE_DeltaM_costhetastar;   //!
   TBranch        *b_LHE_DeltaM_phi1;   //!
   TBranch        *b_LHE_DeltaM_dEtajj;   //!
   TBranch        *b_LHE_DeltaM_dPhijj;   //!
   TBranch        *b_LHE_DeltaM_mjj;   //!
   TBranch        *b_LHE_DeltaM_VBSCentrality;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_pt;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_eta;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_phi;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_M;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_E;   //!
   TBranch        *b_LHE_MothInfo_Wqrk0_Mt;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_pt;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_eta;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_phi;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_M;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_E;   //!
   TBranch        *b_LHE_MothInfo_Wqrk1_Mt;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_pt;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_eta;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_phi;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_E;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_M;   //!
   TBranch        *b_LHE_MothInfo_Iqrk0_Mt;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_pt;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_eta;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_phi;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_E;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_M;   //!
   TBranch        *b_LHE_MothInfo_Iqrk1_Mt;   //!
   TBranch        *b_LHE_MothInfo_mWW;   //!
   TBranch        *b_LHE_MothInfo_mtWW;   //!
   TBranch        *b_LHE_MothInfo_mWLep;   //!
   TBranch        *b_LHE_MothInfo_mtWLep;   //!
   TBranch        *b_LHE_MothInfo_mWHad;   //!
   TBranch        *b_LHE_MothInfo_mtWHad;   //!
   TBranch        *b_LHE_MothInfo_costheta1;   //!
   TBranch        *b_LHE_MothInfo_costheta2;   //!
   TBranch        *b_LHE_MothInfo_phi;   //!
   TBranch        *b_LHE_MothInfo_costhetastar;   //!
   TBranch        *b_LHE_MothInfo_phi1;   //!
   TBranch        *b_LHE_MothInfo_dEtajj;   //!
   TBranch        *b_LHE_MothInfo_dPhijj;   //!
   TBranch        *b_LHE_MothInfo_mjj;   //!
   TBranch        *b_LHE_MothInfo_VBSCentrality;   //!
   TBranch        *b_LHEWeightIDs;   //!
   TBranch        *b_LHEWeights;   //!
   TBranch        *b_ngen_Lept;   //!
   TBranch        *b_gen_LeptPt;   //!
   TBranch        *b_gen_LeptEta;   //!
   TBranch        *b_gen_LeptPhi;   //!
   TBranch        *b_gen_LeptM;   //!
   TBranch        *b_gen_LeptStatus;   //!
   TBranch        *b_gen_LeptId;   //!
   TBranch        *b_gen_LeptMother;   //!
   TBranch        *b_gen_LeptGrandMother;   //!
   TBranch        *b_gen_NuPdgId;   //!
   TBranch        *b_ngen_Nu;   //!
   TBranch        *b_gen_NuPt;   //!
   TBranch        *b_gen_NuEta;   //!
   TBranch        *b_gen_NuPhi;   //!
   TBranch        *b_gen_NuM;   //!
   TBranch        *b_gen_NuQ;   //!
   TBranch        *b_gen_Nustatus;   //!
   TBranch        *b_gen_NuMother;   //!
   TBranch        *b_gen_NuGrandMother;   //!
   TBranch        *b_gen_WJet1_PdgId;   //!
   TBranch        *b_ngen_WJet1_;   //!
   TBranch        *b_gen_WJet1_Pt;   //!
   TBranch        *b_gen_WJet1_Eta;   //!
   TBranch        *b_gen_WJet1_Phi;   //!
   TBranch        *b_gen_WJet1_M;   //!
   TBranch        *b_gen_WJet1_E;   //!
   TBranch        *b_gen_WJet1_Q;   //!
   TBranch        *b_gen_WJet1_status;   //!
   TBranch        *b_gen_WJet1_Mother;   //!
   TBranch        *b_gen_WJet1_GrandMother;   //!
   TBranch        *b_gen_WJet2_PdgId;   //!
   TBranch        *b_ngen_WJet2_;   //!
   TBranch        *b_gen_WJet2_Pt;   //!
   TBranch        *b_gen_WJet2_Eta;   //!
   TBranch        *b_gen_WJet2_Phi;   //!
   TBranch        *b_gen_WJet2_M;   //!
   TBranch        *b_gen_WJet2_E;   //!
   TBranch        *b_gen_WJet2_Q;   //!
   TBranch        *b_gen_WJet2_status;   //!
   TBranch        *b_gen_WJet2_Mother;   //!
   TBranch        *b_gen_WJet2_GrandMother;   //!
   TBranch        *b_gen_VBFjet1_PdgId;   //!
   TBranch        *b_ngen_VBFjet1_;   //!
   TBranch        *b_gen_VBFjet1_Pt;   //!
   TBranch        *b_gen_VBFjet1_Eta;   //!
   TBranch        *b_gen_VBFjet1_Phi;   //!
   TBranch        *b_gen_VBFjet1_M;   //!
   TBranch        *b_gen_VBFjet1_Q;   //!
   TBranch        *b_gen_VBFjet1_status;   //!
   TBranch        *b_gen_VBFjet1_Mother;   //!
   TBranch        *b_gen_VBFjet1_GrandMother;   //!
   TBranch        *b_gen_VBFjet2_PdgId;   //!
   TBranch        *b_ngen_VBFjet2_;   //!
   TBranch        *b_gen_VBFjet2_Pt;   //!
   TBranch        *b_gen_VBFjet2_Eta;   //!
   TBranch        *b_gen_VBFjet2_Phi;   //!
   TBranch        *b_gen_VBFjet2_M;   //!
   TBranch        *b_gen_VBFjet2_Q;   //!
   TBranch        *b_gen_VBFjet2_status;   //!
   TBranch        *b_gen_VBFjet2_Mother;   //!
   TBranch        *b_gen_VBFjet2_GrandMother;   //!
   TBranch        *b_gen_VBFjet1jet2_Pt;   //!
   TBranch        *b_gen_VBFjet1jet2_Eta;   //!
   TBranch        *b_gen_VBFjet1jet2_Phi;   //!
   TBranch        *b_gen_VBFjet1jet2_M;   //!
   TBranch        *b_gen_vbfjet_deltaR;   //!
   TBranch        *b_gen_WHad_Pt;   //!
   TBranch        *b_gen_WHad_M;   //!
   TBranch        *b_gen_WHad_Mt;   //!
   TBranch        *b_gen_WHad_deltaeta;   //!
   TBranch        *b_gen_WHad_deltaphi;   //!
   TBranch        *b_gen_WHad_deltar;   //!
   TBranch        *b_gen_deltaR_LepWHad;   //!
   TBranch        *b_gen_deltaphi_NuWHad;   //!
   TBranch        *b_gen_deltaphi_WlepWHad;   //!
   TBranch        *b_gen_mWW;   //!
   TBranch        *b_gen_mtWW;   //!
   TBranch        *b_gen_mWLep;   //!
   TBranch        *b_gen_mtWLep;   //!
   TBranch        *b_gen_mWHad;   //!
   TBranch        *b_gen_mtWHad;   //!
   TBranch        *b_gen_costheta1;   //!
   TBranch        *b_gen_costheta2;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_costhetastar;   //!
   TBranch        *b_gen_phi1;   //!
   TBranch        *b_gen_dEtajj;   //!
   TBranch        *b_gen_dPhijj;   //!
   TBranch        *b_gen_mjj;   //!
   TBranch        *b_gen_VBSCentrality;   //!
   TBranch        *b_genQuarkStatus;   //!
   TBranch        *b_ngenJet;   //!
   TBranch        *b_nVBFJet;   //!

   Get_Cross_Sec_aQGC(TTree *tree=0);
   virtual ~Get_Cross_Sec_aQGC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outFileName, float Total_Cross_Section);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Get_Cross_Sec_aQGC_cxx
Get_Cross_Sec_aQGC::Get_Cross_Sec_aQGC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/user/rasharma/LHE_GEN_Analyzer_Output/aQGC_WPlepWMhadJJ_EWK_LO_NPle1_mjj200_Pythia8CUEP8M1_13TeV_Madgraph.root");
      }
      f->GetObject("otree",tree);

   }
   Init(tree);
}

Get_Cross_Sec_aQGC::~Get_Cross_Sec_aQGC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Get_Cross_Sec_aQGC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Get_Cross_Sec_aQGC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Get_Cross_Sec_aQGC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LHEWeightIDs = 0;
   LHEWeights = 0;
   genQuarkStatus = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isMuMinus", &isMuMinus, &b_isMuMinus);
   fChain->SetBranchAddress("LHELeptPt", &LHELeptPt, &b_LHELeptPt);
   fChain->SetBranchAddress("LHELeptEta", &LHELeptEta, &b_LHELeptEta);
   fChain->SetBranchAddress("LHELeptPhi", &LHELeptPhi, &b_LHELeptPhi);
   fChain->SetBranchAddress("LHELeptM", &LHELeptM, &b_LHELeptM);
   fChain->SetBranchAddress("LHELeptE", &LHELeptE, &b_LHELeptE);
   fChain->SetBranchAddress("LHENuPt", &LHENuPt, &b_LHENuPt);
   fChain->SetBranchAddress("LHENuEta", &LHENuEta, &b_LHENuEta);
   fChain->SetBranchAddress("LHENuPhi", &LHENuPhi, &b_LHENuPhi);
   fChain->SetBranchAddress("LHENuM", &LHENuM, &b_LHENuM);
   fChain->SetBranchAddress("LHENuE", &LHENuE, &b_LHENuE);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_pt", &LHE_DeltaM_Wqrk0_pt, &b_LHE_DeltaM_Wqrk0_pt);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_eta", &LHE_DeltaM_Wqrk0_eta, &b_LHE_DeltaM_Wqrk0_eta);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_phi", &LHE_DeltaM_Wqrk0_phi, &b_LHE_DeltaM_Wqrk0_phi);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_M", &LHE_DeltaM_Wqrk0_M, &b_LHE_DeltaM_Wqrk0_M);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_E", &LHE_DeltaM_Wqrk0_E, &b_LHE_DeltaM_Wqrk0_E);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk0_Mt", &LHE_DeltaM_Wqrk0_Mt, &b_LHE_DeltaM_Wqrk0_Mt);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_pt", &LHE_DeltaM_Wqrk1_pt, &b_LHE_DeltaM_Wqrk1_pt);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_eta", &LHE_DeltaM_Wqrk1_eta, &b_LHE_DeltaM_Wqrk1_eta);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_phi", &LHE_DeltaM_Wqrk1_phi, &b_LHE_DeltaM_Wqrk1_phi);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_M", &LHE_DeltaM_Wqrk1_M, &b_LHE_DeltaM_Wqrk1_M);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_E", &LHE_DeltaM_Wqrk1_E, &b_LHE_DeltaM_Wqrk1_E);
   fChain->SetBranchAddress("LHE_DeltaM_Wqrk1_Mt", &LHE_DeltaM_Wqrk1_Mt, &b_LHE_DeltaM_Wqrk1_Mt);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_pt", &LHE_DeltaM_Iqrk0_pt, &b_LHE_DeltaM_Iqrk0_pt);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_eta", &LHE_DeltaM_Iqrk0_eta, &b_LHE_DeltaM_Iqrk0_eta);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_phi", &LHE_DeltaM_Iqrk0_phi, &b_LHE_DeltaM_Iqrk0_phi);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_E", &LHE_DeltaM_Iqrk0_E, &b_LHE_DeltaM_Iqrk0_E);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_M", &LHE_DeltaM_Iqrk0_M, &b_LHE_DeltaM_Iqrk0_M);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk0_Mt", &LHE_DeltaM_Iqrk0_Mt, &b_LHE_DeltaM_Iqrk0_Mt);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_pt", &LHE_DeltaM_Iqrk1_pt, &b_LHE_DeltaM_Iqrk1_pt);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_eta", &LHE_DeltaM_Iqrk1_eta, &b_LHE_DeltaM_Iqrk1_eta);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_phi", &LHE_DeltaM_Iqrk1_phi, &b_LHE_DeltaM_Iqrk1_phi);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_E", &LHE_DeltaM_Iqrk1_E, &b_LHE_DeltaM_Iqrk1_E);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_M", &LHE_DeltaM_Iqrk1_M, &b_LHE_DeltaM_Iqrk1_M);
   fChain->SetBranchAddress("LHE_DeltaM_Iqrk1_Mt", &LHE_DeltaM_Iqrk1_Mt, &b_LHE_DeltaM_Iqrk1_Mt);
   fChain->SetBranchAddress("LHE_DeltaM_mWW", &LHE_DeltaM_mWW, &b_LHE_DeltaM_mWW);
   fChain->SetBranchAddress("LHE_DeltaM_mtWW", &LHE_DeltaM_mtWW, &b_LHE_DeltaM_mtWW);
   fChain->SetBranchAddress("LHE_DeltaM_mWLep", &LHE_DeltaM_mWLep, &b_LHE_DeltaM_mWLep);
   fChain->SetBranchAddress("LHE_DeltaM_mtWLep", &LHE_DeltaM_mtWLep, &b_LHE_DeltaM_mtWLep);
   fChain->SetBranchAddress("LHE_DeltaM_mWHad", &LHE_DeltaM_mWHad, &b_LHE_DeltaM_mWHad);
   fChain->SetBranchAddress("LHE_DeltaM_mtWHad", &LHE_DeltaM_mtWHad, &b_LHE_DeltaM_mtWHad);
   fChain->SetBranchAddress("LHE_DeltaM_costheta1", &LHE_DeltaM_costheta1, &b_LHE_DeltaM_costheta1);
   fChain->SetBranchAddress("LHE_DeltaM_costheta2", &LHE_DeltaM_costheta2, &b_LHE_DeltaM_costheta2);
   fChain->SetBranchAddress("LHE_DeltaM_phi", &LHE_DeltaM_phi, &b_LHE_DeltaM_phi);
   fChain->SetBranchAddress("LHE_DeltaM_costhetastar", &LHE_DeltaM_costhetastar, &b_LHE_DeltaM_costhetastar);
   fChain->SetBranchAddress("LHE_DeltaM_phi1", &LHE_DeltaM_phi1, &b_LHE_DeltaM_phi1);
   fChain->SetBranchAddress("LHE_DeltaM_dEtajj", &LHE_DeltaM_dEtajj, &b_LHE_DeltaM_dEtajj);
   fChain->SetBranchAddress("LHE_DeltaM_dPhijj", &LHE_DeltaM_dPhijj, &b_LHE_DeltaM_dPhijj);
   fChain->SetBranchAddress("LHE_DeltaM_mjj", &LHE_DeltaM_mjj, &b_LHE_DeltaM_mjj);
   fChain->SetBranchAddress("LHE_DeltaM_VBSCentrality", &LHE_DeltaM_VBSCentrality, &b_LHE_DeltaM_VBSCentrality);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_pt", &LHE_MothInfo_Wqrk0_pt, &b_LHE_MothInfo_Wqrk0_pt);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_eta", &LHE_MothInfo_Wqrk0_eta, &b_LHE_MothInfo_Wqrk0_eta);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_phi", &LHE_MothInfo_Wqrk0_phi, &b_LHE_MothInfo_Wqrk0_phi);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_M", &LHE_MothInfo_Wqrk0_M, &b_LHE_MothInfo_Wqrk0_M);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_E", &LHE_MothInfo_Wqrk0_E, &b_LHE_MothInfo_Wqrk0_E);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk0_Mt", &LHE_MothInfo_Wqrk0_Mt, &b_LHE_MothInfo_Wqrk0_Mt);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_pt", &LHE_MothInfo_Wqrk1_pt, &b_LHE_MothInfo_Wqrk1_pt);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_eta", &LHE_MothInfo_Wqrk1_eta, &b_LHE_MothInfo_Wqrk1_eta);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_phi", &LHE_MothInfo_Wqrk1_phi, &b_LHE_MothInfo_Wqrk1_phi);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_M", &LHE_MothInfo_Wqrk1_M, &b_LHE_MothInfo_Wqrk1_M);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_E", &LHE_MothInfo_Wqrk1_E, &b_LHE_MothInfo_Wqrk1_E);
   fChain->SetBranchAddress("LHE_MothInfo_Wqrk1_Mt", &LHE_MothInfo_Wqrk1_Mt, &b_LHE_MothInfo_Wqrk1_Mt);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_pt", &LHE_MothInfo_Iqrk0_pt, &b_LHE_MothInfo_Iqrk0_pt);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_eta", &LHE_MothInfo_Iqrk0_eta, &b_LHE_MothInfo_Iqrk0_eta);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_phi", &LHE_MothInfo_Iqrk0_phi, &b_LHE_MothInfo_Iqrk0_phi);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_E", &LHE_MothInfo_Iqrk0_E, &b_LHE_MothInfo_Iqrk0_E);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_M", &LHE_MothInfo_Iqrk0_M, &b_LHE_MothInfo_Iqrk0_M);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk0_Mt", &LHE_MothInfo_Iqrk0_Mt, &b_LHE_MothInfo_Iqrk0_Mt);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_pt", &LHE_MothInfo_Iqrk1_pt, &b_LHE_MothInfo_Iqrk1_pt);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_eta", &LHE_MothInfo_Iqrk1_eta, &b_LHE_MothInfo_Iqrk1_eta);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_phi", &LHE_MothInfo_Iqrk1_phi, &b_LHE_MothInfo_Iqrk1_phi);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_E", &LHE_MothInfo_Iqrk1_E, &b_LHE_MothInfo_Iqrk1_E);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_M", &LHE_MothInfo_Iqrk1_M, &b_LHE_MothInfo_Iqrk1_M);
   fChain->SetBranchAddress("LHE_MothInfo_Iqrk1_Mt", &LHE_MothInfo_Iqrk1_Mt, &b_LHE_MothInfo_Iqrk1_Mt);
   fChain->SetBranchAddress("LHE_MothInfo_mWW", &LHE_MothInfo_mWW, &b_LHE_MothInfo_mWW);
   fChain->SetBranchAddress("LHE_MothInfo_mtWW", &LHE_MothInfo_mtWW, &b_LHE_MothInfo_mtWW);
   fChain->SetBranchAddress("LHE_MothInfo_mWLep", &LHE_MothInfo_mWLep, &b_LHE_MothInfo_mWLep);
   fChain->SetBranchAddress("LHE_MothInfo_mtWLep", &LHE_MothInfo_mtWLep, &b_LHE_MothInfo_mtWLep);
   fChain->SetBranchAddress("LHE_MothInfo_mWHad", &LHE_MothInfo_mWHad, &b_LHE_MothInfo_mWHad);
   fChain->SetBranchAddress("LHE_MothInfo_mtWHad", &LHE_MothInfo_mtWHad, &b_LHE_MothInfo_mtWHad);
   fChain->SetBranchAddress("LHE_MothInfo_costheta1", &LHE_MothInfo_costheta1, &b_LHE_MothInfo_costheta1);
   fChain->SetBranchAddress("LHE_MothInfo_costheta2", &LHE_MothInfo_costheta2, &b_LHE_MothInfo_costheta2);
   fChain->SetBranchAddress("LHE_MothInfo_phi", &LHE_MothInfo_phi, &b_LHE_MothInfo_phi);
   fChain->SetBranchAddress("LHE_MothInfo_costhetastar", &LHE_MothInfo_costhetastar, &b_LHE_MothInfo_costhetastar);
   fChain->SetBranchAddress("LHE_MothInfo_phi1", &LHE_MothInfo_phi1, &b_LHE_MothInfo_phi1);
   fChain->SetBranchAddress("LHE_MothInfo_dEtajj", &LHE_MothInfo_dEtajj, &b_LHE_MothInfo_dEtajj);
   fChain->SetBranchAddress("LHE_MothInfo_dPhijj", &LHE_MothInfo_dPhijj, &b_LHE_MothInfo_dPhijj);
   fChain->SetBranchAddress("LHE_MothInfo_mjj", &LHE_MothInfo_mjj, &b_LHE_MothInfo_mjj);
   fChain->SetBranchAddress("LHE_MothInfo_VBSCentrality", &LHE_MothInfo_VBSCentrality, &b_LHE_MothInfo_VBSCentrality);
   fChain->SetBranchAddress("LHEWeightIDs", &LHEWeightIDs, &b_LHEWeightIDs);
   fChain->SetBranchAddress("LHEWeights", &LHEWeights, &b_LHEWeights);
   fChain->SetBranchAddress("ngen_Lept", &ngen_Lept, &b_ngen_Lept);
   fChain->SetBranchAddress("gen_LeptPt", &gen_LeptPt, &b_gen_LeptPt);
   fChain->SetBranchAddress("gen_LeptEta", &gen_LeptEta, &b_gen_LeptEta);
   fChain->SetBranchAddress("gen_LeptPhi", &gen_LeptPhi, &b_gen_LeptPhi);
   fChain->SetBranchAddress("gen_LeptM", &gen_LeptM, &b_gen_LeptM);
   fChain->SetBranchAddress("gen_LeptStatus", &gen_LeptStatus, &b_gen_LeptStatus);
   fChain->SetBranchAddress("gen_LeptId", &gen_LeptId, &b_gen_LeptId);
   fChain->SetBranchAddress("gen_LeptMother", &gen_LeptMother, &b_gen_LeptMother);
   fChain->SetBranchAddress("gen_LeptGrandMother", &gen_LeptGrandMother, &b_gen_LeptGrandMother);
   fChain->SetBranchAddress("gen_NuPdgId", &gen_NuPdgId, &b_gen_NuPdgId);
   fChain->SetBranchAddress("ngen_Nu", &ngen_Nu, &b_ngen_Nu);
   fChain->SetBranchAddress("gen_NuPt", &gen_NuPt, &b_gen_NuPt);
   fChain->SetBranchAddress("gen_NuEta", &gen_NuEta, &b_gen_NuEta);
   fChain->SetBranchAddress("gen_NuPhi", &gen_NuPhi, &b_gen_NuPhi);
   fChain->SetBranchAddress("gen_NuM", &gen_NuM, &b_gen_NuM);
   fChain->SetBranchAddress("gen_NuQ", &gen_NuQ, &b_gen_NuQ);
   fChain->SetBranchAddress("gen_Nustatus", &gen_Nustatus, &b_gen_Nustatus);
   fChain->SetBranchAddress("gen_NuMother", &gen_NuMother, &b_gen_NuMother);
   fChain->SetBranchAddress("gen_NuGrandMother", &gen_NuGrandMother, &b_gen_NuGrandMother);
   fChain->SetBranchAddress("gen_WJet1_PdgId", &gen_WJet1_PdgId, &b_gen_WJet1_PdgId);
   fChain->SetBranchAddress("ngen_WJet1_", &ngen_WJet1_, &b_ngen_WJet1_);
   fChain->SetBranchAddress("gen_WJet1_Pt", &gen_WJet1_Pt, &b_gen_WJet1_Pt);
   fChain->SetBranchAddress("gen_WJet1_Eta", &gen_WJet1_Eta, &b_gen_WJet1_Eta);
   fChain->SetBranchAddress("gen_WJet1_Phi", &gen_WJet1_Phi, &b_gen_WJet1_Phi);
   fChain->SetBranchAddress("gen_WJet1_M", &gen_WJet1_M, &b_gen_WJet1_M);
   fChain->SetBranchAddress("gen_WJet1_E", &gen_WJet1_E, &b_gen_WJet1_E);
   fChain->SetBranchAddress("gen_WJet1_Q", &gen_WJet1_Q, &b_gen_WJet1_Q);
   fChain->SetBranchAddress("gen_WJet1_status", &gen_WJet1_status, &b_gen_WJet1_status);
   fChain->SetBranchAddress("gen_WJet1_Mother", &gen_WJet1_Mother, &b_gen_WJet1_Mother);
   fChain->SetBranchAddress("gen_WJet1_GrandMother", &gen_WJet1_GrandMother, &b_gen_WJet1_GrandMother);
   fChain->SetBranchAddress("gen_WJet2_PdgId", &gen_WJet2_PdgId, &b_gen_WJet2_PdgId);
   fChain->SetBranchAddress("ngen_WJet2_", &ngen_WJet2_, &b_ngen_WJet2_);
   fChain->SetBranchAddress("gen_WJet2_Pt", &gen_WJet2_Pt, &b_gen_WJet2_Pt);
   fChain->SetBranchAddress("gen_WJet2_Eta", &gen_WJet2_Eta, &b_gen_WJet2_Eta);
   fChain->SetBranchAddress("gen_WJet2_Phi", &gen_WJet2_Phi, &b_gen_WJet2_Phi);
   fChain->SetBranchAddress("gen_WJet2_M", &gen_WJet2_M, &b_gen_WJet2_M);
   fChain->SetBranchAddress("gen_WJet2_E", &gen_WJet2_E, &b_gen_WJet2_E);
   fChain->SetBranchAddress("gen_WJet2_Q", &gen_WJet2_Q, &b_gen_WJet2_Q);
   fChain->SetBranchAddress("gen_WJet2_status", &gen_WJet2_status, &b_gen_WJet2_status);
   fChain->SetBranchAddress("gen_WJet2_Mother", &gen_WJet2_Mother, &b_gen_WJet2_Mother);
   fChain->SetBranchAddress("gen_WJet2_GrandMother", &gen_WJet2_GrandMother, &b_gen_WJet2_GrandMother);
   fChain->SetBranchAddress("gen_VBFjet1_PdgId", &gen_VBFjet1_PdgId, &b_gen_VBFjet1_PdgId);
   fChain->SetBranchAddress("ngen_VBFjet1_", &ngen_VBFjet1_, &b_ngen_VBFjet1_);
   fChain->SetBranchAddress("gen_VBFjet1_Pt", &gen_VBFjet1_Pt, &b_gen_VBFjet1_Pt);
   fChain->SetBranchAddress("gen_VBFjet1_Eta", &gen_VBFjet1_Eta, &b_gen_VBFjet1_Eta);
   fChain->SetBranchAddress("gen_VBFjet1_Phi", &gen_VBFjet1_Phi, &b_gen_VBFjet1_Phi);
   fChain->SetBranchAddress("gen_VBFjet1_M", &gen_VBFjet1_M, &b_gen_VBFjet1_M);
   fChain->SetBranchAddress("gen_VBFjet1_Q", &gen_VBFjet1_Q, &b_gen_VBFjet1_Q);
   fChain->SetBranchAddress("gen_VBFjet1_status", &gen_VBFjet1_status, &b_gen_VBFjet1_status);
   fChain->SetBranchAddress("gen_VBFjet1_Mother", &gen_VBFjet1_Mother, &b_gen_VBFjet1_Mother);
   fChain->SetBranchAddress("gen_VBFjet1_GrandMother", &gen_VBFjet1_GrandMother, &b_gen_VBFjet1_GrandMother);
   fChain->SetBranchAddress("gen_VBFjet2_PdgId", &gen_VBFjet2_PdgId, &b_gen_VBFjet2_PdgId);
   fChain->SetBranchAddress("ngen_VBFjet2_", &ngen_VBFjet2_, &b_ngen_VBFjet2_);
   fChain->SetBranchAddress("gen_VBFjet2_Pt", &gen_VBFjet2_Pt, &b_gen_VBFjet2_Pt);
   fChain->SetBranchAddress("gen_VBFjet2_Eta", &gen_VBFjet2_Eta, &b_gen_VBFjet2_Eta);
   fChain->SetBranchAddress("gen_VBFjet2_Phi", &gen_VBFjet2_Phi, &b_gen_VBFjet2_Phi);
   fChain->SetBranchAddress("gen_VBFjet2_M", &gen_VBFjet2_M, &b_gen_VBFjet2_M);
   fChain->SetBranchAddress("gen_VBFjet2_Q", &gen_VBFjet2_Q, &b_gen_VBFjet2_Q);
   fChain->SetBranchAddress("gen_VBFjet2_status", &gen_VBFjet2_status, &b_gen_VBFjet2_status);
   fChain->SetBranchAddress("gen_VBFjet2_Mother", &gen_VBFjet2_Mother, &b_gen_VBFjet2_Mother);
   fChain->SetBranchAddress("gen_VBFjet2_GrandMother", &gen_VBFjet2_GrandMother, &b_gen_VBFjet2_GrandMother);
   fChain->SetBranchAddress("gen_VBFjet1jet2_Pt", &gen_VBFjet1jet2_Pt, &b_gen_VBFjet1jet2_Pt);
   fChain->SetBranchAddress("gen_VBFjet1jet2_Eta", &gen_VBFjet1jet2_Eta, &b_gen_VBFjet1jet2_Eta);
   fChain->SetBranchAddress("gen_VBFjet1jet2_Phi", &gen_VBFjet1jet2_Phi, &b_gen_VBFjet1jet2_Phi);
   fChain->SetBranchAddress("gen_VBFjet1jet2_M", &gen_VBFjet1jet2_M, &b_gen_VBFjet1jet2_M);
   fChain->SetBranchAddress("gen_vbfjet_deltaR", &gen_vbfjet_deltaR, &b_gen_vbfjet_deltaR);
   fChain->SetBranchAddress("gen_WHad_Pt", &gen_WHad_Pt, &b_gen_WHad_Pt);
   fChain->SetBranchAddress("gen_WHad_M", &gen_WHad_M, &b_gen_WHad_M);
   fChain->SetBranchAddress("gen_WHad_Mt", &gen_WHad_Mt, &b_gen_WHad_Mt);
   fChain->SetBranchAddress("gen_WHad_deltaeta", &gen_WHad_deltaeta, &b_gen_WHad_deltaeta);
   fChain->SetBranchAddress("gen_WHad_deltaphi", &gen_WHad_deltaphi, &b_gen_WHad_deltaphi);
   fChain->SetBranchAddress("gen_WHad_deltar", &gen_WHad_deltar, &b_gen_WHad_deltar);
   fChain->SetBranchAddress("gen_deltaR_LepWHad", &gen_deltaR_LepWHad, &b_gen_deltaR_LepWHad);
   fChain->SetBranchAddress("gen_deltaphi_NuWHad", &gen_deltaphi_NuWHad, &b_gen_deltaphi_NuWHad);
   fChain->SetBranchAddress("gen_deltaphi_WlepWHad", &gen_deltaphi_WlepWHad, &b_gen_deltaphi_WlepWHad);
   fChain->SetBranchAddress("gen_mWW", &gen_mWW, &b_gen_mWW);
   fChain->SetBranchAddress("gen_mtWW", &gen_mtWW, &b_gen_mtWW);
   fChain->SetBranchAddress("gen_mWLep", &gen_mWLep, &b_gen_mWLep);
   fChain->SetBranchAddress("gen_mtWLep", &gen_mtWLep, &b_gen_mtWLep);
   fChain->SetBranchAddress("gen_mWHad", &gen_mWHad, &b_gen_mWHad);
   fChain->SetBranchAddress("gen_mtWHad", &gen_mtWHad, &b_gen_mtWHad);
   fChain->SetBranchAddress("gen_costheta1", &gen_costheta1, &b_gen_costheta1);
   fChain->SetBranchAddress("gen_costheta2", &gen_costheta2, &b_gen_costheta2);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_costhetastar", &gen_costhetastar, &b_gen_costhetastar);
   fChain->SetBranchAddress("gen_phi1", &gen_phi1, &b_gen_phi1);
   fChain->SetBranchAddress("gen_dEtajj", &gen_dEtajj, &b_gen_dEtajj);
   fChain->SetBranchAddress("gen_dPhijj", &gen_dPhijj, &b_gen_dPhijj);
   fChain->SetBranchAddress("gen_mjj", &gen_mjj, &b_gen_mjj);
   fChain->SetBranchAddress("gen_VBSCentrality", &gen_VBSCentrality, &b_gen_VBSCentrality);
   fChain->SetBranchAddress("genQuarkStatus", &genQuarkStatus, &b_genQuarkStatus);
   fChain->SetBranchAddress("ngenJet", &ngenJet, &b_ngenJet);
   fChain->SetBranchAddress("nVBFJet", &nVBFJet, &b_nVBFJet);
   Notify();
}

Bool_t Get_Cross_Sec_aQGC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Get_Cross_Sec_aQGC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Get_Cross_Sec_aQGC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Get_Cross_Sec_aQGC_cxx
