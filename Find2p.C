#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

std::vector<Int_t> GetParticleIndexArray(EventRecord* event, Int_t status, Int_t pdgcode);

//Bool_t ThresholdTest(GHepParticle* particle);
Bool_t ThresholdTest(EventRecord* event);

void Find2p(
    const char* inputFileName = "output_fhc_nue_run1.root",
    const char* outputFileName = "analysis_in_nue_flux_nue_run1.root",
    int startNum = 0,
    int endNum = 1000000)
{
  TFile* myFile = new TFile(inputFileName, "read");
  TTree* myTree = dynamic_cast<TTree *>(myFile->Get("gtree"));

  // initialize event record and connect it to gmcrec (GENIE MC Event Record)
  NtpMCEventRecord* mcrec = new NtpMCEventRecord();
  myTree->SetBranchAddress("gmcrec", &mcrec);

  // Declare and initialize objects and variables
  GHepParticle* particle1 = nullptr;
  GHepParticle* particle2 = nullptr;
  GHepParticle* particle_first_mother_1 = nullptr;
  GHepParticle* particle_last_mother_1 = nullptr;
  GHepParticle* particle_first_mother_2 = nullptr;
  GHepParticle* particle_last_mother_2 = nullptr;
  GHepParticle* p = nullptr;
  GHepParticle* ptemp = nullptr;
  EventRecord* event = nullptr;

  Int_t idx_first_mother_1;
  Int_t idx_last_mother_1;
  Int_t idx_first_mother_2;
  Int_t idx_last_mother_2;
  Int_t pdgid_first_mother_1;
  Int_t pdgid_last_mother_1;
  Int_t pdgid_first_mother_2;
  Int_t pdgid_last_mother_2;
  std::vector<Int_t> idx_target1;
  std::vector<Int_t> idx_target2;

  Float_t weight, prob, xsec, diffxsec, init_nu_e;
  Float_t px1, py1, pz1, e1, vx1, vy1, vz1, vt1;
  Float_t px2, py2, pz2, e2, vx2, vy2, vz2, vt2;
  std::vector<Int_t> pdgcode_etc;
  std::vector<Float_t> px_etc;
  std::vector<Float_t> py_etc;
  std::vector<Float_t> pz_etc;
  std::vector<Float_t> e_etc;
  std::vector<Float_t> vx_etc;
  std::vector<Float_t> vy_etc;
  std::vector<Float_t> vz_etc;
  std::vector<Float_t> vt_etc;

  Bool_t thrCut = false;
  Int_t n_gamma = 0;
  Int_t n_electron = 0;
  Int_t n_positron = 0;
  Int_t n_proton = 0;
  Int_t kPdgTarget1;
  Int_t kPdgTarget2;

  TFile* foutput = new TFile(outputFileName, "recreate");
  TTree* outTree_2g = new TTree("gtree_flat_2g", "gtree_flat_2g");
  outTree_2g->Branch("weight", &weight, "weight/F");
  outTree_2g->Branch("prob", &prob, "prob/F");
  outTree_2g->Branch("xsec", &xsec, "xsec/F");
  outTree_2g->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_2g->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_2g->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_2g->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_2g->Branch("pdgid_first_mother_2", &pdgid_first_mother_2, "pdgid_first_mother_2/I");
  outTree_2g->Branch("pdgid_last_mother_2", &pdgid_last_mother_2, "pdgid_last_mother_2/I");
  outTree_2g->Branch("px1", &px1, "px1/F");
  outTree_2g->Branch("py1", &py1, "py1/F");
  outTree_2g->Branch("pz1", &pz1, "pz1/F");
  outTree_2g->Branch("e1", &e1, "e1/F");
  outTree_2g->Branch("vx1", &vx1, "vx1/F");
  outTree_2g->Branch("vy1", &vy1, "vy1/F");
  outTree_2g->Branch("vz1", &vz1, "vz1/F");
  outTree_2g->Branch("vt1", &vt1, "vt1/F");
  outTree_2g->Branch("px2", &px2, "px2/F");
  outTree_2g->Branch("py2", &py2, "py2/F");
  outTree_2g->Branch("pz2", &pz2, "pz2/F");
  outTree_2g->Branch("e2", &e2, "e2/F");
  outTree_2g->Branch("vx2", &vx2, "vx2/F");
  outTree_2g->Branch("vy2", &vy2, "vy2/F");
  outTree_2g->Branch("vz2", &vz2, "vz2/F");
  outTree_2g->Branch("vt2", &vt2, "vt2/F");

  outTree_2g->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_2g->Branch("px_etc", &px_etc);
  outTree_2g->Branch("py_etc", &py_etc);
  outTree_2g->Branch("pz_etc", &pz_etc);
  outTree_2g->Branch("e_etc", &e_etc);
  outTree_2g->Branch("vx_etc", &vx_etc);
  outTree_2g->Branch("vy_etc", &vy_etc);
  outTree_2g->Branch("vz_etc", &vz_etc);
  outTree_2g->Branch("vt_etc", &vt_etc);

  Long64_t nev = myTree->GetEntries();

  // 2g Loop
  kPdgTarget1 = kPdgGamma;
  kPdgTarget2 = kPdgGamma;

  Int_t n_cnt_2g = 0;
  Int_t n_cnt_thr_2g = 0;
  for(Long64_t i = startNum; i < endNum; ++i)
  {
    if( i % 10000 == 0 )
      std::cout << "[2g]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << std::endl;

		myTree->GetEntry(i);
		event= mcrec->event;

    n_gamma = event->NEntries(kPdgTarget1, kIStStableFinalState);

    // Event Topology Selection
    if( n_gamma != 2 ) continue;
    n_cnt_2g++;

    // Do a threshold test for final state particles
    if( !ThresholdTest(event) ) continue;
    n_cnt_thr_2g++;

    // Find indices of stable target particles
    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();
    particle1 = (GHepParticle*)(*event)[*it_idx_target1];
    particle2 = (GHepParticle*)(*event)[*it_idx_target1+1];

    // Decide leading photon and sub-leading photon
    if( particle1->E() < particle2->E() )
    {
      ptemp = particle1;
      particle1 = particle2;
      particle2 = ptemp;
    }

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();
    idx_first_mother_2 = particle2->FirstMother();
    idx_last_mother_2  = particle2->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];
    particle_first_mother_2 = (GHepParticle*)(*event)[idx_first_mother_2];
    if( idx_last_mother_2 != -1 ) particle_last_mother_2 = (GHepParticle*)(*event)[idx_last_mother_2];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();
    pdgid_first_mother_2 = particle_first_mother_2->Pdg();
    if( idx_last_mother_2 != -1 ) pdgid_last_mother_2 = particle_last_mother_2->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();
    px2 = particle2->Px();
    py2 = particle2->Py();
    pz2 = particle2->Pz();
    e2  = particle2->E();
    vx2 = particle2->Vx();
    vy2 = particle2->Vy();
    vz2 = particle2->Vz();
    vt2 = particle2->Vt();

    outTree_2g->Fill();
  }

  outTree_2g->Write();

  // 1em1ep analysis
  //
  // Initialize variables and objects
  particle1               = nullptr;
  particle2               = nullptr;
  particle_first_mother_1 = nullptr;
  particle_last_mother_1  = nullptr;
  particle_first_mother_2 = nullptr;
  particle_last_mother_2  = nullptr;
  p                       = nullptr;
  ptemp                   = nullptr;
  event                   = nullptr;

  idx_first_mother_1      = -1;
  idx_last_mother_1       = -1;
  idx_first_mother_2      = -1;
  idx_last_mother_2       = -1;
  pdgid_first_mother_1    = -1;
  pdgid_last_mother_1     = -1;
  pdgid_first_mother_2    = -1;
  pdgid_last_mother_2     = -1;
  idx_target1.clear();
  idx_target2.clear();
  weight = 0;
  prob = 0;
  xsec = 0;
  diffxsec = 0;
  init_nu_e = 0;
  px1 = 0;
  py1 = 0;
  pz1 = 0;
  e1  = 0;
  vx1 = 0;
  vy1 = 0;
  vz1 = 0;
  vt1 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  e2  = 0;
  vx2 = 0;
  vy2 = 0;
  vz2 = 0;
  vt2 = 0;
  pdgcode_etc.clear();
  px_etc.clear();
  py_etc.clear();
  pz_etc.clear();
  e_etc.clear();
  vx_etc.clear();
  vy_etc.clear();
  vz_etc.clear();
  vt_etc.clear();

  thrCut = false;
  n_gamma = 0;
  n_electron = 0;
  n_positron = 0;
  n_proton = 0;
  kPdgTarget1 = -999;
  kPdgTarget2 = -999;



  TTree* outTree_1em1ep = new TTree("gtree_flat_1em1ep", "gtree_flat_1em1ep");
  outTree_1em1ep->Branch("weight", &weight, "weight/F");
  outTree_1em1ep->Branch("prob", &prob, "prob/F");
  outTree_1em1ep->Branch("xsec", &xsec, "xsec/F");
  outTree_1em1ep->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1em1ep->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1em1ep->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1em1ep->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1em1ep->Branch("pdgid_first_mother_2", &pdgid_first_mother_2, "pdgid_first_mother_2/I");
  outTree_1em1ep->Branch("pdgid_last_mother_2", &pdgid_last_mother_2, "pdgid_last_mother_2/I");
  outTree_1em1ep->Branch("px1", &px1, "px1/F");
  outTree_1em1ep->Branch("py1", &py1, "py1/F");
  outTree_1em1ep->Branch("pz1", &pz1, "pz1/F");
  outTree_1em1ep->Branch("e1", &e1, "e1/F");
  outTree_1em1ep->Branch("vx1", &vx1, "vx1/F");
  outTree_1em1ep->Branch("vy1", &vy1, "vy1/F");
  outTree_1em1ep->Branch("vz1", &vz1, "vz1/F");
  outTree_1em1ep->Branch("vt1", &vt1, "vt1/F");
  outTree_1em1ep->Branch("px2", &px2, "px2/F");
  outTree_1em1ep->Branch("py2", &py2, "py2/F");
  outTree_1em1ep->Branch("pz2", &pz2, "pz2/F");
  outTree_1em1ep->Branch("e2", &e2, "e2/F");
  outTree_1em1ep->Branch("vx2", &vx2, "vx2/F");
  outTree_1em1ep->Branch("vy2", &vy2, "vy2/F");
  outTree_1em1ep->Branch("vz2", &vz2, "vz2/F");
  outTree_1em1ep->Branch("vt2", &vt2, "vt2/F");

  outTree_1em1ep->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1em1ep->Branch("px_etc", &px_etc);
  outTree_1em1ep->Branch("py_etc", &py_etc);
  outTree_1em1ep->Branch("pz_etc", &pz_etc);
  outTree_1em1ep->Branch("e_etc", &e_etc);
  outTree_1em1ep->Branch("vx_etc", &vx_etc);
  outTree_1em1ep->Branch("vy_etc", &vy_etc);
  outTree_1em1ep->Branch("vz_etc", &vz_etc);
  outTree_1em1ep->Branch("vt_etc", &vt_etc);

  kPdgTarget1 = kPdgElectron;
  kPdgTarget2 = kPdgPositron;

  Int_t n_cnt_1em1ep = 0;
  Int_t n_cnt_thr_1em1ep = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1em1ep]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_electron = event->NEntries(kPdgTarget1, kIStStableFinalState);
    n_positron = event->NEntries(kPdgTarget2, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_electron == 1 && n_positron == 1 ) ) continue;
    n_cnt_1em1ep++;

    // Do a threshold test for final state particles
    if( !ThresholdTest(event) ) continue;
    n_cnt_thr_1em1ep++;

    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);
    idx_target2 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget2);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();
    std::vector<Int_t>::iterator it_idx_target2 = idx_target2.begin();

    particle1 = (GHepParticle*)(*event)[*it_idx_target1];
    particle2 = (GHepParticle*)(*event)[*it_idx_target2];

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();
    idx_first_mother_2 = particle2->FirstMother();
    idx_last_mother_2  = particle2->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];
    particle_first_mother_2 = (GHepParticle*)(*event)[idx_first_mother_2];
    if( idx_last_mother_2 != -1 ) particle_last_mother_2 = (GHepParticle*)(*event)[idx_last_mother_2];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();
    pdgid_first_mother_2 = particle_first_mother_2->Pdg();
    if( idx_last_mother_2 != -1 ) pdgid_last_mother_2 = particle_last_mother_2->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();
    px2 = particle2->Px();
    py2 = particle2->Py();
    pz2 = particle2->Pz();
    e2  = particle2->E();
    vx2 = particle2->Vx();
    vy2 = particle2->Vy();
    vz2 = particle2->Vz();
    vt2 = particle2->Vt();

    outTree_1em1ep->Fill();
  }
  outTree_1em1ep->Write();



  // 1g1em analysis
  //
  // Initialize variables and objects
  particle1               = nullptr;
  particle2               = nullptr;
  particle_first_mother_1 = nullptr;
  particle_last_mother_1  = nullptr;
  particle_first_mother_2 = nullptr;
  particle_last_mother_2  = nullptr;
  p                       = nullptr;
  ptemp                   = nullptr;
  event                   = nullptr;

  idx_first_mother_1      = -1;
  idx_last_mother_1       = -1;
  idx_first_mother_2      = -1;
  idx_last_mother_2       = -1;
  pdgid_first_mother_1    = -1;
  pdgid_last_mother_1     = -1;
  pdgid_first_mother_2    = -1;
  pdgid_last_mother_2     = -1;
  idx_target1.clear();
  idx_target2.clear();
  weight = 0;
  prob = 0;
  xsec = 0;
  diffxsec = 0;
  init_nu_e = 0;
  px1 = 0;
  py1 = 0;
  pz1 = 0;
  e1  = 0;
  vx1 = 0;
  vy1 = 0;
  vz1 = 0;
  vt1 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  e2  = 0;
  vx2 = 0;
  vy2 = 0;
  vz2 = 0;
  vt2 = 0;
  pdgcode_etc.clear();
  px_etc.clear();
  py_etc.clear();
  pz_etc.clear();
  e_etc.clear();
  vx_etc.clear();
  vy_etc.clear();
  vz_etc.clear();
  vt_etc.clear();

  thrCut = false;
  n_gamma = 0;
  n_electron = 0;
  n_positron = 0;
  n_proton = 0;
  kPdgTarget1 = -999;
  kPdgTarget2 = -999;



  TTree* outTree_1g1em = new TTree("gtree_flat_1g1em", "gtree_flat_1g1em");
  outTree_1g1em->Branch("weight", &weight, "weight/F");
  outTree_1g1em->Branch("prob", &prob, "prob/F");
  outTree_1g1em->Branch("xsec", &xsec, "xsec/F");
  outTree_1g1em->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1g1em->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1g1em->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1g1em->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1g1em->Branch("pdgid_first_mother_2", &pdgid_first_mother_2, "pdgid_first_mother_2/I");
  outTree_1g1em->Branch("pdgid_last_mother_2", &pdgid_last_mother_2, "pdgid_last_mother_2/I");
  outTree_1g1em->Branch("px1", &px1, "px1/F");
  outTree_1g1em->Branch("py1", &py1, "py1/F");
  outTree_1g1em->Branch("pz1", &pz1, "pz1/F");
  outTree_1g1em->Branch("e1", &e1, "e1/F");
  outTree_1g1em->Branch("vx1", &vx1, "vx1/F");
  outTree_1g1em->Branch("vy1", &vy1, "vy1/F");
  outTree_1g1em->Branch("vz1", &vz1, "vz1/F");
  outTree_1g1em->Branch("vt1", &vt1, "vt1/F");
  outTree_1g1em->Branch("px2", &px2, "px2/F");
  outTree_1g1em->Branch("py2", &py2, "py2/F");
  outTree_1g1em->Branch("pz2", &pz2, "pz2/F");
  outTree_1g1em->Branch("e2", &e2, "e2/F");
  outTree_1g1em->Branch("vx2", &vx2, "vx2/F");
  outTree_1g1em->Branch("vy2", &vy2, "vy2/F");
  outTree_1g1em->Branch("vz2", &vz2, "vz2/F");
  outTree_1g1em->Branch("vt2", &vt2, "vt2/F");

  outTree_1g1em->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1g1em->Branch("px_etc", &px_etc);
  outTree_1g1em->Branch("py_etc", &py_etc);
  outTree_1g1em->Branch("pz_etc", &pz_etc);
  outTree_1g1em->Branch("e_etc", &e_etc);
  outTree_1g1em->Branch("vx_etc", &vx_etc);
  outTree_1g1em->Branch("vy_etc", &vy_etc);
  outTree_1g1em->Branch("vz_etc", &vz_etc);
  outTree_1g1em->Branch("vt_etc", &vt_etc);

  kPdgTarget1 = kPdgGamma;
  kPdgTarget2 = kPdgElectron;

  Int_t n_cnt_1g1em = 0;
  Int_t n_cnt_thr_1g1em = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1g1em]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_gamma    = event->NEntries(kPdgTarget1, kIStStableFinalState);
    n_electron = event->NEntries(kPdgTarget2, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_gamma == 1 && n_electron == 1 ) ) continue;
    n_cnt_1g1em++;

    // Do a threshold test for final state particles
    if( !ThresholdTest(event) ) continue;
    n_cnt_thr_1g1em++;

    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);
    idx_target2 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget2);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();
    std::vector<Int_t>::iterator it_idx_target2 = idx_target2.begin();

    particle1 = (GHepParticle*)(*event)[*it_idx_target1];
    particle2 = (GHepParticle*)(*event)[*it_idx_target2];

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();
    idx_first_mother_2 = particle2->FirstMother();
    idx_last_mother_2  = particle2->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];
    particle_first_mother_2 = (GHepParticle*)(*event)[idx_first_mother_2];
    if( idx_last_mother_2 != -1 ) particle_last_mother_2 = (GHepParticle*)(*event)[idx_last_mother_2];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();
    pdgid_first_mother_2 = particle_first_mother_2->Pdg();
    if( idx_last_mother_2 != -1 ) pdgid_last_mother_2 = particle_last_mother_2->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();
    px2 = particle2->Px();
    py2 = particle2->Py();
    pz2 = particle2->Pz();
    e2  = particle2->E();
    vx2 = particle2->Vx();
    vy2 = particle2->Vy();
    vz2 = particle2->Vz();
    vt2 = particle2->Vt();

    outTree_1g1em->Fill();
  }
  outTree_1g1em->Write();







  // 1g1ep analysis
  //
  // Initialize variables and objects
  particle1               = nullptr;
  particle2               = nullptr;
  particle_first_mother_1 = nullptr;
  particle_last_mother_1  = nullptr;
  particle_first_mother_2 = nullptr;
  particle_last_mother_2  = nullptr;
  p                       = nullptr;
  ptemp                   = nullptr;
  event                   = nullptr;

  idx_first_mother_1      = -1;
  idx_last_mother_1       = -1;
  idx_first_mother_2      = -1;
  idx_last_mother_2       = -1;
  pdgid_first_mother_1    = -1;
  pdgid_last_mother_1     = -1;
  pdgid_first_mother_2    = -1;
  pdgid_last_mother_2     = -1;
  idx_target1.clear();
  idx_target2.clear();
  weight = 0;
  prob = 0;
  xsec = 0;
  diffxsec = 0;
  init_nu_e = 0;
  px1 = 0;
  py1 = 0;
  pz1 = 0;
  e1  = 0;
  vx1 = 0;
  vy1 = 0;
  vz1 = 0;
  vt1 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  e2  = 0;
  vx2 = 0;
  vy2 = 0;
  vz2 = 0;
  vt2 = 0;
  pdgcode_etc.clear();
  px_etc.clear();
  py_etc.clear();
  pz_etc.clear();
  e_etc.clear();
  vx_etc.clear();
  vy_etc.clear();
  vz_etc.clear();
  vt_etc.clear();

  thrCut = false;
  n_gamma = 0;
  n_electron = 0;
  n_positron = 0;
  n_proton = 0;
  kPdgTarget1 = -999;
  kPdgTarget2 = -999;



  TTree* outTree_1g1ep = new TTree("gtree_flat_1g1ep", "gtree_flat_1g1ep");
  outTree_1g1ep->Branch("weight", &weight, "weight/F");
  outTree_1g1ep->Branch("prob", &prob, "prob/F");
  outTree_1g1ep->Branch("xsec", &xsec, "xsec/F");
  outTree_1g1ep->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1g1ep->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1g1ep->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1g1ep->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1g1ep->Branch("pdgid_first_mother_2", &pdgid_first_mother_2, "pdgid_first_mother_2/I");
  outTree_1g1ep->Branch("pdgid_last_mother_2", &pdgid_last_mother_2, "pdgid_last_mother_2/I");
  outTree_1g1ep->Branch("px1", &px1, "px1/F");
  outTree_1g1ep->Branch("py1", &py1, "py1/F");
  outTree_1g1ep->Branch("pz1", &pz1, "pz1/F");
  outTree_1g1ep->Branch("e1", &e1, "e1/F");
  outTree_1g1ep->Branch("vx1", &vx1, "vx1/F");
  outTree_1g1ep->Branch("vy1", &vy1, "vy1/F");
  outTree_1g1ep->Branch("vz1", &vz1, "vz1/F");
  outTree_1g1ep->Branch("vt1", &vt1, "vt1/F");
  outTree_1g1ep->Branch("px2", &px2, "px2/F");
  outTree_1g1ep->Branch("py2", &py2, "py2/F");
  outTree_1g1ep->Branch("pz2", &pz2, "pz2/F");
  outTree_1g1ep->Branch("e2", &e2, "e2/F");
  outTree_1g1ep->Branch("vx2", &vx2, "vx2/F");
  outTree_1g1ep->Branch("vy2", &vy2, "vy2/F");
  outTree_1g1ep->Branch("vz2", &vz2, "vz2/F");
  outTree_1g1ep->Branch("vt2", &vt2, "vt2/F");

  outTree_1g1ep->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1g1ep->Branch("px_etc", &px_etc);
  outTree_1g1ep->Branch("py_etc", &py_etc);
  outTree_1g1ep->Branch("pz_etc", &pz_etc);
  outTree_1g1ep->Branch("e_etc", &e_etc);
  outTree_1g1ep->Branch("vx_etc", &vx_etc);
  outTree_1g1ep->Branch("vy_etc", &vy_etc);
  outTree_1g1ep->Branch("vz_etc", &vz_etc);
  outTree_1g1ep->Branch("vt_etc", &vt_etc);

  kPdgTarget1 = kPdgGamma;
  kPdgTarget2 = kPdgPositron;

  Int_t n_cnt_1g1ep = 0;
  Int_t n_cnt_thr_1g1ep = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1g1ep]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_gamma    = event->NEntries(kPdgTarget1, kIStStableFinalState);
    n_positron = event->NEntries(kPdgTarget2, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_gamma == 1 && n_positron == 1 ) ) continue;
    n_cnt_1g1ep++;

    // Do a threshold test for final state particles
    if( !ThresholdTest(event) ) continue;
    n_cnt_thr_1g1ep++;

    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);
    idx_target2 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget2);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();
    std::vector<Int_t>::iterator it_idx_target2 = idx_target2.begin();

    particle1 = (GHepParticle*)(*event)[*it_idx_target1];
    particle2 = (GHepParticle*)(*event)[*it_idx_target2];

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();
    idx_first_mother_2 = particle2->FirstMother();
    idx_last_mother_2  = particle2->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];
    particle_first_mother_2 = (GHepParticle*)(*event)[idx_first_mother_2];
    if( idx_last_mother_2 != -1 ) particle_last_mother_2 = (GHepParticle*)(*event)[idx_last_mother_2];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();
    pdgid_first_mother_2 = particle_first_mother_2->Pdg();
    if( idx_last_mother_2 != -1 ) pdgid_last_mother_2 = particle_last_mother_2->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();
    px2 = particle2->Px();
    py2 = particle2->Py();
    pz2 = particle2->Pz();
    e2  = particle2->E();
    vx2 = particle2->Vx();
    vy2 = particle2->Vy();
    vz2 = particle2->Vz();
    vt2 = particle2->Vt();

    outTree_1g1ep->Fill();
  }
  outTree_1g1ep->Write();




  std::cout << "[2g] Events: " << n_cnt_2g << std::endl;
  std::cout << "[2g] Events passed threshold test: " << n_cnt_thr_2g << std::endl;
  std::cout << "[1em1ep] Events: " << n_cnt_1em1ep << std::endl;
  std::cout << "[1em1ep] Events passed threshold test: " << n_cnt_thr_1em1ep << std::endl;
  std::cout << "[1g1ep] Events: " << n_cnt_1g1ep << std::endl;
  std::cout << "[1g1ep] Events passed threshold test: " << n_cnt_thr_1g1ep << std::endl;
  std::cout << "[1g1em] Events: " << n_cnt_1g1em << std::endl;
  std::cout << "[1g1em] Events passed threshold test: " << n_cnt_thr_1g1em << std::endl;


  foutput->Close();
}

std::vector<Int_t> GetParticleIndexArray(EventRecord* event, Int_t status, Int_t pdgcode)
{
  Int_t entries = event->GetEntries();
  GHepParticle* p = nullptr;
  std::vector<Int_t> result;
  for( Int_t i = 0; i < entries; i++ )
  {
    p = (GHepParticle*)(*event)[i];
    if( p->Status() == status && p->Pdg() == pdgcode )
    {
      result.push_back(i);
    }
  }

  return result;
}

Bool_t ThresholdTest(EventRecord* event)
{
  Bool_t ret_value = true;
  Int_t nParticleEntries = event->GetEntries();
  GHepParticle* p = nullptr;

  // Threshold cut values
  const Float_t kEThresholdProton   = 0.05; // unit in GeV
  const Float_t kEThresholdNeutron  = 0.05; // unit in GeV
  const Float_t kEThresholdGamma    = 0.03; // unit in GeV
  const Float_t kEThresholdElectron = 0.03; // unit in GeV
  const Float_t kEThresholdPositron = 0.03; // unit in GeV
  const Float_t kEThresholdMuon     = 0.03; // unit in GeV
  const Float_t kEThresholdPiP      = 0.1;  // unit in GeV
  const Float_t kEThresholdPiM      = 0.1;  // unit in GeV
  const Float_t kEThresholdPi0      = 0.1;  // unit in GeV

  // Loop over for final state particles
  for( Int_t j = 0; j < nParticleEntries; ++j )
  {
    p = (GHepParticle*)(*event)[j];
    if( p->Status() != kIStStableFinalState ) continue;
    switch( p->Pdg() )
    {
      case kPdgProton:
        if( p->E() < kEThresholdProton ) ret_value *= false;
        break;
      case kPdgNeutron:
        if( p->E() < kEThresholdNeutron ) ret_value *= false;
        break;
      case kPdgGamma:
        if( p->E() < kEThresholdGamma ) ret_value *= false;
        break;
      case kPdgElectron:
        if( p->E() < kEThresholdElectron ) ret_value *= false;
        break;
      case kPdgPositron:
        if( p->E() < kEThresholdPositron ) ret_value *= false;
        break;
      case kPdgMuon:
        if( p->E() < kEThresholdMuon ) ret_value *= false;
        break;
      case kPdgPiP:
        if( p->E() < kEThresholdPiP ) ret_value *= false;
        break;
      case kPdgPiM:
        if( p->E() < kEThresholdPiM ) ret_value *= false;
        break;
      case kPdgPi0:
        if( p->E() < kEThresholdPi0 ) ret_value *= false;
        break;
      default:
        ret_value = true;
        break;
    }
  }

  return ret_value;
}
