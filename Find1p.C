#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

std::vector<Int_t> GetParticleIndexArray(EventRecord* event, Int_t status, Int_t pdgcode);

//Bool_t ThresholdTest(GHepParticle* particle);
Bool_t ThresholdTest(EventRecord* event, std::vector<Int_t>* indexArray);

void Find1p(
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
  GHepParticle* particle_first_mother_1 = nullptr;
  GHepParticle* particle_last_mother_1 = nullptr;
  GHepParticle* p = nullptr;
  GHepParticle* ptemp = nullptr;
  GHepParticle* p_etc = nullptr;
  EventRecord* event = nullptr;

  Int_t idx_first_mother_1;
  Int_t idx_last_mother_1;
  Int_t pdgid_first_mother_1;
  Int_t pdgid_last_mother_1;
  std::vector<Int_t> idx_target1;

  Int_t interactionType, scatteringType;
  Float_t weight, prob, xsec, diffxsec, init_nu_e;
  Float_t px1, py1, pz1, e1, vx1, vy1, vz1, vt1;
  std::vector<Int_t> idx_etc;
  std::vector<Int_t>::iterator it_idx_etc;
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

  // 1g analysis


  TFile* foutput = new TFile(outputFileName, "recreate");
  TTree* outTree_1g = new TTree("gtree_flat_1g", "gtree_flat_1g");
  outTree_1g->Branch("interactionType", &interactionType, "interactionType/I");
  outTree_1g->Branch("scatteringType", &scatteringType, "scatteringType/I");
  outTree_1g->Branch("weight", &weight, "weight/F");
  outTree_1g->Branch("prob", &prob, "prob/F");
  outTree_1g->Branch("xsec", &xsec, "xsec/F");
  outTree_1g->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1g->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1g->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1g->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1g->Branch("px1", &px1, "px1/F");
  outTree_1g->Branch("py1", &py1, "py1/F");
  outTree_1g->Branch("pz1", &pz1, "pz1/F");
  outTree_1g->Branch("e1", &e1, "e1/F");
  outTree_1g->Branch("vx1", &vx1, "vx1/F");
  outTree_1g->Branch("vy1", &vy1, "vy1/F");
  outTree_1g->Branch("vz1", &vz1, "vz1/F");
  outTree_1g->Branch("vt1", &vt1, "vt1/F");

  outTree_1g->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1g->Branch("px_etc", &px_etc);
  outTree_1g->Branch("py_etc", &py_etc);
  outTree_1g->Branch("pz_etc", &pz_etc);
  outTree_1g->Branch("e_etc", &e_etc);
  outTree_1g->Branch("vx_etc", &vx_etc);
  outTree_1g->Branch("vy_etc", &vy_etc);
  outTree_1g->Branch("vz_etc", &vz_etc);
  outTree_1g->Branch("vt_etc", &vt_etc);

  Long64_t nev = myTree->GetEntries();

  kPdgTarget1 = kPdgGamma;

  Int_t n_cnt_1g = 0;
  Int_t n_cnt_thr_1g = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1g]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_gamma = event->NEntries(kPdgTarget1, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_gamma == 1 ) ) continue;
    n_cnt_1g++;

    // Do a threshold test for final state particles
    idx_etc.clear();
    thrCut = ThresholdTest(event, &idx_etc);
    if( !thrCut ) continue;
    n_cnt_thr_1g++;

    // Find indices of stable target particles
    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);
    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();
    particle1 = (GHepParticle*)(*event)[*it_idx_target1];

    // Fill the variables
    interactionType = event->Summary()->ProcInfo().InteractionTypeId();
    scatteringType  = event->Summary()->ProcInfo().ScatteringTypeId();
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();

    if( thrCut )
    {
      std::cout << "Event: " << i << " passed threshold test. Attemping to add etc particles." << std::endl;
      pdgcode_etc.clear();
      px_etc.clear();
      py_etc.clear();
      pz_etc.clear();
      e_etc.clear();
      vx_etc.clear();
      vy_etc.clear();
      vz_etc.clear();
      vt_etc.clear();
      for( it_idx_etc = idx_etc.begin(); it_idx_etc != idx_etc.end(); ++it_idx_etc)
      {
        p_etc = (GHepParticle*)(*event)[*it_idx_etc];
        pdgcode_etc.push_back(p_etc->Pdg());
        px_etc.push_back(p_etc->Px());
        py_etc.push_back(p_etc->Py());
        pz_etc.push_back(p_etc->Pz());
        e_etc.push_back(p_etc->E());
        vx_etc.push_back(p_etc->Vx());
        vy_etc.push_back(p_etc->Vy());
        vz_etc.push_back(p_etc->Vz());
        vt_etc.push_back(p_etc->Vt());
      }
    }

    outTree_1g->Fill();
  }
  outTree_1g->Write();







  // 1ep analysis
  //
  // Initialize variables and objects
  particle1               = nullptr;
  particle_first_mother_1 = nullptr;
  particle_last_mother_1  = nullptr;
  p                       = nullptr;
  p_etc                   = nullptr;
  ptemp                   = nullptr;
  event                   = nullptr;

  idx_first_mother_1      = -1;
  idx_last_mother_1       = -1;
  pdgid_first_mother_1    = -1;
  pdgid_last_mother_1     = -1;
  idx_target1.clear();
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



  TTree* outTree_1ep = new TTree("gtree_flat_1ep", "gtree_flat_1ep");
  outTree_1ep->Branch("weight", &weight, "weight/F");
  outTree_1ep->Branch("prob", &prob, "prob/F");
  outTree_1ep->Branch("xsec", &xsec, "xsec/F");
  outTree_1ep->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1ep->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1ep->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1ep->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1ep->Branch("px1", &px1, "px1/F");
  outTree_1ep->Branch("py1", &py1, "py1/F");
  outTree_1ep->Branch("pz1", &pz1, "pz1/F");
  outTree_1ep->Branch("e1", &e1, "e1/F");
  outTree_1ep->Branch("vx1", &vx1, "vx1/F");
  outTree_1ep->Branch("vy1", &vy1, "vy1/F");
  outTree_1ep->Branch("vz1", &vz1, "vz1/F");
  outTree_1ep->Branch("vt1", &vt1, "vt1/F");

  outTree_1ep->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1ep->Branch("px_etc", &px_etc);
  outTree_1ep->Branch("py_etc", &py_etc);
  outTree_1ep->Branch("pz_etc", &pz_etc);
  outTree_1ep->Branch("e_etc", &e_etc);
  outTree_1ep->Branch("vx_etc", &vx_etc);
  outTree_1ep->Branch("vy_etc", &vy_etc);
  outTree_1ep->Branch("vz_etc", &vz_etc);
  outTree_1ep->Branch("vt_etc", &vt_etc);

  kPdgTarget1 = kPdgPositron;

  Int_t n_cnt_1ep = 0;
  Int_t n_cnt_thr_1ep = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1ep]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_positron = event->NEntries(kPdgTarget1, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_positron == 1 ) ) continue;
    n_cnt_1ep++;

    // Do a threshold test for final state particles
    idx_etc.clear();
    thrCut = ThresholdTest(event, &idx_etc);
    if( !thrCut ) continue;
    n_cnt_thr_1ep++;

    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();

    particle1 = (GHepParticle*)(*event)[*it_idx_target1];

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();

    if( thrCut )
    {
      std::cout << "Event: " << i << " passed threshold test. Attemping to add etc particles." << std::endl;
      pdgcode_etc.clear();
      px_etc.clear();
      py_etc.clear();
      pz_etc.clear();
      e_etc.clear();
      vx_etc.clear();
      vy_etc.clear();
      vz_etc.clear();
      vt_etc.clear();
      for( it_idx_etc = idx_etc.begin(); it_idx_etc != idx_etc.end(); ++it_idx_etc)
      {
        p_etc = (GHepParticle*)(*event)[*it_idx_etc];
        pdgcode_etc.push_back(p_etc->Pdg());
        px_etc.push_back(p_etc->Px());
        py_etc.push_back(p_etc->Py());
        pz_etc.push_back(p_etc->Pz());
        e_etc.push_back(p_etc->E());
        vx_etc.push_back(p_etc->Vx());
        vy_etc.push_back(p_etc->Vy());
        vz_etc.push_back(p_etc->Vz());
        vt_etc.push_back(p_etc->Vt());
      }
    }

    outTree_1ep->Fill();
  }
  outTree_1ep->Write();






  // 1em analysis
  //
  // Initialize variables and objects
  particle1               = nullptr;
  particle_first_mother_1 = nullptr;
  particle_last_mother_1  = nullptr;
  p                       = nullptr;
  p_etc                   = nullptr;
  ptemp                   = nullptr;
  event                   = nullptr;

  idx_first_mother_1      = -1;
  idx_last_mother_1       = -1;
  pdgid_first_mother_1    = -1;
  pdgid_last_mother_1     = -1;
  idx_target1.clear();
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



  TTree* outTree_1em = new TTree("gtree_flat_1em", "gtree_flat_1em");
  outTree_1em->Branch("weight", &weight, "weight/F");
  outTree_1em->Branch("prob", &prob, "prob/F");
  outTree_1em->Branch("xsec", &xsec, "xsec/F");
  outTree_1em->Branch("diffxsec", &diffxsec, "diffxsec/F");
  outTree_1em->Branch("init_nu_e", &init_nu_e, "init_nu_e/F");
  outTree_1em->Branch("pdgid_first_mother_1", &pdgid_first_mother_1, "pdgid_first_mother_1/I");
  outTree_1em->Branch("pdgid_last_mother_1", &pdgid_last_mother_1, "pdgid_last_mother_1/I");
  outTree_1em->Branch("px1", &px1, "px1/F");
  outTree_1em->Branch("py1", &py1, "py1/F");
  outTree_1em->Branch("pz1", &pz1, "pz1/F");
  outTree_1em->Branch("e1", &e1, "e1/F");
  outTree_1em->Branch("vx1", &vx1, "vx1/F");
  outTree_1em->Branch("vy1", &vy1, "vy1/F");
  outTree_1em->Branch("vz1", &vz1, "vz1/F");
  outTree_1em->Branch("vt1", &vt1, "vt1/F");

  outTree_1em->Branch("pdgcode_etc", &pdgcode_etc);
  outTree_1em->Branch("px_etc", &px_etc);
  outTree_1em->Branch("py_etc", &py_etc);
  outTree_1em->Branch("pz_etc", &pz_etc);
  outTree_1em->Branch("e_etc", &e_etc);
  outTree_1em->Branch("vx_etc", &vx_etc);
  outTree_1em->Branch("vy_etc", &vy_etc);
  outTree_1em->Branch("vz_etc", &vz_etc);
  outTree_1em->Branch("vt_etc", &vt_etc);

  kPdgTarget1 = kPdgElectron;

  Int_t n_cnt_1em = 0;
  Int_t n_cnt_thr_1em = 0;
  for(Long64_t i = startNum; i < endNum; i++)
  {
    if( i % 10000 == 0 ) cout << "[1em]: Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << endl;
		myTree->GetEntry(i);
		event= mcrec->event;

    n_electron = event->NEntries(kPdgTarget1, kIStStableFinalState);

    // ================== Event Topology Selection ===========================
    if( !( n_electron == 1 ) ) continue;
    n_cnt_1em++;

    // Do a threshold test for final state particles
    idx_etc.clear();
    thrCut = ThresholdTest(event, &idx_etc);
    if( !thrCut ) continue;
    n_cnt_thr_1em++;

    idx_target1 = GetParticleIndexArray(event, kIStStableFinalState, kPdgTarget1);

    std::vector<Int_t>::iterator it_idx_target1 = idx_target1.begin();

    particle1 = (GHepParticle*)(*event)[*it_idx_target1];

    // Fill the variables
    weight    = event->Weight();
    prob      = event->Probability();
    xsec      = event->XSec();
    diffxsec  = event->DiffXSec();
    init_nu_e = ((GHepParticle*)(*event)[0])->E();

    idx_first_mother_1 = particle1->FirstMother();
    idx_last_mother_1  = particle1->LastMother();

    particle_first_mother_1 = (GHepParticle*)(*event)[idx_first_mother_1];
    if( idx_last_mother_1 != -1 ) particle_last_mother_1 = (GHepParticle*)(*event)[idx_last_mother_1];

    pdgid_first_mother_1 = particle_first_mother_1->Pdg();
    if( idx_last_mother_1 != -1 ) pdgid_last_mother_1 = particle_last_mother_1->Pdg();

    px1 = particle1->Px();
    py1 = particle1->Py();
    pz1 = particle1->Pz();
    e1  = particle1->E();
    vx1 = particle1->Vx();
    vy1 = particle1->Vy();
    vz1 = particle1->Vz();
    vt1 = particle1->Vt();

    if( thrCut )
    {
      std::cout << "Event: " << i << " passed threshold test. Attemping to add etc particles." << std::endl;
      pdgcode_etc.clear();
      px_etc.clear();
      py_etc.clear();
      pz_etc.clear();
      e_etc.clear();
      vx_etc.clear();
      vy_etc.clear();
      vz_etc.clear();
      vt_etc.clear();
      for( it_idx_etc = idx_etc.begin(); it_idx_etc != idx_etc.end(); ++it_idx_etc)
      {
        p_etc = (GHepParticle*)(*event)[*it_idx_etc];
        pdgcode_etc.push_back(p_etc->Pdg());
        px_etc.push_back(p_etc->Px());
        py_etc.push_back(p_etc->Py());
        pz_etc.push_back(p_etc->Pz());
        e_etc.push_back(p_etc->E());
        vx_etc.push_back(p_etc->Vx());
        vy_etc.push_back(p_etc->Vy());
        vz_etc.push_back(p_etc->Vz());
        vt_etc.push_back(p_etc->Vt());
      }
    }

    outTree_1em->Fill();
  }
  outTree_1em->Write();




  std::cout << "[1g] Events: " << n_cnt_1g << std::endl;
  std::cout << "[1g] Events passed threshold test: " << n_cnt_thr_1g << std::endl;
  std::cout << "[1ep] Events: " << n_cnt_1ep << std::endl;
  std::cout << "[1ep] Events passed threshold test: " << n_cnt_thr_1ep << std::endl;
  std::cout << "[1em] Events: " << n_cnt_1em << std::endl;
  std::cout << "[1em] Events passed threshold test: " << n_cnt_thr_1em << std::endl;


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

Bool_t ThresholdTest(EventRecord* event, std::vector<Int_t>* indexArray)
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
        else indexArray->push_back(j);
        break;
      case kPdgNeutron:
        if( p->E() < kEThresholdNeutron ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgGamma:
        if( p->E() < kEThresholdGamma ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgElectron:
        if( p->E() < kEThresholdElectron ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgPositron:
        if( p->E() < kEThresholdPositron ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgMuon:
        if( p->E() < kEThresholdMuon ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgPiP:
        if( p->E() < kEThresholdPiP ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgPiM:
        if( p->E() < kEThresholdPiM ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      case kPdgPi0:
        if( p->E() < kEThresholdPi0 ) ret_value *= false;
        else indexArray->push_back(j);
        break;
      default:
        ret_value = true;
        break;
    }
  }

  return ret_value;
}
