#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

std::vector<Int_t> GetParticleIndexArray(EventRecord* event, Int_t status, Int_t pdgcode);
Bool_t ThresholdTest(EventRecord* event, std::vector<Int_t>* indexArray);
Bool_t Is2gEvent(EventRecord* event, std::vector<Int_t>* idx_target);
Bool_t Is1g1emEvent(EventRecord* event, Int_t* idx_gamma, Int_t* idx_em);
Bool_t Is1g1epEvent(EventRecord* event, Int_t* idx_gamma, Int_t* idx_ep);
Bool_t Is1em1epEvent(EventRecord* event, Int_t* idx_em, Int_t* idx_ep);
Bool_t Is1emEvent(EventRecord* event, Int_t* idx_em);
Bool_t Is1epEvent(EventRecord* event, Int_t* idx_ep);
Bool_t Is1gEvent(EventRecord* event, Int_t* idx_gamma);
Bool_t IsEventHasOtherParticle(EventRecord* event);

// Threshold cut values
const Float_t kEThresholdProton   = 0.05; // unit in GeV
const Float_t kEThresholdNeutron  = 0.05; // unit in GeV
const Float_t kEThresholdGamma    = 0.03; // unit in GeV
const Float_t kEThresholdElectron = 0.03; // unit in GeV
const Float_t kEThresholdPositron = 0.03; // unit in GeV
const Float_t kEThresholdMuon     = 0.03; // unit in GeV
const Float_t kEThresholdPiP      = 0.10; // unit in GeV
const Float_t kEThresholdPiM      = 0.10; // unit in GeV
const Float_t kEThresholdPi0      = 0.10; // unit in GeV

void make_flat_ttree(
    const char* inputFileName = "fhc/nue/output_fhc_nue_run1.root",
    const char* outputFileName = "analysis_in_nue_flux_nue_run1.root",
    int startNum = 0,
    int endNum = 1000000
    )
{
  TFile* myFile = new TFile(inputFileName, "READ");
  TTree* myTree = dynamic_cast<TTree*>(myFile->Get("gtree"));

  // initialize event record and connect it to gmcrec (GENIE MC Event Record)
  NtpMCEventRecord* mcrec = new NtpMCEventRecord();
  myTree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev = myTree->GetEntries();

  // Declare and initialize objects and variables
  GHepParticle* particle1 = nullptr;
  GHepParticle* particle2 = nullptr;
  GHepParticle* ptemp     = nullptr;
  EventRecord* event = nullptr;
  std::vector<Int_t> idx_particle;
  Float_t weight_2g, prob_2g, xsec_2g, diffxsec_2g, init_nu_e_2g;
  Float_t px1_2g, py1_2g, pz1_2g, e1_2g, vx1_2g, vy1_2g, vz1_2g, vt1_2g;
  Float_t px2_2g, py2_2g, pz2_2g, e2_2g, vx2_2g, vy2_2g, vz2_2g, vt2_2g;
  Float_t weight_1g1em, prob_1g1em, xsec_1g1em, diffxsec_1g1em, init_nu_e_1g1em;
  Float_t px1_1g1em, py1_1g1em, pz1_1g1em, e1_1g1em, vx1_1g1em, vy1_1g1em, vz1_1g1em, vt1_1g1em;
  Float_t px2_1g1em, py2_1g1em, pz2_1g1em, e2_1g1em, vx2_1g1em, vy2_1g1em, vz2_1g1em, vt2_1g1em;
  Float_t weight_1g1ep, prob_1g1ep, xsec_1g1ep, diffxsec_1g1ep, init_nu_e_1g1ep;
  Float_t px1_1g1ep, py1_1g1ep, pz1_1g1ep, e1_1g1ep, vx1_1g1ep, vy1_1g1ep, vz1_1g1ep, vt1_1g1ep;
  Float_t px2_1g1ep, py2_1g1ep, pz2_1g1ep, e2_1g1ep, vx2_1g1ep, vy2_1g1ep, vz2_1g1ep, vt2_1g1ep;
  Float_t weight_1em1ep, prob_1em1ep, xsec_1em1ep, diffxsec_1em1ep, init_nu_e_1em1ep;
  Float_t px1_1em1ep, py1_1em1ep, pz1_1em1ep, e1_1em1ep, vx1_1em1ep, vy1_1em1ep, vz1_1em1ep, vt1_1em1ep;
  Float_t px2_1em1ep, py2_1em1ep, pz2_1em1ep, e2_1em1ep, vx2_1em1ep, vy2_1em1ep, vz2_1em1ep, vt2_1em1ep;
  Float_t weight_1em, prob_1em, xsec_1em, diffxsec_1em, init_nu_e_1em;
  Float_t px1_1em, py1_1em, pz1_1em, e1_1em, vx1_1em, vy1_1em, vz1_1em, vt1_1em;
  Float_t px2_1em, py2_1em, pz2_1em, e2_1em, vx2_1em, vy2_1em, vz2_1em, vt2_1em;
  Float_t weight_1ep, prob_1ep, xsec_1ep, diffxsec_1ep, init_nu_e_1ep;
  Float_t px1_1ep, py1_1ep, pz1_1ep, e1_1ep, vx1_1ep, vy1_1ep, vz1_1ep, vt1_1ep;
  Float_t px2_1ep, py2_1ep, pz2_1ep, e2_1ep, vx2_1ep, vy2_1ep, vz2_1ep, vt2_1ep;
  Float_t weight_1g, prob_1g, xsec_1g, diffxsec_1g, init_nu_e_1g;
  Float_t px1_1g, py1_1g, pz1_1g, e1_1g, vx1_1g, vy1_1g, vz1_1g, vt1_1g;
  Float_t px2_1g, py2_1g, pz2_1g, e2_1g, vx2_1g, vy2_1g, vz2_1g, vt2_1g;

  TFile* fOutput = new TFile(outputFileName, "RECREATE");
  TTree* outTree_2g = new TTree("gtree_flat_2g", "gtree_flat_2g");
  TTree* outTree_1g1em = new TTree("gtree_flat_1g1em", "gtree_flat_1g1em");
  TTree* outTree_1g1ep = new TTree("gtree_flat_1g1ep", "gtree_flat_1g1ep");
  TTree* outTree_1em1ep = new TTree("gtree_flat_1em1ep", "gtree_flat_1em1ep");
  TTree* outTree_1em = new TTree("gtree_flat_1em", "gtree_flat_1em");
  TTree* outTree_1ep = new TTree("gtree_flat_1ep", "gtree_flat_1ep");
  TTree* outTree_1g = new TTree("gtree_flat_1g", "gtree_flat_1g");

  outTree_2g->Branch("weight", &weight_2g, "weight/F");
  outTree_2g->Branch("prob", &prob_2g, "prob/F");
  outTree_2g->Branch("xsec", &xsec_2g, "xsec/F");
  outTree_2g->Branch("diffxsec", &diffxsec_2g, "diffxsec/F");
  outTree_2g->Branch("init_nu_e", &init_nu_e_2g, "init_nu_e/F");
  outTree_2g->Branch("px1", &px1_2g, "px1/F");
  outTree_2g->Branch("py1", &py1_2g, "py1/F");
  outTree_2g->Branch("pz1", &pz1_2g, "pz1/F");
  outTree_2g->Branch("e1", &e1_2g, "e1/F");
  outTree_2g->Branch("vx1", &vx1_2g, "vx1/F");
  outTree_2g->Branch("vy1", &vy1_2g, "vy1/F");
  outTree_2g->Branch("vz1", &vz1_2g, "vz1/F");
  outTree_2g->Branch("vt1", &vt1_2g, "vt1/F");
  outTree_2g->Branch("px2", &px2_2g, "px2/F");
  outTree_2g->Branch("py2", &py2_2g, "py2/F");
  outTree_2g->Branch("pz2", &pz2_2g, "pz2/F");
  outTree_2g->Branch("e2", &e2_2g, "e2/F");
  outTree_2g->Branch("vx2", &vx2_2g, "vx2/F");
  outTree_2g->Branch("vy2", &vy2_2g, "vy2/F");
  outTree_2g->Branch("vz2", &vz2_2g, "vz2/F");
  outTree_2g->Branch("vt2", &vt2_2g, "vt2/F");
  outTree_1g1em->Branch("weight", &weight_1g1em, "weight/F");
  outTree_1g1em->Branch("prob", &prob_1g1em, "prob/F");
  outTree_1g1em->Branch("xsec", &xsec_1g1em, "xsec/F");
  outTree_1g1em->Branch("diffxsec", &diffxsec_1g1em, "diffxsec/F");
  outTree_1g1em->Branch("init_nu_e", &init_nu_e_1g1em, "init_nu_e/F");
  outTree_1g1em->Branch("pxg", &px1_1g1em, "pxg/F");
  outTree_1g1em->Branch("pyg", &py1_1g1em, "pyg/F");
  outTree_1g1em->Branch("pzg", &pz1_1g1em, "pzg/F");
  outTree_1g1em->Branch("eg", &e1_1g1em, "eg/F");
  outTree_1g1em->Branch("vxg", &vx1_1g1em, "vxg/F");
  outTree_1g1em->Branch("vyg", &vy1_1g1em, "vyg/F");
  outTree_1g1em->Branch("vzg", &vz1_1g1em, "vzg/F");
  outTree_1g1em->Branch("vtg", &vt1_1g1em, "vtg/F");
  outTree_1g1em->Branch("pxem", &px2_1g1em, "pxem/F");
  outTree_1g1em->Branch("pyem", &py2_1g1em, "pyem/F");
  outTree_1g1em->Branch("pzem", &pz2_1g1em, "pzem/F");
  outTree_1g1em->Branch("eem", &e2_1g1em, "eem/F");
  outTree_1g1em->Branch("vxem", &vx2_1g1em, "vxem/F");
  outTree_1g1em->Branch("vyem", &vy2_1g1em, "vyem/F");
  outTree_1g1em->Branch("vzem", &vz2_1g1em, "vzem/F");
  outTree_1g1em->Branch("vtem", &vt2_1g1em, "vtem/F");
  outTree_1g1ep->Branch("weight", &weight_1g1ep, "weight/F");
  outTree_1g1ep->Branch("prob", &prob_1g1ep, "prob/F");
  outTree_1g1ep->Branch("xsec", &xsec_1g1ep, "xsec/F");
  outTree_1g1ep->Branch("diffxsec", &diffxsec_1g1ep, "diffxsec/F");
  outTree_1g1ep->Branch("init_nu_e", &init_nu_e_1g1ep, "init_nu_e/F");
  outTree_1g1ep->Branch("pxg", &px1_1g1ep, "pxg/F");
  outTree_1g1ep->Branch("pyg", &py1_1g1ep, "pyg/F");
  outTree_1g1ep->Branch("pzg", &pz1_1g1ep, "pzg/F");
  outTree_1g1ep->Branch("eg", &e1_1g1ep, "eg/F");
  outTree_1g1ep->Branch("vxg", &vx1_1g1ep, "vxg/F");
  outTree_1g1ep->Branch("vyg", &vy1_1g1ep, "vyg/F");
  outTree_1g1ep->Branch("vzg", &vz1_1g1ep, "vzg/F");
  outTree_1g1ep->Branch("vtg", &vt1_1g1ep, "vtg/F");
  outTree_1g1ep->Branch("pxep", &px2_1g1ep, "pxep/F");
  outTree_1g1ep->Branch("pyep", &py2_1g1ep, "pyep/F");
  outTree_1g1ep->Branch("pzep", &pz2_1g1ep, "pzep/F");
  outTree_1g1ep->Branch("eep", &e2_1g1ep, "eep/F");
  outTree_1g1ep->Branch("vxep", &vx2_1g1ep, "vxep/F");
  outTree_1g1ep->Branch("vyep", &vy2_1g1ep, "vyep/F");
  outTree_1g1ep->Branch("vzep", &vz2_1g1ep, "vzep/F");
  outTree_1g1ep->Branch("vtep", &vt2_1g1ep, "vtep/F");
  outTree_1em1ep->Branch("weight", &weight_1em1ep, "weight/F");
  outTree_1em1ep->Branch("prob", &prob_1em1ep, "prob/F");
  outTree_1em1ep->Branch("xsec", &xsec_1em1ep, "xsec/F");
  outTree_1em1ep->Branch("diffxsec", &diffxsec_1em1ep, "diffxsec/F");
  outTree_1em1ep->Branch("init_nu_e", &init_nu_e_1em1ep, "init_nu_e/F");
  outTree_1em1ep->Branch("pxem", &px1_1em1ep, "pxem/F");
  outTree_1em1ep->Branch("pyem", &py1_1em1ep, "pyem/F");
  outTree_1em1ep->Branch("pzem", &pz1_1em1ep, "pzem/F");
  outTree_1em1ep->Branch("eem", &e1_1em1ep, "eem/F");
  outTree_1em1ep->Branch("vxem", &vx1_1em1ep, "vxem/F");
  outTree_1em1ep->Branch("vyem", &vy1_1em1ep, "vyem/F");
  outTree_1em1ep->Branch("vzem", &vz1_1em1ep, "vzem/F");
  outTree_1em1ep->Branch("vtem", &vt1_1em1ep, "vtem/F");
  outTree_1em1ep->Branch("pxep", &px2_1em1ep, "pxep/F");
  outTree_1em1ep->Branch("pyep", &py2_1em1ep, "pyep/F");
  outTree_1em1ep->Branch("pzep", &pz2_1em1ep, "pzep/F");
  outTree_1em1ep->Branch("eep", &e2_1em1ep, "eep/F");
  outTree_1em1ep->Branch("vxep", &vx2_1em1ep, "vxep/F");
  outTree_1em1ep->Branch("vyep", &vy2_1em1ep, "vyep/F");
  outTree_1em1ep->Branch("vzep", &vz2_1em1ep, "vzep/F");
  outTree_1em1ep->Branch("vtep", &vt2_1em1ep, "vtep/F");
  outTree_1em->Branch("weight", &weight_1em, "weight/F");
  outTree_1em->Branch("prob", &prob_1em, "prob/F");
  outTree_1em->Branch("xsec", &xsec_1em, "xsec/F");
  outTree_1em->Branch("diffxsec", &diffxsec_1em, "diffxsec/F");
  outTree_1em->Branch("init_nu_e", &init_nu_e_1em, "init_nu_e/F");
  outTree_1em->Branch("pxem", &px1_1em, "pxem/F");
  outTree_1em->Branch("pyem", &py1_1em, "pyem/F");
  outTree_1em->Branch("pzem", &pz1_1em, "pzem/F");
  outTree_1em->Branch("eem", &e1_1em, "eem/F");
  outTree_1em->Branch("vxem", &vx1_1em, "vxem/F");
  outTree_1em->Branch("vyem", &vy1_1em, "vyem/F");
  outTree_1em->Branch("vzem", &vz1_1em, "vzem/F");
  outTree_1em->Branch("vtem", &vt1_1em, "vtem/F");
  outTree_1ep->Branch("weight", &weight_1ep, "weight/F");
  outTree_1ep->Branch("prob", &prob_1ep, "prob/F");
  outTree_1ep->Branch("xsec", &xsec_1ep, "xsec/F");
  outTree_1ep->Branch("diffxsec", &diffxsec_1ep, "diffxsec/F");
  outTree_1ep->Branch("init_nu_e", &init_nu_e_1ep, "init_nu_e/F");
  outTree_1ep->Branch("pxep", &px1_1ep, "pxep/F");
  outTree_1ep->Branch("pyep", &py1_1ep, "pyep/F");
  outTree_1ep->Branch("pzep", &pz1_1ep, "pzep/F");
  outTree_1ep->Branch("eep", &e1_1ep, "eep/F");
  outTree_1ep->Branch("vxep", &vx1_1ep, "vxep/F");
  outTree_1ep->Branch("vyep", &vy1_1ep, "vyep/F");
  outTree_1ep->Branch("vzep", &vz1_1ep, "vzep/F");
  outTree_1ep->Branch("vtep", &vt1_1ep, "vtep/F");
  outTree_1g->Branch("weight", &weight_1g, "weight/F");
  outTree_1g->Branch("prob", &prob_1g, "prob/F");
  outTree_1g->Branch("xsec", &xsec_1g, "xsec/F");
  outTree_1g->Branch("diffxsec", &diffxsec_1g, "diffxsec/F");
  outTree_1g->Branch("init_nu_e", &init_nu_e_1g, "init_nu_e/F");
  outTree_1g->Branch("pxg", &px1_1g, "pxg/F");
  outTree_1g->Branch("pyg", &py1_1g, "pyg/F");
  outTree_1g->Branch("pzg", &pz1_1g, "pzg/F");
  outTree_1g->Branch("eg", &e1_1g, "eg/F");
  outTree_1g->Branch("vxg", &vx1_1g, "vxg/F");
  outTree_1g->Branch("vyg", &vy1_1g, "vyg/F");
  outTree_1g->Branch("vzg", &vz1_1g, "vzg/F");
  outTree_1g->Branch("vtg", &vt1_1g, "vtg/F");

  for(Long64_t i = startNum; i < endNum; ++i)
  {
    idx_particle.clear();
    if( i % 10000 == 0 )
      std::cout << "Processing event " << i << " out of " << nev << "(" << (Double_t)i/(Double_t)nev * 100. << "%)." << std::endl;

    myTree->GetEntry(i);
    event = mcrec->event;

    Int_t idx_gamma, idx_ep, idx_em;
    std::cout << "Event: " << i << ", "
      << Is2gEvent(event, &idx_particle) << " / ";
    for( std::vector<Int_t>::iterator it = idx_particle.begin(); it != idx_particle.end(); ++it) cout << *it << ", ";
    std::cout << "Is1g1emEvent(): " << Is1g1emEvent(event, &idx_gamma, &idx_em) << " / ";
    std::cout << "Is1g1epEvent(): " << Is1g1epEvent(event, &idx_gamma, &idx_ep) << " / ";
    std::cout << "Is1em1epEvent(): " << Is1em1epEvent(event, &idx_em, &idx_ep) << " / ";
    std::cout << "Is1emEvent(): " << Is1emEvent(event, &idx_em) << " / ";
    std::cout << "Is1epEvent(): " << Is1epEvent(event, &idx_ep) << std::endl;
    cout << endl;

    if( Is2gEvent(event, &idx_particle) )
    {
      std::vector<Int_t>::iterator it_idx_particle = idx_particle.begin();
      particle1 = (GHepParticle*)(*event)[*it_idx_particle];
      particle2 = (GHepParticle*)(*event)[*it_idx_particle+1];

      // decide leading photon and sub-leading photon
      if( particle1->E() < particle2->E() )
      {
        ptemp = particle1;
        particle1 = particle2;
        particle2 = ptemp;
      }
      weight_2g    = event->Weight();
      prob_2g      = event->Probability();
      xsec_2g      = event->XSec();
      diffxsec_2g  = event->DiffXSec();
      init_nu_e_2g = ((GHepParticle*)(*event)[0])->E();

      px1_2g = particle1->Px();
      py1_2g = particle1->Py();
      pz1_2g = particle1->Pz();
      e1_2g  = particle1->E();
      vx1_2g = particle1->Vx();
      vy1_2g = particle1->Vy();
      vz1_2g = particle1->Vz();
      vt1_2g = particle1->Vt();
      px2_2g = particle2->Px();
      py2_2g = particle2->Py();
      pz2_2g = particle2->Pz();
      e2_2g  = particle2->E();
      vx2_2g = particle2->Vx();
      vy2_2g = particle2->Vy();
      vz2_2g = particle2->Vz();
      vt2_2g = particle2->Vt();
      outTree_2g->Fill();
    }
    if( Is1g1emEvent(event, &idx_gamma, &idx_em) )
    {
      particle1 = (GHepParticle*)(*event)[idx_gamma];
      particle2 = (GHepParticle*)(*event)[idx_em];

      weight_1g1em    = event->Weight();
      prob_1g1em      = event->Probability();
      xsec_1g1em      = event->XSec();
      diffxsec_1g1em  = event->DiffXSec();
      init_nu_e_1g1em = ((GHepParticle*)(*event)[0])->E();

      px1_1g1em = particle1->Px();
      py1_1g1em = particle1->Py();
      pz1_1g1em = particle1->Pz();
      e1_1g1em  = particle1->E();
      vx1_1g1em = particle1->Vx();
      vy1_1g1em = particle1->Vy();
      vz1_1g1em = particle1->Vz();
      vt1_1g1em = particle1->Vt();
      px2_1g1em = particle2->Px();
      py2_1g1em = particle2->Py();
      pz2_1g1em = particle2->Pz();
      e2_1g1em  = particle2->E();
      vx2_1g1em = particle2->Vx();
      vy2_1g1em = particle2->Vy();
      vz2_1g1em = particle2->Vz();
      vt2_1g1em = particle2->Vt();
      outTree_1g1em->Fill();
    }
    if( Is1g1epEvent(event, &idx_gamma, &idx_ep) )
    {
      particle1 = (GHepParticle*)(*event)[idx_gamma];
      particle2 = (GHepParticle*)(*event)[idx_ep];

      weight_1g1ep    = event->Weight();
      prob_1g1ep      = event->Probability();
      xsec_1g1ep      = event->XSec();
      diffxsec_1g1ep  = event->DiffXSec();
      init_nu_e_1g1ep = ((GHepParticle*)(*event)[0])->E();

      px1_1g1ep = particle1->Px();
      py1_1g1ep = particle1->Py();
      pz1_1g1ep = particle1->Pz();
      e1_1g1ep  = particle1->E();
      vx1_1g1ep = particle1->Vx();
      vy1_1g1ep = particle1->Vy();
      vz1_1g1ep = particle1->Vz();
      vt1_1g1ep = particle1->Vt();
      px2_1g1ep = particle2->Px();
      py2_1g1ep = particle2->Py();
      pz2_1g1ep = particle2->Pz();
      e2_1g1ep  = particle2->E();
      vx2_1g1ep = particle2->Vx();
      vy2_1g1ep = particle2->Vy();
      vz2_1g1ep = particle2->Vz();
      vt2_1g1ep = particle2->Vt();
      outTree_1g1ep->Fill();
    }
    if( Is1em1epEvent(event, &idx_em, &idx_ep) )
    {
      particle1 = (GHepParticle*)(*event)[idx_em];
      particle2 = (GHepParticle*)(*event)[idx_ep];

      weight_1em1ep    = event->Weight();
      prob_1em1ep      = event->Probability();
      xsec_1em1ep      = event->XSec();
      diffxsec_1em1ep  = event->DiffXSec();
      init_nu_e_1em1ep = ((GHepParticle*)(*event)[0])->E();

      px1_1em1ep = particle1->Px();
      py1_1em1ep = particle1->Py();
      pz1_1em1ep = particle1->Pz();
      e1_1em1ep  = particle1->E();
      vx1_1em1ep = particle1->Vx();
      vy1_1em1ep = particle1->Vy();
      vz1_1em1ep = particle1->Vz();
      vt1_1em1ep = particle1->Vt();
      px2_1em1ep = particle2->Px();
      py2_1em1ep = particle2->Py();
      pz2_1em1ep = particle2->Pz();
      e2_1em1ep  = particle2->E();
      vx2_1em1ep = particle2->Vx();
      vy2_1em1ep = particle2->Vy();
      vz2_1em1ep = particle2->Vz();
      vt2_1em1ep = particle2->Vt();
      outTree_1em1ep->Fill();
    }
    if( Is1emEvent(event, &idx_em) )
    {
      particle1 = (GHepParticle*)(*event)[idx_em];

      weight_1em    = event->Weight();
      prob_1em      = event->Probability();
      xsec_1em      = event->XSec();
      diffxsec_1em  = event->DiffXSec();
      init_nu_e_1em = ((GHepParticle*)(*event)[0])->E();

      px1_1em = particle1->Px();
      py1_1em = particle1->Py();
      pz1_1em = particle1->Pz();
      e1_1em  = particle1->E();
      vx1_1em = particle1->Vx();
      vy1_1em = particle1->Vy();
      vz1_1em = particle1->Vz();
      vt1_1em = particle1->Vt();
      outTree_1em->Fill();
    }
    if( Is1epEvent(event, &idx_ep) )
    {
      particle1 = (GHepParticle*)(*event)[idx_em];

      weight_1ep    = event->Weight();
      prob_1ep      = event->Probability();
      xsec_1ep      = event->XSec();
      diffxsec_1ep  = event->DiffXSec();
      init_nu_e_1ep = ((GHepParticle*)(*event)[0])->E();

      px1_1ep = particle1->Px();
      py1_1ep = particle1->Py();
      pz1_1ep = particle1->Pz();
      e1_1ep  = particle1->E();
      vx1_1ep = particle1->Vx();
      vy1_1ep = particle1->Vy();
      vz1_1ep = particle1->Vz();
      vt1_1ep = particle1->Vt();
      outTree_1ep->Fill();
    }
    if( Is1gEvent(event, &idx_gamma) )
    {
      particle1 = (GHepParticle*)(*event)[idx_gamma];

      weight_1g    = event->Weight();
      prob_1g      = event->Probability();
      xsec_1g      = event->XSec();
      diffxsec_1g  = event->DiffXSec();
      init_nu_e_1g = ((GHepParticle*)(*event)[0])->E();

      px1_1g = particle1->Px();
      py1_1g = particle1->Py();
      pz1_1g = particle1->Pz();
      e1_1g  = particle1->E();
      vx1_1g = particle1->Vx();
      vy1_1g = particle1->Vy();
      vz1_1g = particle1->Vz();
      vt1_1g = particle1->Vt();
      outTree_1g->Fill();
    }
  }
  outTree_2g->Write();
  outTree_1g1em->Write();
  outTree_1g1ep->Write();
  outTree_1em->Write();
  outTree_1ep->Write();
  outTree_1g->Write();
}

Bool_t Is2gEvent(EventRecord* event, std::vector<Int_t>* idx_target)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgGamma);
  Int_t n_cnt_2g = idx.size();
  if( idx.size() > 0 )
    for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
    {
      if( *it > idx.size() ) break;
      particle = (GHepParticle*)(*event)[*it];
      if( particle != nullptr )
        if( particle->E() < kEThresholdGamma )
        {
          idx.erase(it);
        }
    }

  result = ( idx.size() == 2 && !IsEventHasOtherParticle(event) ) ? true : false;
  if( result )
    std::copy(idx.begin(), idx.end(), std::back_inserter(*idx_target));

  return result;
}

Bool_t Is1g1emEvent(EventRecord* event, Int_t* idx_gamma, Int_t* idx_em)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgGamma);
  Int_t n_cnt_gamma_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdGamma) )
      idx.erase(it);
  }
  n_cnt_gamma_thr = idx.size();
  if( n_cnt_gamma_thr == 1 ) *idx_gamma = *(idx.begin());

  idx.clear();
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgElectron);
  Int_t n_cnt_electron_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdElectron) )
      idx.erase(it);
  }
  n_cnt_electron_thr = idx.size();
  if( n_cnt_electron_thr == 1 ) *idx_em = *(idx.begin());

  result = (n_cnt_gamma_thr == 1 && n_cnt_electron_thr == 1 && !IsEventHasOtherParticle(event)) ? true : false;

  return result;
}

Bool_t Is1g1epEvent(EventRecord* event, Int_t* idx_gamma, Int_t* idx_ep)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgGamma);
  Int_t n_cnt_gamma_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdGamma) )
      idx.erase(it);
  }
  n_cnt_gamma_thr = idx.size();
  if( n_cnt_gamma_thr == 1 ) *idx_gamma = *(idx.begin());

  idx.clear();
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgPositron);
  Int_t n_cnt_electron_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPositron) )
      idx.erase(it);
  }
  n_cnt_electron_thr = idx.size();
  if( n_cnt_electron_thr == 1 ) *idx_ep = *(idx.begin());

  result = (n_cnt_gamma_thr == 1 && n_cnt_electron_thr == 1 && IsEventHasOtherParticle(event)) ? true : false;
  return result;
}

Bool_t Is1em1epEvent(EventRecord* event, Int_t* idx_em, Int_t* idx_ep)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgElectron);
  Int_t n_cnt_electron_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdElectron) )
      idx.erase(it);
  }
  n_cnt_electron_thr = idx.size();
  if( n_cnt_electron_thr == 1 ) *idx_em = *(idx.begin());

  idx.clear();
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgPositron);
  Int_t n_cnt_positron_thr = 0;
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPositron) )
      idx.erase(it);
  }
  n_cnt_positron_thr = idx.size();
  if( n_cnt_positron_thr == 1 ) *idx_ep = *(idx.begin());

  result = (n_cnt_electron_thr == 1 && n_cnt_positron_thr == 1 && IsEventHasOtherParticle(event)) ? true : false;

  return result;
}

Bool_t Is1emEvent(EventRecord* event, Int_t* idx_em)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgElectron);
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdElectron) )
      idx.erase(it);
  }
  if( idx.size() == 1 ) *idx_em = *(idx.begin());
  result = (idx.size() == 1 && IsEventHasOtherParticle(event)) ? true : false;

  return result;
}

Bool_t Is1epEvent(EventRecord* event, Int_t* idx_ep)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgPositron);
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPositron) )
      idx.erase(it);
  }
  if( idx.size() == 1 ) *idx_ep = *(idx.begin());
  result = (idx.size() == 1 && IsEventHasOtherParticle(event)) ? true : false;

  return result;
}

Bool_t Is1gEvent(EventRecord* event, Int_t* idx_gamma)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  std::vector<Int_t> idx;
  idx = GetParticleIndexArray(event, kIStStableFinalState, kPdgGamma);
  for( std::vector<Int_t>::iterator it = idx.begin(); it != idx.end(); ++it)
  {
    if( *it > idx.size() ) break;
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdGamma) )
      idx.erase(it);
  }
  if( idx.size() == 1 ) *idx_gamma = *(idx.begin());
  result = (idx.size() == 1 && IsEventHasOtherParticle(event)) ? true : false;

  return result;
}

Bool_t IsEventHasOtherParticle(EventRecord* event)
{
  Bool_t result = false;
  GHepParticle* particle = nullptr;
  Int_t n_cnt_proton_thr = 0;
  Int_t n_cnt_neutron_thr = 0;
  Int_t n_cnt_muonp_thr = 0;
  Int_t n_cnt_muonm_thr = 0;
  Int_t n_cnt_pip_thr = 0;
  Int_t n_cnt_pim_thr = 0;
  Int_t n_cnt_pi0_thr = 0;

  // Proton?
  std::vector<Int_t> idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgProton);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdProton) )
      idx_target.erase(it);
  }
  n_cnt_proton_thr = idx_target.size();

  // Neutron?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgNeutron);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdNeutron) )
      idx_target.erase(it);
  }
  n_cnt_neutron_thr = idx_target.size();

  // mu+?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgAntiMuon);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdMuon) )
      idx_target.erase(it);
  }
  n_cnt_muonp_thr = idx_target.size();

  // mu-?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgMuon);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdMuon) )
      idx_target.erase(it);
  }
  n_cnt_muonm_thr = idx_target.size();

  // Pi+?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgPiP);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPiP) )
      idx_target.erase(it);
  }
  n_cnt_pip_thr = idx_target.size();

  // Pi-?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgPiM);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPiM) )
      idx_target.erase(it);
  }
  n_cnt_pim_thr = idx_target.size();

  // Pi0?
  idx_target = GetParticleIndexArray(event, kIStStableFinalState, kPdgPi0);
  for( std::vector<Int_t>::iterator it = idx_target.begin(); it != idx_target.end(); ++it)
  {
    particle = (GHepParticle*)(*event)[*it];
    if( !(particle->E() > kEThresholdPi0) )
      idx_target.erase(it);
  }
  n_cnt_pi0_thr = idx_target.size();

  Int_t ncount = n_cnt_proton_thr + 
    n_cnt_neutron_thr +
    n_cnt_muonm_thr +
    n_cnt_muonp_thr +
    n_cnt_pip_thr +
    n_cnt_pim_thr +
    n_cnt_pi0_thr;

  ncount == 0 ? result = true : result = false;
  return result;
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
