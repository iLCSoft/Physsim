// *****************************************************
// e+e- ------> eeH
// Example Processor to use Matrix Element Method
//                        ----Junping
// *****************************************************
#include "MEMExampleEEHProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "Utilities.h"

//#include "physsim/LCMEEEH.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace lcme;

MEMExampleEEHProcessor aMEMExampleEEHProcessor ;


MEMExampleEEHProcessor::MEMExampleEEHProcessor() : Processor("MEMExampleEEHProcessor") {
  
  // modify processor description
  _description = "MEMExampleEEHProcessor calculates the differential cross section ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticlesSkimmed") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "InputMCTruthLinkCollection" , 
			   "Name of the MCTruthLink collection"  ,
			   _colMCTL ,
			   std::string("RecoMCTruthLink") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPandoraPFOsCollection" , 
			   "Name of the PandoraPFOs collection"  ,
			   _colPFOs ,
			   std::string("PandoraPFOs") ) ;

}

void MEMExampleEEHProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize LCMEEEH with Higgs mass of 125 GeV and beam polarisations
  // by default P(e-,e+) = (0.,0.) like this case
  _eeh = new LCMEEEH("LCMEEEH","EEH",125.,0.,0.);
  
}

void MEMExampleEEHProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void MEMExampleEEHProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

  cerr << endl << "Hello, Matrix Elemment Example Processor!" << " Event No. " << _nEvt << endl;

  // -- Get the MCTruth Linker --
  //  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  //  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

  // -- Read out MC information --  
  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nMCP = colMC->getNumberOfElements();

  static TNtupleD *hGen = 0;
  if (!hGen) {
    stringstream tupstr_gen;
    tupstr_gen << "sigma:sigmall:sigmarl:sigmalr:sigmarr" << ":"
	       << "mee:meeh:costheta:phi:costhetaf:phif"
	       << ends;
    hGen = new TNtupleD("hGen","",tupstr_gen.str().data());
  }
  TLorentzVector lortzLep1MC,lortzLep2MC,lortzHMC;
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    Double_t energy = mcPart->getEnergy();
    TVector3 pv = TVector3(mcPart->getMomentum());
    TLorentzVector lortz = TLorentzVector(pv,energy);
    // get the lorentz vectors of final states
    if (i == 0) {
      lortzLep1MC += lortz;
    }
    if (i == 1) {
      lortzLep2MC += lortz;
    }
    if (i == 2) {
      lortzHMC += lortz;
    }
  }

  // -----------------------------------
  // calculate the matrix element
  // -----------------------------------
  // put four-momenta of final states to an array
  TLorentzVector vLortzMC[3] = {lortzLep1MC, lortzLep2MC, lortzHMC};
  // pass kinematics to ME object
  _eeh->SetMomentumFinal(vLortzMC);
  // matrix element can be given for each combination of initial and final helicities
  Int_t vHelLLLL[2] = {-2,-2};
  Int_t vHelLRLR[2] = {-1,-1};
  Int_t vHelRLLR[2] = {1,-1};
  Int_t vHelRRRR[2] = {2,2};
  Double_t dSigmaLL = _eeh->GetMatrixElement2(vHelLLLL);
  Double_t dSigmaLR = _eeh->GetMatrixElement2(vHelLRLR);
  Double_t dSigmaRL = _eeh->GetMatrixElement2(vHelRLLR);
  Double_t dSigmaRR = _eeh->GetMatrixElement2(vHelRRRR);
  // if no combination of helicities specified, final combinations are summed 
  // and initial combinations are weighted by beam polarisations
  Double_t dSigma   = _eeh->GetMatrixElement2();
  // that's all need to do to get matrix elment for each event
  // -----------------------------------

  Double_t data_gen[20];
  data_gen[0] = dSigma;
  data_gen[1] = dSigmaLL;
  data_gen[2] = dSigmaRL;
  data_gen[3] = dSigmaLR;
  data_gen[4] = dSigmaRR;
  data_gen[5] = TMath::Sqrt(_eeh->GetQ2EE());
  data_gen[6] = TMath::Sqrt(_eeh->GetQ2EEH());
  data_gen[7] = _eeh->GetCosTheta();
  data_gen[8] = _eeh->GetPhi();
  data_gen[9] = _eeh->GetCosThetaF();
  data_gen[10]= _eeh->GetPhiF();
  hGen->Fill(data_gen);

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  //  _nEvt ++ ;

}



void MEMExampleEEHProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MEMExampleEEHProcessor::end(){ 

  cerr << "MEMExampleEEHProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  
}
