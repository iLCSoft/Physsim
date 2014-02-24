#ifndef ZHHBASES_H
#define ZHHBASES_H
//*****************************************************************************
//* =====================
//*  ZHHBASES
//* =====================
//*  
//* (Description)
//*    e+e- -> ZHH generator
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*    2013/10/03  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "HELLib.h"
#include "GENLib.h"

//_______________________________________________________________________
// =====================
//  class ZHHBases
// =====================
//-----------------------------------------------------------------------
class ZHHBases{
public:
  // --------------------------------------------------------------------
  //  Member Functions
  // --------------------------------------------------------------------
  // ----------------------
  //  C-tor and D-tor
  // ----------------------
  ZHHBases(const char *name  = "ZHHBases", 
           const char *title = "ZHH Bases",
  	   TLorentzVector lortzZf1 = TLorentzVector(0.,0.,0.,0.),
  	   TLorentzVector lortzZf2 = TLorentzVector(0.,0.,0.,0.),
  	   TLorentzVector lortzH1  = TLorentzVector(0.,0.,0.,0.),
  	   TLorentzVector lortzH2  = TLorentzVector(0.,0.,0.,0.)	   );
  virtual ~ZHHBases();

  // ----------------------
  //  Getters and Setters
  // ----------------------
  Double_t GetMass     ()           const { return fMass;      }
  Double_t GetEcmInit  ()           const { return fEcmInit;   }

  Double_t GetQ2ZHH    ()           const { return fQ2ZHH;     }
  Double_t GetCosTheta ()           const { return fCosTheta;  }
  Double_t GetPhi      ()           const { return fPhi;       }
  Double_t GetQ2Z      ()           const { return fQ2Z;       }
  Double_t GetCosThetaF()           const { return fCosThetaF; }
  Double_t GetPhiF     ()           const { return fPhiF;      }
  Double_t GetQ2HH     ()           const { return fQ2HH;      }
  Double_t GetCosThetaH()           const { return fCosThetaH; }
  Double_t GetPhiH     ()           const { return fPhiH;      }
  Double_t GetZBoost   ()           const { return fZBoost;    }
  Double_t GetEcmIP    ()           const { return fEcmIP;     }

  void     SetMass     (Double_t m      ) { fMass      = m;    }
  void     SetEcmInit  (Double_t ecm    ) { fEcmInit   = ecm;  }
  void     SetISR      (Bool_t b = kTRUE) { fISR       = b;    }
  void     SetBeamStr  (Bool_t b = kTRUE) { fBeamStr   = b;    }
  void     SetBeamWidth(Double_t w      ) { fBeamWidth = w;    }
  void     SetPole     (Double_t p      ) { fPole      = p;    }

  // ----------------------
  //   Base class methods
  // ----------------------
  void    SetHelicities(Int_t iHelInit = -1, Int_t iHelFnl = -1); // Helicities for intial and final particles
  virtual void Initialize();        // Bases initialization
  Double_t GetMatrixElement2();     // Bases matrix element squared

  // ----------------------
  //   Utility methods
  // ----------------------
  private:
  void     SelectHelicities(Double_t &weight);

  Double_t  DSigmaDX     (GENBranch &cmbranch);
  Double_t  AmpSquared   (GENBranch &cmbranch);
  Complex_t FullAmplitude();
  Complex_t AmpEEtoZHH   (const HELFermion &em,
                          const HELFermion &ep,
                          const HELScalar  &h1,
                          const HELScalar  &h2,
                          const HELVector  &zf);

private:
  // --------------------------------------------------------------------
  //  Data Members
  // --------------------------------------------------------------------
  // ----------------
  //  Lorentz Vector of final state particle
  TLorentzVector fLortzZf1,fLortzZf2; // two fermions from Z
  TLorentzVector fLortzH1,fLortzH2; // two H

  // ----------------
  // ----------------
  //  Job parameters
  // ----------------
  Double_t fMass;           // m_h    : mass  of H

  Double_t fEcmInit;        // Initial Ecm
  Int_t    fISR;            // ISR on?
  Int_t    fBeamStr;        // Beamstrahlung on?
  Double_t fBeamWidth;      // Beam width relative to Ebm(nominal)
  Double_t fPole;           // electron polarization
  Int_t    fZModesLo;       // Z decay mode lo;
  Int_t    fZModesHi;       // Z decay mode hi;

  // ----------------
  //  Particle Data
  // ----------------
  GENPDTZBoson *fZBosonPtr;       //! PD table entry of "Z"

  // ----------------
  //  Event info
  // ----------------
  Double_t       fZBoost;         // p_z(cm)      in lab frame
  Double_t       fEcmIP;          // Ecm after B-strahlung
  Double_t       fQ2ZHH;          // q^2 of ZHH system
  Double_t       fQ2Z;            // q^2 of final state Z
  Double_t       fQ2HH;           // q^2 of final state HH
  GENDecayMode  *fZModePtr;       // pointer to Z decay mode
  GENPDTEntry   *f3Ptr;           // point to 1st Z daughter
  GENPDTEntry   *f4Ptr;           // point to 2nd Z daughter
  Int_t          fHelInitial[2];  // initial state helicities
  Int_t          fHelFinal  [4];  // final   state helicities
  Int_t          fJCombI;         // intial helicity combination
  Int_t          fJCombF;         // final  helicity combination
  Int_t          fZMode;          // Z decay mode
  ANL4DVector    fK[2];           // [0,1] = [e-, e+]
  ANL4DVector    fP[4];           // [0,1,2,3] = [h1, h2, fz1, fz2]
  Double_t       fM[4];           // [0,1,2,4] = [mh, mh, m3 , m4 ]

  // ----------------
  //  Integ. vars
  // ----------------
  Double_t       fHelCombInitial; // initial state helicity combination
  Double_t       fHelCombFinal;   // final   state helicity combination
  Double_t       fZDecayMode;     // decay mode selector for Z
  Double_t       fCosTheta;       // cos(theta_x) in cm  frame
  Double_t       fPhi;            // phi_x        in cm  frame
  Double_t       fXQ2Z;           // q^2 of final state Z
  Double_t       fCosThetaF;      // cos(theta_H) in HH  frame
  Double_t       fPhiF;           // phi_H        in HH  frame
  Double_t       fXQ2HH;          // q^2 of final state Z
  Double_t       fCosThetaH;      // cos(theta_H) in HH  frame
  Double_t       fPhiH;           // phi_H        in HH  frame

  Double_t       fR_BW_m;         //! naturl Ebeam spread for e-
  Double_t       fR_BW_p;         //!                     for e+
  Double_t       fR_BS_m;         //! beamstrahlung       for e-
  Double_t       fR_BS_p;         //!                     for e+

  Double_t       fR_ISR_var;      //! initial state radiation
  Double_t       fR_ISR_side;     //! from which side ? (e-/e+) 

  ClassDef(ZHHBases, 1) // Bases for e+e- -> XX process
};

#endif
