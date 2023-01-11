#ifndef LCMEZZH_H
#define LCMEZZH_H
//*****************************************************************************
//* =====================
//*  LCMEZZH
//* =====================
//*  
//* (Description)
//*    e+e- -> ZZH Matrix Element
//*
//* (Update Record)
//*    2007/03/29  K.Fujii	Original version.
//*    2014/12/08  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "physsim/HELLib.h"
#include "physsim/GENLib.h"

#include "physsim/LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEZZH
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEZZH : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEZZH(const char *name, const char *title,
	    Double_t massHiggs = 125.,
	    Double_t polE = 0.,
	    Double_t polP = 0.);
    virtual ~LCMEZZH();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetMass       ()           const { return fMass;        }
    Double_t GetCosTheta   ()           const { return fCosTheta;    }
    Double_t GetPhi        ()           const { return fPhi;         }
    Double_t GetQ2ZZ       ()           const { return fQ2ZZ;        }
    Double_t GetCosThetaZ  ()           const { return fCosThetaZ;   }
    Double_t GetPhiZ       ()           const { return fPhiZ;        }
    Double_t GetQ2Z1       ()           const { return fQ2Z1;        }
    Double_t GetCosThetaZ1F()           const { return fCosThetaZ1F; }
    Double_t GetPhiZ1F     ()           const { return fPhiZ1F;      }
    Double_t GetQ2Z2       ()           const { return fQ2Z2;        }
    Double_t GetCosThetaZ2F()           const { return fCosThetaZ2F; }
    Double_t GetPhiZ2F     ()           const { return fPhiZ2F;      }
    Double_t GetQ2ZZH      ()           const { return fQ2ZZH;       }
    Bool_t   GetPropagator ()           const { return fPropagator;  }
    
    // ----------------------
    //   Base class methods
    // ----------------------
    Double_t GetMatrixElement2();                 // matrix element squared with weighted helicities
    Double_t GetMatrixElement2(Int_t vHel[]);     // matrix element squared with specified helicities
    void     SetMass     (Double_t m      ) { fMass      = m;    }
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetZDecayMode(Int_t iDecayMode1, Int_t iDecayMode2);             // set Z decay mode
    void     SetPropagator(Bool_t i) {fPropagator = i;};
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoZZH   (const HELFermion &em,
			    const HELFermion &ep,
			    const HELVector  &z1,
			    const HELVector  &z2,
			    const HELScalar  &hs);
    void     SetHelicities(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzZ1f1,fLortzZ1f2; // two fermions from Z1
    TLorentzVector fLortzZ2f1,fLortzZ2f2; // two fermions from Z1
    TLorentzVector fLortzZ1,fLortzZ2; // Z1, Z2
    TLorentzVector fLortzH; // H
    
    // ----------------
    //  Z decay mode
    // ----------------
    Int_t    fZ1DecayMode,fZ2DecayMode;           // Z decay mode;
    GENDecayMode  *fZ1ModePtr,*fZ2ModePtr;       // point to Z decay mode 
    GENPDTEntry   *f1Ptr;           // point to 1st Z1 daughter
    GENPDTEntry   *f2Ptr;           // point to 2nd Z1 daughter
    GENPDTEntry   *f3Ptr;           // point to 1st Z2 daughter
    GENPDTEntry   *f4Ptr;           // point to 2nd Z2 daughter
    
    // ----------------
    //  Particle Data
    // ----------------
    GENPDTZBoson *fZ1BosonPtr,*fZ2BosonPtr;       //! PD table entry of "Z"
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [5];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[5];           // [0,1,2,3,4] = [h , fz11,fz12, fz21,fz22]
    Double_t       fM[5];           // [0,1,2,3,4] = [mh,   m1,  m2,   m3,  m4]
    
    Double_t       fMass;           // m_h    : mass  of H
    Double_t       fQ2ZZH;          // q^2 of ZZH system
    Double_t       fQ2Z1;           // q^2 of final state Z1
    Double_t       fQ2Z2;           // q^2 of final state Z2
    Double_t       fQ2ZZ;           // q^2 of final state ZZ
    Double_t       fQ2H;            // q^2 of final state H
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fCosThetaZ1F;    // cos(theta_f) in Z1   frame
    Double_t       fPhiZ1F;         // phi_f        in Z1   frame
    Double_t       fCosThetaZ2F;    // cos(theta_f) in Z2   frame
    Double_t       fPhiZ2F;         // phi_f        in Z2   frame
    Double_t       fCosThetaZ;      // cos(theta_Z1) in ZZ  frame
    Double_t       fPhiZ;           // phi_Z1        in ZZ  frame
    
    Bool_t         fPropagator;     // including B-W of Higgs
    
    ClassDef(LCMEZZH, 1) // Matrix Element for e+e- -> ZZH process
  };
}
#endif
