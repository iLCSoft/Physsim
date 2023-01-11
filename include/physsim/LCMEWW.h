#ifndef LCMEWW_H
#define LCMEWW_H
//*****************************************************************************
//* =====================
//*  LCMEWW
//* =====================
//*  
//* (Description)
//*    e+e- -> WW Matrix Element
//*
//* (Update Record)
//*    2012/03/30  K.Fujii	Original version.
//*    2014/10/28  J.Tian       modified for MEM in Marlin
//*****************************************************************************

#include "TNamed.h"
#include "TLorentzVector.h"

#include "physsim/HELLib.h"
#include "physsim/GENLib.h"

#include "physsim/LCMEBase.h"
//_______________________________________________________________________
// =====================
//  class LCMEWW
// =====================
//-----------------------------------------------------------------------
namespace lcme{
  class LCMEWW : public LCMEBase {
  public:
    // --------------------------------------------------------------------
    //  Member Functions
    // --------------------------------------------------------------------
    // ----------------------
    //  C-tor and D-tor
    // ----------------------
    LCMEWW(const char *name, const char *title,
	    Double_t polE = 0.,
	    Double_t polP = 0.,
	    Int_t   iNoDecay = 0);
    virtual ~LCMEWW();
    
    // ----------------------
    //  Getters and Setters
    // ----------------------
    Double_t GetQ2WW       ()           const { return fQ2WW;        }
    Double_t GetCosTheta   ()           const { return fCosTheta;    }
    Double_t GetPhi        ()           const { return fPhi;         }
    Double_t GetQ2Wm       ()           const { return fQ2Wm;        }
    Double_t GetCosThetaWmF()           const { return fCosThetaWmF; }
    Double_t GetPhiWmF     ()           const { return fPhiWmF;      }
    Double_t GetQ2Wp       ()           const { return fQ2Wp;        }
    Double_t GetCosThetaWpF()           const { return fCosThetaWpF; }
    Double_t GetPhiWpF     ()           const { return fPhiWpF;      }
    
    // ----------------------
    //   Base class methods
    // ----------------------
    Double_t GetMatrixElement2();                 // matrix element squared with weighted helicities
    Double_t GetMatrixElement2(Int_t vHel[]);     // matrix element squared with specified helicities
    void     SetMomentumFinal(TLorentzVector vLortz[]); // set four-momenta of final states
    void     SetWDecayMode(Int_t iDecayMode1, Int_t iDecayMode2);             // set W decay mode
    
    // ----------------------
    //   Utility methods
    // ----------------------
  private:
    void Initialize();        // Bases initialization
    Double_t  DSigmaDX     ();
    Complex_t FullAmplitude();
    Complex_t AmpEEtoWW (const HELFermion &em,
			 const HELFermion &ep,
			 const HELVector  &wm,
			 const HELVector  &wp);
    void     SetHelicities(Int_t vHel[]);
    
  private:
    // --------------------------------------------------------------------
    //  Data Members
    // --------------------------------------------------------------------
    // ----------------
    //  Lorentz Vector of final state particle
    TLorentzVector fLortzWm,fLortzWp;     // Wm,Wp
    TLorentzVector fLortzWmf1,fLortzWmf2; // two fermions from Wm decay
    TLorentzVector fLortzWpf1,fLortzWpf2; // two fermions from Wp decay
    
    // ----------------
    //  Particle Data
    // ----------------
    Int_t    fWmDecayMode;           // Wm decay mode;
    Int_t    fWpDecayMode;           // Wp decay mode;
    GENPDTWBoson  *fWmBosonPtr;     //! PD table entry of "Wm"
    GENDecayMode  *fWmModePtr;      // pointer to Wm decay mode
    GENPDTEntry   *f1Ptr;           // point to 1st Wm daughter
    GENPDTEntry   *f2Ptr;           // point to 2nd Wm daughter
    GENPDTWBoson  *fWpBosonPtr;     //! PD table entry of "Wp"
    GENDecayMode  *fWpModePtr;      // pointer to Wp decay mode
    GENPDTEntry   *f3Ptr;           // point to 1st Wp daughter
    GENPDTEntry   *f4Ptr;           // point to 2nd Wp daughter
    GENPDTZBoson  *fZBosonPtr;       //! PD table entry of "Z0"
    
    // ----------------
    //  Event info
    // ----------------
    Int_t          fHelInitial[2];  // initial state helicities
    Int_t          fHelFinal  [4];  // final   state helicities
    ANL4DVector    fK[2];           // [0,1] = [e-, e+]
    ANL4DVector    fP[4];           // [0,1,2,3] = [fwm1,fwm2, fwp1,fwp2]
    Double_t       fM[4];           // [0,1,2,3] = [  m1,  m2,   m3,  m4]
    
    Double_t       fQ2WW;           // q^2 of WW system
    Double_t       fCosTheta;       // cos(theta_x) in cm  frame
    Double_t       fPhi;            // phi_x        in cm  frame
    Double_t       fQ2Wm;          // q^2 of final state Wm
    Double_t       fCosThetaWmF;    // cos(theta_f) in Wm  frame
    Double_t       fPhiWmF;         // phi_f        in Wm  frame
    Double_t       fQ2Wp;          // q^2 of final state Wp
    Double_t       fCosThetaWpF;    // cos(theta_f) in Wp  frame
    Double_t       fPhiWpF;         // phi_f        in Wp  frame

    ClassDef(LCMEWW, 1) // Matrix Element for e+e- -> WW process
  };
}
#endif
