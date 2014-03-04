#ifndef GENLIB_H
#define GENLIB_H
//*****************************************************************************
//* =====================
//*  GENLib
//* =====================
//*  
//* (Description)
//*    Class library for event generators.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "ANL4DVector.h"
#include "TObjArray.h"
#include "TAttLockable.h"

class GENPDTEntry;
class GENPDTWBoson;
class GENPDTZBoson;
class GENPDTPhoton;
class GENPDTGluon;
class GENBranch;

//_______________________________________________________________________
// =====================
//  class GENDecayMode
// =====================
//-----------------------------------------------------------------------
class GENDecayMode : public TObjArray, public lcme::TAttLockable {
friend class GENModePicker;
public:
   GENDecayMode(Double_t gm = 0.) : fGamma(gm), fBR(0.), fCumBR(0.) {}
   virtual ~GENDecayMode() {}

   inline  Double_t GetGamma()             { return fGamma; }
   inline  Double_t GetBR   ()             { return fBR;    }
   inline  void     SetGamma(Double_t gm ) { fGamma = gm;   }
   inline  void     SetBR   (Double_t br ) { fBR    = br;   }

   void     DebugPrint(const Option_t *opt ="");

private:
   Double_t fGamma;		// partial width
   Double_t fBR;		// branching fraction
   Double_t fCumBR;		// cumulative branching fraction

   ClassDef(GENDecayMode, 1) 	// Decay mode class
};

//_______________________________________________________________________
// =====================
//  class GENModePicker
// =====================
//-----------------------------------------------------------------------
class GENModePicker : public TObjArray {
public:
   GENModePicker() : fGamma(0.), fBRsum(0.), fDone(kFALSE) {}
   virtual ~GENModePicker() {}

   using TObjArray::Add;
   virtual void     Add          (GENDecayMode *mp);

           GENDecayMode *PickMode(Double_t x, 
			          Double_t &weight,
				  Int_t    &mode);

   const   GENDecayMode *GetMode (Int_t m) const;
           GENDecayMode *GetMode (Int_t m);
           Double_t      GetWidth() { if(!fDone) Update(); return fGamma; }

protected:
   virtual void     Update();

private:
   Double_t  fGamma;		// total width [GeV]
   Double_t  fBRsum;		// BR sum of unlocked modes
   Bool_t    fDone;		// true if updated

   ClassDef(GENModePicker, 1) 	// Decay mode picker class
};

//_______________________________________________________________________
// =====================
//  class GENPDTEntry
// =====================
//-----------------------------------------------------------------------
class GENPDTEntry: public GENModePicker {
public:
   GENPDTEntry() {}
   GENPDTEntry(const Char_t     *name,
                     Int_t       pid,
                     Double_t    charge,
                     Double_t    spin,
                     Double_t    mass,
                     Int_t       gen   = 0,
                     Double_t    ispin = 0.,
                     Double_t    color = 1.);
   virtual ~GENPDTEntry();

   inline  TString & GetName  () { return fName;    }
   inline  Int_t     GetPID   () { return fPID;     }
   inline  Double_t  GetCharge() { return fCharge;  }
   inline  Double_t  GetMass  () { return fMass;    }
   inline  Int_t     GetGenNo () { return fGen;     }
   inline  Double_t  GetISpin () { return fIsoSpin; }
   inline  Double_t  GetColor () { return fColor;   }

   Double_t  GetQ2BW  (Double_t    qmin,   // Q_min
                       Double_t    qmax,   // Q_max
                       Double_t       x,   // integration variable
                       Double_t &weight);  // Jacobian weight

   void      SetQ2BW  (Double_t    qmin,   // Q_min
                       Double_t    qmax,   // Q_max
                       Double_t      q2,   // Q2
                       Double_t &weight);  // Jacobian weight

   void      DebugPrint(const Option_t *opt = "");

protected:
   TString       fName;		//  name
   Int_t         fPID;		//  PDG ID code
   Double_t      fCharge;	//  charge
   Double_t      fSpin;		//  spin
   Double_t      fMass;		//  mass  [GeV]
   Int_t         fGen;		//  generation
   Double_t      fIsoSpin;	//  (0, 1, 2) = (0, up, down)
   Double_t      fColor;	//  color factor

   ClassDef(GENPDTEntry, 1) 	// PD table entry class
};

//_______________________________________________________________________
// =====================
//  class GENPDTWBoson
// =====================
//-----------------------------------------------------------------------
class GENPDTWBoson: public GENPDTEntry {
public:
   GENPDTWBoson();
   ~GENPDTWBoson();

private:
   void     Initialize();   
   Double_t GamToFF(Double_t mu, Double_t md, Double_t vff, Double_t color);

   ClassDef(GENPDTWBoson, 1)  // W boson class
};

//_______________________________________________________________________
// =====================
//  class GENPDTZBoson
// =====================
//-----------------------------------------------------------------------
class GENPDTZBoson: public GENPDTEntry {
public:
   GENPDTZBoson();
   ~GENPDTZBoson();

private:
   void     Initialize();   
   Double_t GamToFF(Double_t t3, Double_t qf, Double_t cf, Double_t m);

   ClassDef(GENPDTZBoson, 1)  // Z boson class
};

//_______________________________________________________________________
// =====================
//  class GENPDTPhoton
// =====================
//-----------------------------------------------------------------------
class GENPDTPhoton: public GENPDTEntry {
public:
   GENPDTPhoton();
   ~GENPDTPhoton();

private:
   void     Initialize();   

   ClassDef(GENPDTPhoton, 1)  // Photon class
};

//_______________________________________________________________________
// =====================
//  class GENPDTGluon
// =====================
//-----------------------------------------------------------------------
class GENPDTGluon: public GENPDTEntry {
public:
   GENPDTGluon();
   ~GENPDTGluon();

private:
   void     Initialize();   

   ClassDef(GENPDTGluon, 1)  // Gluon class
};

//_______________________________________________________________________
// =====================
//  class GENBranch
// =====================
//-----------------------------------------------------------------------
class GENBranch {
public:
   GENBranch(Double_t q2    = 0.,
             Double_t costh = 0.,
             Double_t phi   = 0.,
             Double_t m12   = 0.,
	     Double_t m22   = 0.);

   GENBranch(Double_t   q2,
             Double_t   costh,
             Double_t   phi,
             GENBranch *br1p,
             GENBranch *br2p);

   GENBranch(Double_t   q2,
             Double_t   costh,
             Double_t   phi,
             GENBranch *br1p,
             Double_t   m22);

   GENBranch(Double_t   q2,
             Double_t   costh,
             Double_t   phi,
             Double_t   m12,
             GENBranch *br2p);

   virtual ~GENBranch() {}

   inline  Double_t GetQ2      ()            { return fQ2;       }
   inline  Double_t GetCosTheta()            { return fCosTheta; }
   inline  Double_t GetPhi     ()            { return fPhi;      }
   inline  Double_t GetM12     ()            { return fM12;      }
   inline  Double_t GetM22     ()            { return fM22;      }
   inline  Double_t GetBetaBar ()            { return fBetaBar;  }

   inline  void     SetQ2      (Double_t q2) { fQ2       = q2;   }
   inline  void     SetCosTheta(Double_t cs) { fCosTheta = cs;   }
   inline  void     SetPhi     (Double_t fi) { fPhi      = fi;   }

   inline  GENBranch * GetBranchPtr(Int_t i) { return i ? fBR2Ptr
                                                       : fBR1Ptr; }

private:
   Double_t   fQ2;        // q^2
   Double_t   fCosTheta;  // cos(theta)
   Double_t   fPhi;       // phi
   Double_t   fM12;       // m1*m1
   Double_t   fM22;       // m2*m2
   GENBranch *fBR1Ptr;    // 1st daughter branch if any
   GENBranch *fBR2Ptr;    // 2nd daughter branch if any
   Double_t   fBetaBar;   // beta_bar

   ClassDef(GENBranch, 1)  // Branch class
};

//_______________________________________________________________________
// =====================
//  class GENFrame
// =====================
//-----------------------------------------------------------------------
class GENFrame {
public:
   GENFrame();
   GENFrame(const ANL4DVector &q, const GENFrame &eb);
   virtual ~GENFrame() {}

   ANL4DVector Transform(const ANL4DVector &pb);

private:
   TVector3 fEV[3];	// axis vectors

   ClassDef(GENFrame, 1)  // Reference frame class
};

//_______________________________________________________________________
// =====================
//  class GENPhase2
// =====================
//-----------------------------------------------------------------------
class GENPhase2 {
public:
   GENPhase2() {}
   GENPhase2(const ANL4DVector &q,
                   Double_t     m12,
                   Double_t     m22,
             const GENFrame    &eb,
                   Double_t     costh,
                   Double_t     phi,
                   Int_t        mode = 0);
   virtual ~GENPhase2() {}

   ANL4DVector GetFourMomentum(Int_t i);
   GENFrame    GetFrame(Int_t i = 1);
   Double_t    GetBetaBar();

private:
   Double_t       Beta2(Double_t x1, Double_t x2);
   void           Update();

private:
   ANL4DVector  fQ;		// parent 4-momentum
   Double_t     fM12;		// 1st dauter mass^2
   Double_t     fM22;		// 2nd dauter mass^2
   GENFrame     fEb;		// Eb: original frame
   GENFrame     fEa;		// Ea: parent's helicity frame
   Double_t     fCosTheta;   	// cos(theta_1) in Ea
   Double_t     fPhi;        	// phi_1        in Ea
   ANL4DVector  fP1;		// 1st daughter 4-momentum in Eb
   ANL4DVector  fP2;		// 2nd daughter 4-momentum in Eb
   Double_t     fBetaBar;	// beta_bar
   Int_t        fMode;		// (0,1)=(no transf, transf)
   Bool_t       fDone;		// true if updated

   ClassDef(GENPhase2, 1)  // 2-body phase space class
};
#endif
