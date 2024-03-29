#ifndef HELLIB_H
#define HELLIB_H
//*****************************************************************************
//* =====================
//*  HELLib Classes
//* =====================
//*  
//* (Description)
//*    Class library for helicity amplitude calculations.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version.
//*
//*****************************************************************************

#include "physsim/ANL4DVector.h"

#include <complex>
typedef std::complex<Double_t>  Complex_t;

class HELFermion;
class HELScalar;
class HELVector;

//_______________________________________________________________________
// =====================
//  class TVectorC
// =====================
//-----------------------------------------------------------------------
class TVectorC : public TObject {
public:
#if 0
  TVectorC(Int_t n = 4) : fData(n) {}
#else
  TVectorC(Int_t n = 4) { fData[0]=0.; fData[1]=0.; fData[2]=0.; fData[3]=0.; }
#endif
  //  TVectorC(const TVectorC &src) : fData(src.fData) {}
  // fg: the above does not compile, e.g. on clang
  TVectorC(const TVectorC &src) { fData[0]=src[0]; fData[1]=src[1]; fData[2]=src[2]; fData[3]=src[3]; }

  Complex_t  operator[](Int_t i) const { return fData[i]; }
  Complex_t &operator[](Int_t i)       { return fData[i]; }

private:
#if 0
  std::vector<Complex_t> fData; // data
#else
  Complex_t   fData[4]; // data
#endif

  ClassDef(TVectorC, 1) 	// Complex vector class
};

//_______________________________________________________________________
// =====================
//  class HELFermion
// =====================
//-----------------------------------------------------------------------
class HELFermion : public TVectorC {
friend class HELVector;
friend class HELScalar;
friend class HELVertex;
public:
#if 0
   HELFermion() : TVectorC(4) {}
#else
   HELFermion() {}
#endif
   HELFermion(const ANL4DVector &p,
                    Double_t     m,
                    Int_t        hel,
                    Int_t        nsf = 1,
		    Bool_t       isincom = kFALSE);
   HELFermion(const HELFermion  &f,
              const HELVector   &v,
                    Double_t     gl,             
                    Double_t     gr,             
                    Double_t     m,
                    Double_t     gam);
   virtual ~HELFermion() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Double_t     GetMass        () const  { return fM;   }
   inline       Int_t        GetHelicity    () const  { return fHel; }
   inline       Int_t        GetNSF         () const  { return fNSF; }
   
private:
   ANL4DVector fP;                // 4-momentum * fNSF
   Double_t    fM;                // mass
   Int_t       fHel;              // (-1,+1) = (-1/2,+1/2)
   Int_t       fNSF;              // (-1,+1) = (antiparticle, particle)
   Bool_t      fIsIncoming;       // true if incoming

   ClassDef(HELFermion, 1)  // Incoming fermion class
};

//_______________________________________________________________________
// =====================
//  class HELVector
// =====================
//-----------------------------------------------------------------------
class HELVector: public TVectorC {
friend class HELFermion;
friend class HELScalar;
friend class HELVertex;
public:
#if 0
   HELVector() : TVectorC(4) {}
#else
   HELVector() {}
#endif
   HELVector(const ANL4DVector &p,
                   Double_t     m,
	           Int_t        hel,
		   Int_t        nsv = 1);
   HELVector(const HELFermion &fin, 
             const HELFermion &fout,
	           Double_t    gl,
		   Double_t    gr,
		   Double_t    m,
		   Double_t    gm);
   HELVector(const HELFermion &fin, 
             const HELFermion &fout,
	           Double_t    gla,
		   Double_t    gra,
	           Double_t    glz,
		   Double_t    grz,
		   Double_t    m,
		   Double_t    gm);
   HELVector(const HELVector   &v,
             const HELScalar   &sc,
		   Double_t    g,
		   Double_t    m,
		   Double_t    gm);
   HELVector(const HELVector   &v1,
             const HELVector   &v2,
		   Double_t    g,
		   Double_t    m,
		   Double_t    gm);
   HELVector (Double_t ebm, 
              Double_t eef,
	      Double_t sh,
	      Double_t ch,
	      Double_t fi,
	      Int_t    helbm,
	      Int_t    helef,
	      Int_t    nsf,
	      Double_t ge = TMath::Sqrt(4.*TMath::Pi()/128.),
	      Double_t me = 0.510998902e-3);
   HELVector (Complex_t v0,
              Complex_t v1,
              Complex_t v2,
              Complex_t v3,
	      const ANL4DVector &p);

   virtual ~HELVector() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Double_t     GetMass        () const  { return fM;   }
   inline       Int_t        GetHelicity    () const  { return fHel; }
   inline       Int_t        GetNSV         () const  { return fNSV; }
#if 1
   void DebugPrint() const;
#endif
   
private:
   ANL4DVector fP;                // 4-momentum * fNSV
   Double_t    fM;                // mass
   Int_t       fHel;              // helicity
   Int_t       fNSV;              // (-1,1) = (initial, final)

   ClassDef(HELVector, 1)  // Vector boson class
};

//_______________________________________________________________________
// =====================
//  class HELScalar
// =====================
//-----------------------------------------------------------------------
class HELScalar: public Complex_t {
friend class HELFermion;
friend class HELVector;
friend class HELVertex;
public:
   HELScalar(const ANL4DVector &p,
		   Int_t        nss = 1);
   HELScalar(const HELVector  &vc,
             const HELScalar  &sc,
                   Double_t    g,
                   Double_t    m,
                   Double_t    gm);
   HELScalar(const HELScalar  &s1,
             const HELScalar  &s2,
                   Double_t    g,
                   Double_t    m,
                   Double_t    gm);
   HELScalar(const HELVector  &v1,
             const HELVector  &v2,
                   Double_t    g,
                   Double_t    m,
                   Double_t    gm);
   HELScalar(const HELVector  &v1,
             const HELVector  &v2,
                   Double_t    g1,
                   Double_t    g2,
                   Double_t    g3,
                   Double_t    m,
                   Double_t    gm);
   HELScalar(const HELFermion &in,
             const HELFermion &out,
                   Complex_t   gl,
                   Complex_t   gr,
                   Double_t    m,
                   Double_t    gm);

   virtual ~HELScalar() {}

   inline const ANL4DVector &GetFourMomentum() const  { return fP;   }
   inline       Int_t        GetNSS         () const  { return fNSS; }
   
private:
   ANL4DVector fP;                // 4-momentum * fNSS
   Int_t       fNSS;              // (-1,1) = (initial, final)

   ClassDef(HELScalar, 1)  // Scalar boson class
};

//_______________________________________________________________________
// =====================
//  class HELVertex
// =====================
//-----------------------------------------------------------------------
class HELVertex: public Complex_t {
public:
   HELVertex(Complex_t val = 0.) : Complex_t(val) {}
   HELVertex(const HELFermion &in,
             const HELFermion &out,
             const HELVector  &v,
	           Double_t    gl,
		   Double_t    gr);
   HELVertex(const HELVector  &v1,
             const HELVector  &v2,
             const HELScalar  &sc,
                   Double_t    g);
   HELVertex(const HELVector  &v1,
             const HELVector  &v2,
             const HELScalar  &sc,
                   Double_t    g1,
                   Double_t    g2,
                   Double_t    g3);
   HELVertex(const HELVector  &v1,
             const HELVector  &v2,
             const HELVector  &v3,
                   Double_t    g);
   HELVertex(const HELVector  &wm,
             const HELVector  &w31,
             const HELVector  &wp,
             const HELVector  &w32,
                   Double_t    g31,
                   Double_t    g32,
                   Double_t    mw,
                   Double_t    gamw,
                   Bool_t      mode);
   HELVertex(const HELVector  &wm1,
             const HELVector  &wp1,
             const HELVector  &wm2,
             const HELVector  &wp2,
                   Double_t    gwwa,
                   Double_t    gwwz,
                   Double_t    mz,
                   Double_t    gamz);
   HELVertex(const HELVector  &v1,
             const HELVector  &v2,
             const HELVector  &v3,
             const HELVector  &v4,
                   Double_t     g);
   HELVertex(const HELVector  &vc,
             const HELScalar  &s1,
             const HELScalar  &s2,
	           Double_t    g);
   HELVertex(const HELVector  &v1,
             const HELVector  &v2,
             const HELScalar  &s1,
             const HELScalar  &s2,
	           Double_t    g);
   HELVertex(const HELFermion &in,
             const HELFermion &out,
             const HELScalar  &sc,
	           Complex_t   gl,
		   Complex_t   gr);
            
   virtual ~HELVertex() {}

   ClassDef(HELVertex, 1)  // Vertex class
};
#endif
