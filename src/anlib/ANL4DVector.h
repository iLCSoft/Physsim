#ifndef __ANLLVECTOR__
#define __ANLLVECTOR__
//*************************************************************************
//* ===================
//*  ANL4DVector Class
//* ===================
//*
//* (Description)
//*    A very primitive lockable Lorentz vector class.
//* (Requires)
//*	class TLorentzVector
//* 	class TAttLockable
//*     class ANL3DVector
//* (Provides)
//* 	class ANL4DVector
//* (Update Recored)
//*    1999/09/05  K.Ikematsu	Original version.
//*    2000/03/23  K.Ikematsu	Added Get3D method.
//*    2000/03/23  K.Ikematsu	Added GetTheta method.
//*    2000/03/28  K.Ikematsu	Added Acol method.
//*
//*************************************************************************
//
#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TAttLockable.h"
//#include "ANL3DVector.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------
//  Lockable TLorentzVector
//  -----------------------
//
class ANL4DVector : public TLorentzVector, public lcme::TAttLockable {
public:
   ANL4DVector(Double_t e=0., Double_t px=0., Double_t py=0., Double_t pz=0.) 
   			       : TLorentzVector(px,py,pz,e) {}
   ANL4DVector(Float_t e, Float_t px=0., Float_t py=0., Float_t pz=0.) {
     TLorentzVector::operator[](kT) = e;
     TLorentzVector::operator[](kX) = px;
     TLorentzVector::operator[](kY) = py;
     TLorentzVector::operator[](kZ) = pz;
   }
   //   ANL4DVector(const TVector &q) {
   //     TLorentzVector::operator[](kT) = q(0);
   //     TLorentzVector::operator[](kX) = q(1);
   //     TLorentzVector::operator[](kY) = q(2);
   //     TLorentzVector::operator[](kZ) = q(3);
   //   }

   ANL4DVector(const TLorentzVector &q) : TLorentzVector(q) {}
   ANL4DVector(const ANL4DVector &q) : TLorentzVector(q), lcme::TAttLockable(q) {}

   virtual ~ANL4DVector() {}

   inline Double_t & E()  { return TLorentzVector::operator[](kT); }
   inline Double_t & Px() { return TLorentzVector::operator[](kX); }
   inline Double_t & Py() { return TLorentzVector::operator[](kY); }
   inline Double_t & Pz() { return TLorentzVector::operator[](kZ); }
   inline Double_t & T()  { return TLorentzVector::operator[](kT); }
   inline Double_t & X()  { return TLorentzVector::operator[](kX); }
   inline Double_t & Y()  { return TLorentzVector::operator[](kY); }
   inline Double_t & Z()  { return TLorentzVector::operator[](kZ); }
   inline Double_t & operator()(Int_t i) {
     return TLorentzVector::operator[]((i>0 ? i-1 : 3));
   }

   inline Double_t E()  const { return TLorentzVector::operator[](kT); }
   inline Double_t Px() const { return TLorentzVector::operator[](kX); }
   inline Double_t Py() const { return TLorentzVector::operator[](kY); }
   inline Double_t Pz() const { return TLorentzVector::operator[](kZ); }
   inline Double_t T()  const { return TLorentzVector::operator[](kT); }
   inline Double_t X()  const { return TLorentzVector::operator[](kX); }
   inline Double_t Y()  const { return TLorentzVector::operator[](kY); }
   inline Double_t Z()  const { return TLorentzVector::operator[](kZ); }
   inline Double_t operator()(Int_t i) const {
     return TLorentzVector::operator[]((i>0 ? i-1 : 3));
   }

   inline friend ANL4DVector operator+ (const ANL4DVector &q1,
				       const ANL4DVector &q2) {
     ANL4DVector ans = q1; ans += q2; return ans;
   }
   inline friend ANL4DVector operator- (const ANL4DVector &q1,
				       const ANL4DVector &q2) {
     ANL4DVector ans = q1; ans -= q2; return ans;
   }
   inline friend ANL4DVector operator+ (const ANL4DVector &q1) { return q1; }
   inline friend ANL4DVector operator- (const ANL4DVector &q1) {
     ANL4DVector ans(-q1(0),-q1(1),-q1(2),-q1(3)); return ans;
   }
   inline friend Double_t operator* (const ANL4DVector &q1,
				     const ANL4DVector &q2) {
     return ( q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) );
   }

   //   inline ANL3DVector Get3D() const {
   //     ANL3DVector vec(operator()(1),operator()(2),operator()(3));
   //     return vec;
   //   }

   inline Double_t GetPt2()   const { return ( operator()(1)*operator()(1) 
                                       + operator()(2)*operator()(2) ); }
   inline Double_t GetMag2()  const { return ( GetPt2() 
                                       + operator()(3)*operator()(3) ); }
   inline Double_t GetPt()    const { return TMath::Sqrt( GetPt2() );       }
   inline Double_t GetMag()   const { return TMath::Sqrt( GetMag2() );      }
   inline Double_t GetMass2() const { return ( operator()(0)*operator()(0) 
                                       - GetMag2() );   }
   inline Double_t GetMass() const {
     return ( GetMass2() < 0  ? (-TMath::Sqrt(-GetMass2()))
                               : ( TMath::Sqrt( GetMass2())) );
   }
   inline Double_t GetTheta() const { return (180.*(TMath::ACos(CosTheta()))/TMath::Pi()); }
   inline Double_t GetTheta(const ANL4DVector &q) const {
     return (180.*(TMath::ACos(CosTheta(q)))/TMath::Pi());
   }
   inline Double_t Acol(const ANL4DVector &q) const {
     return (180.*(TMath::Pi()-TMath::ACos(CosTheta(q)))/TMath::Pi());
   }
   inline Double_t Acop(const ANL4DVector &q) const {
     Double_t c = (operator()(1)*q(1)+operator()(2)*q(2))
                  /(GetPt()*q.GetPt());
     return (180.*(TMath::Pi()-TMath::ACos(c))/TMath::Pi());
   }
#if 0
   inline Double_t CosTheta() const { return TLorentzVector::CosTheta(); }
#else
   inline Double_t CosTheta() const { 
     return (GetMag() == 0. ? 1.0 : operator()(3)/GetMag());
   }
#endif
   inline Double_t CosTheta(const ANL4DVector &q) const {
     return (operator()(1)*q(1)+operator()(2)*q(2)+operator()(3)*q(3))
                  /(GetMag()*q.GetMag());
   }

   inline virtual void DebugPrint(const Char_t *opt = "Brief") const { 
     cerr << "p    = " << operator()(0) << " " << operator()(1) << " "
                       << operator()(2) << " " << operator()(3) << endl;
     if (TString(opt).Contains("Detailed")) {
       cerr << "pt   = " << GetPt()   << endl;
       cerr << "ap   = " << GetMag()  << endl;
       cerr << "mass = " << GetMass() << endl;
     }
   }

   ClassDef(ANL4DVector,1)  // Lockable Lorentz vector class
};

#endif
