//*****************************************************************************
//* =====================
//*  HELLib 
//* =====================
//*  
//* (Description)
//*    Class library for helicity amplitude calculations.
//*
//* (Update Record)
//*    2007/01/27  K.Fujii	Original version based on HELAS.
//*
//*****************************************************************************

#include "HELLib.h"

#include <sstream>
#include <iomanip>
//#define __DEBUG__

  
ClassImp(HELFermion)
ClassImp(HELVector)
ClassImp(HELScalar)
ClassImp(HELVertex)

//-----------------------------------------------------------------------------
// ==============================
//  class HELFermion
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------------------
// IXXXXX() and OXXXXX()
//----------------------
HELFermion::HELFermion(const ANL4DVector &p,
                             Double_t     m,
                             Int_t        hel,
                             Int_t        nsf,
                             Bool_t       isincom)
          : TVectorC(4),
            fP(nsf*p.E(), nsf*p.Px(), nsf*p.Py(), nsf*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSF(nsf),
            fIsIncoming(isincom)
{
   Double_t  sf[2], sfomeg[2], omega[2], pp, pp3, sqpop3, sqm;
   Complex_t chi[2];
   Int_t nh = hel*nsf;
   Int_t in = isincom ? 1 : -1;
   if (m == 0.) {
      sqpop3 = nsf*TMath::Sqrt(TMath::Max(p(0)+p(3), 0.));
      chi[0] = sqpop3;
      if (sqpop3 == 0.) chi[1] = -hel * TMath::Sqrt(2.*p(0));
      else              chi[1] = Complex_t(nh*p(1), in*p(2))/sqpop3;
      Int_t iu = (1-in)/2;
      Int_t id = (1+in)/2;
      if (nh == 1*in) {
         (*this)[0] = 0.;
         (*this)[1] = 0.;
         (*this)[2] = chi[iu];
         (*this)[3] = chi[id];
      } else {
         (*this)[0] = chi[id];
         (*this)[1] = chi[iu];
         (*this)[2] = 0.;
         (*this)[3] = 0.;
      }
   } else {
      pp = TMath::Min(p(0), p.Vect().Mag());
      if (pp == 0.) {
         sqm = TMath::Sqrt(m);
         Int_t ip =    (1+in*nh)/2;
         Int_t im = in*(1-in*nh)/2;
         (*this)[0] = ip       * sqm;
         (*this)[1] = im * nsf * sqm;
         (*this)[2] = ip * nsf * sqm;
         (*this)[3] = im       * sqm;
      } else {
         sf[0] = (1+nsf+(1-nsf)*nh)*0.5;
         sf[1] = (1+nsf-(1-nsf)*nh)*0.5;
         omega[0] = TMath::Sqrt(p(0)+pp);
         omega[1] = m/omega[0];
         Int_t ip = (1+nh)/2;
         Int_t im = (1-nh)/2;
         sfomeg[0] = sf[0] * omega[ip];
         sfomeg[1] = sf[1] * omega[im];
         pp3    = TMath::Max(pp+p(3), 0.);
         chi[0] = TMath::Sqrt(pp3*0.5/pp);
         if (pp3 == 0.) chi[1] = -nh;
         else           chi[1] = Complex_t(nh*p(1), in*p(2))/TMath::Sqrt(2.*pp*pp3);
         Int_t iu = (1-in)/2;
         Int_t id = (1+in)/2;
         (*this)[0] = sfomeg[iu] * chi[im];
         (*this)[1] = sfomeg[iu] * chi[ip];
         (*this)[2] = sfomeg[id] * chi[im];
         (*this)[3] = sfomeg[id] * chi[ip];
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (isincom ? "fin" : "fot")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " nsf*fP = (" << nsf*fP.E () << ", "
                                 << nsf*fP.Px() << ", "
                                 << nsf*fP.Py() << ", "
                                 << nsf*fP.Pz() << ") " << endl;
#endif
}

//----------------------
// FVIXXX() or FVOXXX()
//----------------------
HELFermion::HELFermion(const HELFermion &f,
                       const HELVector  &v,
                             Double_t    gl,
                             Double_t    gr,
                             Double_t    m,
                             Double_t    gam)
          : TVectorC(4),
            fP((f.fIsIncoming ? f.fP - v.fP : f.fP + v.fP)),
	    fM(m),
            fHel(0),
            fNSF(0),
            fIsIncoming(f.fIsIncoming)
{
   static const Complex_t kI(0., 1.);
   
   Double_t  pf2 = fP.Mag2();
   Complex_t d   = -1./Complex_t(pf2-m*m, TMath::Max(TMath::Sign(m*gam, pf2), 0.));
   if (fIsIncoming) {
      Complex_t sl1 = (v[0] +    v[3])*f[0]
                    + (v[1] - kI*v[2])*f[1];
      Complex_t sl2 = (v[1] + kI*v[2])*f[0]
                    + (v[0] -    v[3])*f[1];
      if (gr == 0.) {
         (*this)[0] = gl * ((fP(0) - fP(3))*sl1 - (fP(1) - fP(2)*kI)*sl2)*d;
         (*this)[1] = gl * ((fP(0) + fP(3))*sl2 - (fP(1) + fP(2)*kI)*sl1)*d;
         (*this)[2] = gl * m * sl1 * d;
         (*this)[3] = gl * m * sl2 * d;
      } else {
         Complex_t sr1 =   (v[0] -    v[3])*f[2]
                         - (v[1] - kI*v[2])*f[3];
         Complex_t sr2 = - (v[1] + kI*v[2])*f[2]
                         + (v[0] +    v[3])*f[3];
         (*this)[0] = ( gl * ((fP(0) - fP(3))*sl1 - (fP(1) - fP(2)*kI)*sl2)
                      + gr * m * sr1) * d;
         (*this)[1] = ( gl * ((fP(0) + fP(3))*sl2 - (fP(1) + fP(2)*kI)*sl1)
                      + gr * m * sr2) * d;
         (*this)[2] = ( gr * ((fP(0) + fP(3))*sr1 + (fP(1) - fP(2)*kI)*sr2)
                      + gl * m * sl1) * d;
         (*this)[3] = ( gr * ((fP(0) - fP(3))*sr2 + (fP(1) + fP(2)*kI)*sr1)
                      + gl * m * sl2) * d;
      }
   } else {
      Complex_t sl1 = (v[0] +    v[3])*f[2]
                    + (v[1] + kI*v[2])*f[3];
      Complex_t sl2 = (v[1] - kI*v[2])*f[2]
                    + (v[0] -    v[3])*f[3];
      if (gr == 0.) {
         (*this)[0] = gl * m * sl1 * d;
         (*this)[1] = gl * m * sl2 * d;
         (*this)[2] = gl * ((fP(0) - fP(3))*sl1 - (fP(1) + fP(2)*kI)*sl2)*d;
         (*this)[3] = gl * ((fP(0) + fP(3))*sl2 - (fP(1) - fP(2)*kI)*sl1)*d;
      } else {
         Complex_t sr1 =   (v[0] -    v[3])*f[0]
                         - (v[1] + kI*v[2])*f[1];
         Complex_t sr2 = - (v[1] - kI*v[2])*f[0]
                         + (v[0] +    v[3])*f[1];
         (*this)[0] = ( gr * ((fP(0) + fP(3))*sr1 + (fP(1) + fP(2)*kI)*sr2)
                      + gl * m * sl1) * d;
         (*this)[1] = ( gr * ((fP(0) - fP(3))*sr2 + (fP(1) - fP(2)*kI)*sr1)
                      + gl * m * sl2) * d;
         (*this)[2] = ( gl * ((fP(0) - fP(3))*sl1 - (fP(1) + fP(2)*kI)*sl2)
                      + gr * m * sr1) * d;
         (*this)[3] = ( gl * ((fP(0) + fP(3))*sl2 - (fP(1) - fP(2)*kI)*sl1)
                      + gr * m * sr2) * d;
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (fIsIncoming ? "fvi" : "fvo")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELVector
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// VXXXXX()
//----------
HELVector::HELVector(const ANL4DVector &p,
                           Double_t     m,
                           Int_t        hel,
                           Int_t        nsv)
          : TVectorC(4),
            fP(nsv*p.E(), nsv*p.Px(), nsv*p.Py(), nsv*p.Pz()),
	    fM(m),
            fHel(hel),
            fNSV(nsv)
{
   static const Double_t kSqh = TMath::Sqrt(0.5);
   Double_t  hel0, pt, pt2, pp, pzpt, emp;
   Int_t     nsvahl;
   nsvahl = nsv*TMath::Abs(hel);
   pt2 = p.GetPt2();
   pp  = TMath::Min(p(0), p.Vect().Mag());
   pt  = TMath::Min(pp, TMath::Sqrt(pt2));
   if (m == 0.) {
      pp = p(0);
      pt = p.GetPt();
      (*this)[0] = 0.;
      (*this)[3] = hel*pt/pp*kSqh;
      if (pt != 0.) {
         pzpt = p(3)/(pp*pt)*kSqh*hel;
         (*this)[1] = Complex_t(-p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
         (*this)[2] = Complex_t(-p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
      } else {
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
      }
   } else {
      hel0 = 1. - TMath::Abs(hel);
      if (pp == 0.) {
         (*this)[0] = 0.;
         (*this)[1] = -hel*kSqh;
         (*this)[2] = Complex_t(0., nsvahl*kSqh);
         (*this)[3] = hel0;
      } else {
         emp = p(0)/(m*pp);
	 (*this)[0] = hel0*pp/m;
         (*this)[3] = hel0*p(3)*emp + hel*pt/pp*kSqh;
	 if (pt != 0.) {
            pzpt = p(3)/(pp*pt)*kSqh*hel;
	    (*this)[1] = Complex_t(hel0*p(1)*emp - p(1)*pzpt, -nsvahl*p(2)/pt*kSqh);
	    (*this)[2] = Complex_t(hel0*p(2)*emp - p(2)*pzpt,  nsvahl*p(1)/pt*kSqh);
	 } else {
	    (*this)[1] = -hel*kSqh;
	    (*this)[2] = Complex_t(0., nsvahl*TMath::Sign(kSqh,p(3)));
	 }
      }
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << (nsv < 0 ? "vin" : "vot")
		<<  " =(" << (*this)[0] << ", "
                          << (*this)[1] << ", "
                          << (*this)[2] << ", "
                          << (*this)[3] << ") " << endl;
           cerr << " nsv*fP = (" << nsv*fP.E () << ", "
                                 << nsv*fP.Px() << ", "
                                 << nsv*fP.Py() << ", "
                                 << nsv*fP.Pz() << ") " << endl;
#endif
}

HELVector::HELVector(Complex_t v0,
                     Complex_t v1,
                     Complex_t v2,
                     Complex_t v3,
                     const ANL4DVector &p)
	  : fM(0),
            fHel(0),
            fNSV(0)
{
   (*this)[0] = v0;
   (*this)[1] = v1;
   (*this)[2] = v2;
   (*this)[3] = v3;
   fP = p;
}

//----------
// JIOXXX()
//----------
HELVector::HELVector(const HELFermion &fin,
                     const HELFermion &fout,
                           Double_t    glv,
                           Double_t    grv,
                           Double_t    mv,
                           Double_t    gamv)
          : TVectorC(4),
            fP(fout.fP - fin.fP),
	    fM(mv),
            fHel(0),
            fNSV(0)
{
   Complex_t c0, c1, c2, c3, cs, d;
   Double_t  q2, vm2, dd;
   q2  = fP.Mag2();
   vm2 = mv*mv;
   if (mv == 0.) {
      dd = 1./q2;
      if (grv == 0.) { // purely left-handed
         dd *= glv;
	 (*this)[0] = ( fout[2]*fin[0] + fout[3]*fin[1]) * dd;
	 (*this)[1] = (-fout[2]*fin[1] - fout[3]*fin[0]) * dd;
	 (*this)[2] = ( fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., dd);
	 (*this)[3] = (-fout[2]*fin[0] + fout[3]*fin[1]) * dd;
      } else {
         (*this)[0] = (  glv * ( fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] + fout[1]*fin[3])) * dd;
         (*this)[1] = (- glv * ( fout[2]*fin[1] + fout[3]*fin[0])
                       + grv * ( fout[0]*fin[3] + fout[1]*fin[2])) * dd;
         (*this)[2] = (  glv * ( fout[2]*fin[1] - fout[3]*fin[0])
                       + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., dd);
         (*this)[3] = (  glv * (-fout[2]*fin[0] + fout[3]*fin[1])
                       + grv * ( fout[0]*fin[2] - fout[1]*fin[3])) * dd;
      }
   } else {
      d = 1./Complex_t(q2 - vm2, TMath::Max(TMath::Sign(mv*gamv, q2), 0.));
      if (grv == 0.) { // purely left-handed
         d *= glv;
	 c0 =  fout[2]*fin[0] + fout[3]*fin[1];
	 c1 = -fout[2]*fin[1] - fout[3]*fin[0];
	 c2 = (fout[2]*fin[1] - fout[3]*fin[0]) * Complex_t(0., 1.);
	 c3 = -fout[2]*fin[0] + fout[3]*fin[1];
      } else {
         c0 =   glv * ( fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] + fout[1]*fin[3]);
         c1 = - glv * ( fout[2]*fin[1] + fout[3]*fin[0])
              + grv * ( fout[0]*fin[3] + fout[1]*fin[2]);
         c2 =  (glv * ( fout[2]*fin[1] - fout[3]*fin[0])
              + grv * (-fout[0]*fin[3] + fout[1]*fin[2])) * Complex_t(0., 1.);
         c3 =   glv * (-fout[2]*fin[0] + fout[3]*fin[1])
              + grv * ( fout[0]*fin[2] - fout[1]*fin[3]);
      }
      cs = (fP(0)*c0 - fP(1)*c1 - fP(2)*c2 - fP(3)*c3)/vm2;
      (*this)[0] = (c0 - cs*fP(0))*d;
      (*this)[1] = (c1 - cs*fP(1))*d;
      (*this)[2] = (c2 - cs*fP(2))*d;
      (*this)[3] = (c3 - cs*fP(3))*d;
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " jio =(" << (*this)[0] << ", "
                             << (*this)[1] << ", "
                             << (*this)[2] << ", "
                             << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}

//----------
// J3XXXX()
//----------
HELVector::HELVector(const HELFermion &fin,
                     const HELFermion &fout,
                           Double_t    gla,
                           Double_t    gra,
                           Double_t    glz,
                           Double_t    grz,
                           Double_t    mz,
                           Double_t    gamz)
          : TVectorC(4),
            fP(fout.fP - fin.fP),
	    fM(mz),
            fHel(0),
            fNSV(0)
{
   Double_t q2  = fP.Mag2();
   Double_t zm2 = mz*mz;
   Double_t zmw = mz*gamz;

   Double_t  da   = 1./q2;
   Double_t  ww   = TMath::Max(TMath::Sign(zmw, q2), 0.);
   Complex_t dz   = 1./Complex_t(q2 - zm2, ww);
   Complex_t ddif = Complex_t(-zm2, ww) * da * dz;
   Double_t  cw   = 1./TMath::Sqrt(1.+ TMath::Power(grz/gra,2));
   Double_t  sw   = TMath::Sqrt((1.-cw)*(1+cw));
   Double_t  gn   = gra * sw;
   Double_t  gz3l = glz * cw;
   Double_t  ga3l = gla * sw;

   static const Complex_t kI(0., 1.);
   Complex_t c0l  =   fout[2] * fin[0] + fout[3] * fin[1];
   Complex_t c0r  =   fout[0] * fin[2] + fout[1] * fin[3];
   Complex_t c1l  = -(fout[2] * fin[1] + fout[3] * fin[0]);
   Complex_t c1r  =   fout[0] * fin[3] + fout[1] * fin[2];
   Complex_t c2l  =  (fout[2] * fin[1] - fout[3] * fin[0])*kI;
   Complex_t c2r  = (-fout[0] * fin[3] + fout[1] * fin[2])*kI;
   Complex_t c3l  =  -fout[2] * fin[0] + fout[3] * fin[1];
   Complex_t c3r  =   fout[0] * fin[2] - fout[1] * fin[3];
   Complex_t csl  = (-fP(0)*c0l+fP(1)*c1l+fP(2)*c2l+fP(3)*c3l)/zm2; 
   Complex_t csr  = (-fP(0)*c0r+fP(1)*c1r+fP(2)*c2r+fP(3)*c3r)/zm2; 

   (*this)[0] = gz3l * dz * (c0l  + csl * fP(0)) + ga3l * c0l *da
               + gn * (c0r * ddif - csr * fP(0)*dz);
   (*this)[1] = gz3l * dz * (c1l  + csl * fP(1)) + ga3l * c1l *da
               + gn * (c1r * ddif - csr * fP(1)*dz);
   (*this)[2] = gz3l * dz * (c2l  + csl * fP(2)) + ga3l * c2l *da
               + gn * (c2r * ddif - csr * fP(2)*dz);
   (*this)[3] = gz3l * dz * (c3l  + csl * fP(3)) + ga3l * c3l *da
               + gn * (c3r * ddif - csr * fP(3)*dz);
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " j3 =(" << (*this)[0] << ", "
                            << (*this)[1] << ", "
                            << (*this)[2] << ", "
                            << (*this)[3] << ") " << endl;
           cerr << " fP = (" << fP.E () << ", "
                             << fP.Px() << ", "
                             << fP.Py() << ", "
                             << fP.Pz() << ") " << endl;
#endif
}
//----------
// JVSXXX()
//----------
HELVector::HELVector(const HELVector  &vc,
                     const HELScalar  &sc,
                           Double_t    g,
                           Double_t    mv,
                           Double_t    gmv)
          : TVectorC(4),
            fP(vc.fP + sc.fP),
	    fM(mv),
            fHel(0),
            fNSV(0)
{
   Double_t  q2  = fP.Mag2();
   if (mv == 0.) {
      Complex_t dg = g*sc/q2;
      (*this)[0] = dg * vc[0];
      (*this)[1] = dg * vc[1];
      (*this)[2] = dg * vc[2];
      (*this)[3] = dg * vc[3];
   } else {
      Double_t  vm2 = mv*mv;
      Double_t  mg  = TMath::Max(TMath::Sign(mv*gmv, q2), 0.);
      Complex_t dg  = g*sc/Complex_t(q2 - vm2, mg);
      Complex_t vk  = (-fP(0)*vc[0] + fP(1)*vc[1] + fP(2)*vc[2] + fP(3)*vc[3])/vm2;
      (*this)[0] = dg * (fP(0)*vk + vc[0]);
      (*this)[1] = dg * (fP(1)*vk + vc[1]);
      (*this)[2] = dg * (fP(2)*vk + vc[2]);
      (*this)[3] = dg * (fP(3)*vk + vc[3]);
   } 
}
//----------
// JVVXXX()
//----------
HELVector::HELVector(const HELVector  &v1,
                     const HELVector  &v2,
                           Double_t    g,
                           Double_t    mv,
                           Double_t    gmv)
          : TVectorC(4),
            fP(v1.fP + v2.fP),
	    fM(mv),
            fHel(0),
            fNSV(0)
{
   Double_t  s   = fP.Mag2();
   Double_t  vm2 = mv*mv;
   Complex_t v12 = v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3];
   Complex_t sv1 =  (v2.fP(0)+fP(0))*v1[0] -(v2.fP(1)+fP(1))*v1[1]
                                           -(v2.fP(2)+fP(2))*v1[2] 
                                           -(v2.fP(3)+fP(3))*v1[3];
   Complex_t sv2 = -(v1.fP(0)+fP(0))*v2[0] +(v1.fP(1)+fP(1))*v2[1]
                                           +(v1.fP(2)+fP(2))*v2[2] 
                                           +(v1.fP(3)+fP(3))*v2[3];
   Complex_t j12[4];
   j12[0] = (v1.fP(0)-v2.fP(0))*v12 + sv1*v2[0] + sv2*v1[0];
   j12[1] = (v1.fP(1)-v2.fP(1))*v12 + sv1*v2[1] + sv2*v1[1];
   j12[2] = (v1.fP(2)-v2.fP(2))*v12 + sv1*v2[2] + sv2*v1[2];
   j12[3] = (v1.fP(3)-v2.fP(3))*v12 + sv1*v2[3] + sv2*v1[3];

   if (mv == 0.) {
      Double_t gs = -g/s;
      (*this)[0] = gs * j12[0];
      (*this)[1] = gs * j12[1];
      (*this)[2] = gs * j12[2];
      (*this)[3] = gs * j12[3];
   } else {
      Double_t  m1  = v1.fP.Mag2();
      Double_t  m2  = v2.fP.Mag2();
      Complex_t s11 = v1.fP(0)*v1[0] - v1.fP(1)*v1[1] - v1.fP(2)*v1[2] - v1.fP(3)*v1[3];
      Complex_t s12 = v1.fP(0)*v2[0] - v1.fP(1)*v2[1] - v1.fP(2)*v2[2] - v1.fP(3)*v2[3];
      Complex_t s21 = v2.fP(0)*v1[0] - v2.fP(1)*v1[1] - v2.fP(2)*v1[2] - v2.fP(3)*v1[3];
      Complex_t s22 = v2.fP(0)*v2[0] - v2.fP(1)*v2[1] - v2.fP(2)*v2[2] - v2.fP(3)*v2[3];
      Complex_t js  = (v12*(-m1+m2) + s11*s12 - s21*s22)/vm2;
      Double_t  mg  = TMath::Max(TMath::Sign(mv*gmv, s), 0.);
      Complex_t dg  = -g/Complex_t(s - vm2, mg);
      (*this)[0] = dg * (j12[0] + fP(0)*js);
      (*this)[1] = dg * (j12[1] + fP(1)*js);
      (*this)[2] = dg * (j12[2] + fP(2)*js);
      (*this)[3] = dg * (j12[3] + fP(3)*js);
   }
}

//----------
// JEEXXX()
//----------
HELVector::HELVector (Double_t ebm, 
                      Double_t eef,
	              Double_t sh,
	              Double_t ch,
	              Double_t fi,
	              Int_t    helbm,
	              Int_t    helef,
	              Int_t    nsf,
		      Double_t ge,
		      Double_t me)
          : TVectorC(4),
	    fM(0.),
            fHel(0),
            fNSV(0)
{
   Double_t hi  = helbm;
   Double_t sf  = nsf;
   Double_t sfh = helbm*nsf;
   Double_t cs[2];
   cs[(1+nsf)/2] = sh;
   cs[(1-nsf)/2] = ch;

   Double_t x   = eef/ebm;
   Double_t me2 = me*me;
   Double_t q2  = -4.*cs[1]*cs[1]*(eef*ebm-me2)
	          + sf*(1.-x)*(1-x)/x*(sh+ch)*(sh-ch)*me2;
   Double_t rfp = 1 + nsf;
   Double_t rfm = 1 - nsf;
   Double_t cfi = TMath::Cos(fi);
   Double_t sfi = TMath::Sin(fi);

   if (helbm == helef) {
      Double_t  rxc   = 2.*x/(1.-x)*cs[0]*cs[0];
      Complex_t coeff = ge*2.*ebm*TMath::Sqrt(x)*cs[1]/q2
	               *(Complex_t(rfp,0.)-rfm*Complex_t(cfi, -sfi*hi))*0.5;
      (*this)[0] =  Complex_t(0., 0.);
      (*this)[1] =  coeff*Complex_t((1.+rxc)*cfi, -sfh*sfi);
      (*this)[2] =  coeff*Complex_t((1.+rxc)*sfi,  sfh*cfi);
      (*this)[3] =  coeff*(-sf*rxc/cs[0]*cs[1]);
   } else {
      Complex_t coeff = ge*me/q2/TMath::Sqrt(x)
	               *(Complex_t(rfp,0.)+rfm*Complex_t(cfi, sfi*hi))*0.5*hi;
      (*this)[0] = -coeff*(1.+x)*cs[1]*Complex_t(cfi, sfh*sfi);
      (*this)[1] =  coeff*(1.-x)*cs[0];
      (*this)[2] =  (*this)[1]*Complex_t(0., sfh);
      (*this)[3] =  (*this)[0]*sf*(1.-x)/(1.+x);
   }
   Double_t cth = (ch+sh)*(ch-sh);
   Double_t sth  = 2.*sh*ch;

   fP(0) = -ebm*(1. - x);
   fP(1) = ebm*x*sth*cfi;
   fP(2) = ebm*x*sth*sfi;
   fP(3) = -ebm*(sf - x*cth);
}

void HELVector::DebugPrint() const
{
   cerr << "p=(" << fP(0) << ", "
                 << fP(1) << ", "
                 << fP(2) << ", "
                 << fP(3) << ") " << endl;
   cerr << "v=(" << (*this)[0] << ", "
	         << (*this)[1] << ", "
	         << (*this)[2] << ", "
	         << (*this)[3] << ") " << endl;
}


//-----------------------------------------------------------------------------
// ==============================
//  class HELScalar
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// SXXXXX()
//----------
HELScalar::HELScalar(const ANL4DVector &p,
                           Int_t        nss)
          : Complex_t(1., 0.),
            fP(nss*p.E(), nss*p.Px(), nss*p.Py(), nss*p.Pz()),
            fNSS(nss)
{
}

//----------
// HVSXXX()
//----------
HELScalar::HELScalar(const HELVector  &vc,
                     const HELScalar  &sc,
	                   Double_t    g,
                           Double_t    m,
                           Double_t    gm)
          : fP(vc.fP + sc.fP),
            fNSS(1)
{
   Double_t  q2  = fP.Mag2();
   Double_t  mg  = TMath::Max(TMath::Sign(m*gm, q2), 0.);
   Complex_t dg  = -g/Complex_t(q2 - m*m, mg);
   Complex_t qvv = vc.fP(0)*vc[0] - vc.fP(1)*vc[1] - vc.fP(2)*vc[2] - vc.fP(3)*vc[3];
   Complex_t qpv = sc.fP(0)*vc[0] - sc.fP(1)*vc[1] - sc.fP(2)*vc[2] - sc.fP(3)*vc[3];
   this->Complex_t::operator=(dg*(2.*qpv+qvv)*sc);
}

//----------
// HSSXXX()
//----------
HELScalar::HELScalar(const HELScalar  &s1,
                     const HELScalar  &s2,
	                   Double_t    g,
                           Double_t    m,
                           Double_t    gm)
          : fP(s1.fP + s2.fP),
            fNSS(1)
{
   Double_t  q2  = fP.Mag2();
   Double_t  mg  = TMath::Max(TMath::Sign(m*gm, q2), 0.);
   Complex_t dg  = -g/Complex_t(q2 - m*m, mg);
   this->Complex_t::operator=(dg*s1*s2);
}

//----------
// HVVXXX()
//----------
HELScalar::HELScalar(const HELVector  &v1,
                     const HELVector  &v2,
                           Double_t    g,
                           Double_t    m,
                           Double_t    gm)
          : fP(v1.fP + v2.fP),
            fNSS(1)
{
   Double_t  q2  = fP.Mag2();
   Double_t  mg  = TMath::Max(TMath::Sign(m*gm, q2), 0.);
   Complex_t dg  = -g/Complex_t(q2 - m*m, mg);

   this->Complex_t::operator=(dg * (v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3]));
}
// For anomalous VVh coupling
HELScalar::HELScalar(const HELVector  &v1,
                     const HELVector  &v2,
                           Double_t    g1, // gvh + 2*mv^2*(a/Lamda)
                           Double_t    g2, // -2*(b/Lambda)
                           Double_t    g3, // -4*(btilde/Lambda)
                           Double_t    m,
                           Double_t    gm)
          : fP(v1.fP + v2.fP),
            fNSS(1)
{
   Double_t  q2  = fP.Mag2();
   Double_t  mg  = TMath::Max(TMath::Sign(m*gm, q2), 0.);
   Complex_t d  = Complex_t(-1.,0.)/Complex_t(q2 - m*m, mg);

   Double_t  p1p2 = v1.fP*v2.fP;
   Complex_t p1v2 = v1.fP(0)*v2[0]-v1.fP(1)*v2[1]-v1.fP(2)*v2[2]-v1.fP(3)*v2[3];
   Complex_t p2v1 = v2.fP(0)*v1[0]-v2.fP(1)*v1[1]-v2.fP(2)*v1[2]-v2.fP(3)*v1[3];
   Complex_t v1v2 = v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];

   Complex_t ans  = g1*v1v2;
             ans += g2*(p1p2*v1v2-p1v2*p2v1);
             ans += g3*(v1.fP(0)*v1[1]*v2.fP(2)*v2[3]
                       -v1.fP(0)*v1[1]*v2.fP(3)*v2[2]
                       -v1.fP(0)*v1[2]*v2.fP(1)*v2[3]
                       +v1.fP(0)*v1[2]*v2.fP(3)*v2[1]
                       +v1.fP(0)*v1[3]*v2.fP(1)*v2[2]
                       -v1.fP(0)*v1[3]*v2.fP(2)*v2[1]
                       -v1.fP(1)*v1[0]*v2.fP(2)*v2[3]
                       +v1.fP(1)*v1[0]*v2.fP(3)*v2[2]
                       +v1.fP(1)*v1[2]*v2.fP(0)*v2[3]
                       -v1.fP(1)*v1[2]*v2.fP(3)*v2[0]
                       -v1.fP(1)*v1[3]*v2.fP(0)*v2[2]
                       +v1.fP(1)*v1[3]*v2.fP(2)*v2[0]
                       +v1.fP(2)*v1[0]*v2.fP(1)*v2[3]
                       -v1.fP(2)*v1[0]*v2.fP(3)*v2[1]
                       -v1.fP(2)*v1[1]*v2.fP(0)*v2[3]
                       +v1.fP(2)*v1[1]*v2.fP(3)*v2[0]
                       +v1.fP(2)*v1[3]*v2.fP(0)*v2[1]
                       -v1.fP(2)*v1[3]*v2.fP(1)*v2[0]
                       -v1.fP(3)*v1[0]*v2.fP(1)*v2[2]
                       +v1.fP(3)*v1[0]*v2.fP(2)*v2[1]
                       +v1.fP(3)*v1[1]*v2.fP(0)*v2[2]
                       -v1.fP(3)*v1[1]*v2.fP(2)*v2[0]
                       -v1.fP(3)*v1[2]*v2.fP(0)*v2[1]
                       +v1.fP(3)*v1[2]*v2.fP(1)*v2[0]);
             ans *= d;
   this->Complex_t::operator=(ans);
}

//----------
// HIOXXX()
//----------
HELScalar::HELScalar(const HELFermion &fin,
                     const HELFermion &fout,
                           Complex_t   gl,
                           Complex_t   gr,
                           Double_t    m,
                           Double_t    gm)
          : fP(fout.fP - fin.fP),
	    fNSS(1)
{
   Double_t  q2  = fP.Mag2();
   Double_t  mg  = TMath::Max(TMath::Sign(m*gm, q2), 0.);
   Complex_t dn  = Complex_t(q2 - m*m, mg);

   this->Complex_t::operator=((gl*(fout[0]*fin[0]+fout[1]*fin[1]) 
		              +gr*(fout[2]*fin[2]+fout[3]*fin[3]))/dn);
}

//-----------------------------------------------------------------------------
// ==============================
//  class HELVertex
// ==============================
//_____________________________________________________________________________
//_____________________________________________________________________________
// --------------------------
//  C-tor
// --------------------------
//----------
// IOVXXX()
//----------
HELVertex::HELVertex(const HELFermion &in,
                     const HELFermion &out,
                     const HELVector  &v,
                           Double_t    gl,
                           Double_t    gr)
{
   
   *this    =  gl * ( (out[2]*in[0] + out[3]*in[1]) * v[0]
                     +(out[2]*in[1] + out[3]*in[0]) * v[1]
                     -(out[2]*in[1] - out[3]*in[0]) * v[2] * Complex_t(0., 1.)
                     +(out[2]*in[0] - out[3]*in[1]) * v[3] );
   if (gr != 0.) {
      *this += gr * ( (out[0]*in[2] + out[1]*in[3]) * v[0]
                     -(out[0]*in[3] + out[1]*in[2]) * v[1]
                     +(out[0]*in[3] - out[1]*in[2]) * v[2] * Complex_t(0., 1.)
                     -(out[0]*in[2] - out[1]*in[3]) * v[3] );
   }
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " iov =(" << (*this) << endl;
#endif
}

//----------
// VVSXXX()
//----------
HELVertex::HELVertex(const HELVector &v1,
                     const HELVector &v2,
                     const HELScalar &sc,
                           Double_t    g)
{
   *this = g * sc * (v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3]);
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " vvs =(" << (*this) << endl;
#endif
}
// For anomalous VVh coupling
HELVertex::HELVertex(const HELVector  &v1,
                     const HELVector  &v2,
                     const HELScalar  &sc,
                           Double_t    g1, // gvh + 2*mv^2*(a/Lamda)
                           Double_t    g2, // -2*(b/Lambda)
                           Double_t    g3) // -4*(btilde/Lambda)
{
   Double_t  p1p2 = v1.fP*v2.fP;
   Complex_t p1v2 = v1.fP(0)*v2[0]-v1.fP(1)*v2[1]-v1.fP(2)*v2[2]-v1.fP(3)*v2[3];
   Complex_t p2v1 = v2.fP(0)*v1[0]-v2.fP(1)*v1[1]-v2.fP(2)*v1[2]-v2.fP(3)*v1[3];
   Complex_t v1v2 = v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];

   Complex_t ans  = g1*v1v2;
             ans += g2*(p1p2*v1v2-p1v2*p2v1);
             ans += g3*(v1.fP(0)*v1[1]*v2.fP(2)*v2[3]
                       -v1.fP(0)*v1[1]*v2.fP(3)*v2[2]
                       -v1.fP(0)*v1[2]*v2.fP(1)*v2[3]
                       +v1.fP(0)*v1[2]*v2.fP(3)*v2[1]
                       +v1.fP(0)*v1[3]*v2.fP(1)*v2[2]
                       -v1.fP(0)*v1[3]*v2.fP(2)*v2[1]
                       -v1.fP(1)*v1[0]*v2.fP(2)*v2[3]
                       +v1.fP(1)*v1[0]*v2.fP(3)*v2[2]
                       +v1.fP(1)*v1[2]*v2.fP(0)*v2[3]
                       -v1.fP(1)*v1[2]*v2.fP(3)*v2[0]
                       -v1.fP(1)*v1[3]*v2.fP(0)*v2[2]
                       +v1.fP(1)*v1[3]*v2.fP(2)*v2[0]
                       +v1.fP(2)*v1[0]*v2.fP(1)*v2[3]
                       -v1.fP(2)*v1[0]*v2.fP(3)*v2[1]
                       -v1.fP(2)*v1[1]*v2.fP(0)*v2[3]
                       +v1.fP(2)*v1[1]*v2.fP(3)*v2[0]
                       +v1.fP(2)*v1[3]*v2.fP(0)*v2[1]
                       -v1.fP(2)*v1[3]*v2.fP(1)*v2[0]
                       -v1.fP(3)*v1[0]*v2.fP(1)*v2[2]
                       +v1.fP(3)*v1[0]*v2.fP(2)*v2[1]
                       +v1.fP(3)*v1[1]*v2.fP(0)*v2[2]
                       -v1.fP(3)*v1[1]*v2.fP(2)*v2[0]
                       -v1.fP(3)*v1[2]*v2.fP(0)*v2[1]
                       +v1.fP(3)*v1[2]*v2.fP(1)*v2[0]);
   this->Complex_t::operator=(ans);
}

//----------
// VVVXXX()
//----------
HELVertex::HELVertex(const HELVector &wm,
                     const HELVector &wp,
                     const HELVector &w3,
                           Double_t    g)
{
   Complex_t v12 = wm[0]*wp[0] - wm[1]*wp[1] - wm[2]*wp[2] - wm[3]*wp[3];
   Complex_t v23 = wp[0]*w3[0] - wp[1]*w3[1] - wp[2]*w3[2] - wp[3]*w3[3];
   Complex_t v31 = w3[0]*wm[0] - w3[1]*wm[1] - w3[2]*wm[2] - w3[3]*wm[3];
   Complex_t xv1 = 0.;
   Complex_t xv2 = 0.;
   Complex_t xv3 = 0.;
   if (abs(wm[0]) != 0. 
       && abs(wm[0]) >= TMath::Max(TMath::Max(abs(wm[1]),
                                              abs(wm[2])),
                                              abs(wm[3]))*1.e-1) {
      xv1 = wm.fP(0)/wm[0];
   }
   if (abs(wp[0]) != 0. 
       && abs(wp[0]) >= TMath::Max(TMath::Max(abs(wp[1]),
                                              abs(wp[2])),
                                              abs(wp[3]))*1.e-1) {
      xv2 = wp.fP(0)/wp[0];
   }
   if (abs(w3[0]) != 0. 
       && abs(w3[0]) >= TMath::Max(TMath::Max(abs(w3[1]),
                                              abs(w3[2])),
                                              abs(w3[3]))*1.e-1) {
      xv3 = w3.fP(0)/w3[0];
   }
   Complex_t p12 = (wm.fP(0) - xv1*wm[0])*wp[0] - (wm.fP(1) - xv1*wm[1])*wp[1]
                 - (wm.fP(2) - xv1*wm[2])*wp[2] - (wm.fP(3) - xv1*wm[3])*wp[3];
   Complex_t p13 = (wm.fP(0) - xv1*wm[0])*w3[0] - (wm.fP(1) - xv1*wm[1])*w3[1]
                 - (wm.fP(2) - xv1*wm[2])*w3[2] - (wm.fP(3) - xv1*wm[3])*w3[3];
   Complex_t p21 = (wp.fP(0) - xv2*wp[0])*wm[0] - (wp.fP(1) - xv2*wp[1])*wm[1]
                 - (wp.fP(2) - xv2*wp[2])*wm[2] - (wp.fP(3) - xv2*wp[3])*wm[3];
   Complex_t p23 = (wp.fP(0) - xv2*wp[0])*w3[0] - (wp.fP(1) - xv2*wp[1])*w3[1]
                 - (wp.fP(2) - xv2*wp[2])*w3[2] - (wp.fP(3) - xv2*wp[3])*w3[3];
   Complex_t p31 = (w3.fP(0) - xv3*w3[0])*wm[0] - (w3.fP(1) - xv3*w3[1])*wm[1]
                 - (w3.fP(2) - xv3*w3[2])*wm[2] - (w3.fP(3) - xv3*w3[3])*wm[3];
   Complex_t p32 = (w3.fP(0) - xv3*w3[0])*wp[0] - (w3.fP(1) - xv3*w3[1])*wp[1]
                 - (w3.fP(2) - xv3*w3[2])*wp[2] - (w3.fP(3) - xv3*w3[3])*wp[3];
   *this = - (v12*(p13-p23) + v23*(p21-p31) + v31*(p32-p12)) * g;
#ifdef __DEBUG__ 
           cerr << " ---- " << endl;
	   cerr << " vvv =(" << (*this) << endl;
#endif
}

//----------
// W3W3XX()
//----------
// The possible sets of the inputs are as follows:                       
//    -------------------------------------------                         
//    |  WM  |  W31 |  WP  |  W32 |  G31 |  G32 |                         
//    -------------------------------------------                         
//    |  W-  |  W3  |  W+  |  W3  |  GW  |  GW  |                         
//    |  W-  |  W3  |  W+  |  Z   |  GW  | GWWZ |                         
//    |  W-  |  W3  |  W+  |  A   |  GW  | GWWA |                         
//    |  W-  |  Z   |  W+  |  Z   | GWWZ | GWWZ |                         
//    |  W-  |  Z   |  W+  |  A   | GWWZ | GWWA |                         
//    |  W-  |  A   |  W+  |  A   | GWWA | GWWA |                         
//    -------------------------------------------                         
// where all the bosons are defined by the flowing-OUT quantum number.   

HELVertex::HELVertex(const HELVector  &wm,
                     const HELVector  &w31,
                     const HELVector  &wp,
                     const HELVector  &w32,
                           Double_t    g31,
                           Double_t    g32,
                           Double_t    mw,
                           Double_t    gamw,
                           Bool_t      mode)
{
   mode = mode;
   Complex_t v12 = wm [0]*w31[0] - wm [1]*w31[1] - wm [2]*w31[2] - wm [3]*w31[3];
   Complex_t v13 = wm [0]*wp [0] - wm [1]*wp [1] - wm [2]*wp [2] - wm [3]*wp [3];
   Complex_t v14 = wm [0]*w32[0] - wm [1]*w32[1] - wm [2]*w32[2] - wm [3]*w32[3];
   Complex_t v23 = w31[0]*wp [0] - w31[1]*wp [1] - w31[2]*wp [2] - w31[3]*wp [3];
   Complex_t v24 = w31[0]*w32[0] - w31[1]*w32[1] - w31[2]*w32[2] - w31[3]*w32[3];
   Complex_t v34 = wp [0]*w32[0] - wp [1]*w32[1] - wp [2]*w32[2] - wp [3]*w32[3];

   ANL4DVector q = wm.fP + w31.fP;
   ANL4DVector k = wm.fP + w32.fP;

   Double_t s    = q*q;
   Double_t t    = k*k;
   Double_t mw2  = mw*mw;
   Complex_t dws = -1./Complex_t(s-mw2, TMath::Max(TMath::Sign(mw*gamw,s), 0.));
   Complex_t dwt = -1./Complex_t(t-mw2, TMath::Max(TMath::Sign(mw*gamw,t), 0.));

   Complex_t sv1 =   (w31.fP(0)+q(0))*wm [0] - (w31.fP(1)+q(1))*wm [1]
                   - (w31.fP(2)+q(2))*wm [2] - (w31.fP(3)+q(3))*wm [3];
   Complex_t sv2 = - (wm .fP(0)+q(0))*w31[0] + (wm .fP(1)+q(1))*w31[1]
                   + (wm .fP(2)+q(2))*w31[2] + (wm .fP(3)+q(3))*w31[3];
   Complex_t sv3 =   (w32.fP(0)-q(0))*wp [0] - (w32.fP(1)-q(1))*wp [1]
                   - (w32.fP(2)-q(2))*wp [2] - (w32.fP(3)-q(3))*wp [3];
   Complex_t sv4 = - (wp .fP(0)-q(0))*w32[0] + (wp .fP(1)-q(1))*w32[1]
                   + (wp .fP(2)-q(2))*w32[2] + (wp .fP(3)-q(3))*w32[3];

   Complex_t tv1 =   (w32.fP(0)+k(0))*wm [0] - (w32.fP(1)+k(1))*wm [1]
                   - (w32.fP(2)+k(2))*wm [2] - (w32.fP(3)+k(3))*wm [3];
   Complex_t tv2 = - (wp .fP(0)-k(0))*w31[0] + (wp .fP(1)-k(1))*w31[1]
                   + (wp .fP(2)-k(2))*w31[2] + (wp .fP(3)-k(3))*w31[3];
   Complex_t tv3 =   (w31.fP(0)-k(0))*wp [0] - (w31.fP(1)-k(1))*wp [1]
                   - (w31.fP(2)-k(2))*wp [2] - (w31.fP(3)-k(3))*wp [3];
   Complex_t tv4 = - (wm .fP(0)+k(0))*w32[0] + (wm .fP(1)+k(1))*w32[1]
                   + (wm .fP(2)+k(2))*w32[2] + (wm .fP(3)+k(3))*w32[3];

   TVectorC j12, j34, j14, j32;
   j12[0] = (wm .fP(0)-w31.fP(0))*v12 + sv1*w31[0] + sv2*wm [0];
   j12[1] = (wm .fP(1)-w31.fP(1))*v12 + sv1*w31[1] + sv2*wm [1];
   j12[2] = (wm .fP(2)-w31.fP(2))*v12 + sv1*w31[2] + sv2*wm [2];
   j12[3] = (wm .fP(3)-w31.fP(3))*v12 + sv1*w31[3] + sv2*wm [3];

   j34[0] = (wp .fP(0)-w32.fP(0))*v34 + sv3*w32[0] + sv4*wp [0];
   j34[1] = (wp .fP(1)-w32.fP(1))*v34 + sv3*w32[1] + sv4*wp [1];
   j34[2] = (wp .fP(2)-w32.fP(2))*v34 + sv3*w32[2] + sv4*wp [2];
   j34[3] = (wp .fP(3)-w32.fP(3))*v34 + sv3*w32[3] + sv4*wp [3];

   j14[0] = (wm .fP(0)-w32.fP(0))*v14 + tv1*w32[0] + tv4*wm [0];
   j14[1] = (wm .fP(1)-w32.fP(1))*v14 + tv1*w32[1] + tv4*wm [1];
   j14[2] = (wm .fP(2)-w32.fP(2))*v14 + tv1*w32[2] + tv4*wm [2];
   j14[3] = (wm .fP(3)-w32.fP(3))*v14 + tv1*w32[3] + tv4*wm [3];

   j32[0] = (wp .fP(0)-w31.fP(0))*v23 + tv3*w31[0] + tv2*wp [0];
   j32[1] = (wp .fP(1)-w31.fP(1))*v23 + tv3*w31[1] + tv2*wp [1];
   j32[2] = (wp .fP(2)-w31.fP(2))*v23 + tv3*w31[2] + tv2*wp [2];
   j32[3] = (wp .fP(3)-w31.fP(3))*v23 + tv3*w31[3] + tv2*wp [3];

   Complex_t js12 = q(0)*j12[0] - q(1)*j12[1] - q(2)*j12[2] - q(3)*j12[3];
   Complex_t js34 = q(0)*j34[0] - q(1)*j34[1] - q(2)*j34[2] - q(3)*j34[3];
   Complex_t js14 = k(0)*j14[0] - k(1)*j14[1] - k(2)*j14[2] - k(3)*j14[3];
   Complex_t js32 = k(0)*j32[0] - k(1)*j32[1] - k(2)*j32[2] - k(3)*j32[3];

   Complex_t js = j12[0]*j34[0] - j12[1]*j34[1] - j12[2]*j34[2] - j12[3]*j34[3];
   Complex_t jt = j14[0]*j32[0] - j14[1]*j32[1] - j14[2]*j32[2] - j14[3]*j32[3];

   Complex_t vertex = v12*v34 + v14*v23 - 2.*v13*v24
                    + dws*(js - js12*js34/mw2) + dwt*(jt - js14*js32/mw2);
   *this = vertex * (g31*g32);
}

//----------
// WWWWXX()
//----------
HELVertex::HELVertex(const HELVector  &wm1,
                     const HELVector  &wp1,
                     const HELVector  &wm2,
                     const HELVector  &wp2,
                           Double_t    gwwa,
                           Double_t    gwwz,
                           Double_t    mz,
                           Double_t    gamz)
{
   Complex_t v12 = wm1[0]*wp1[0] - wm1[1]*wp1[1] - wm1[2]*wp1[2] - wm1[3]*wp1[3];
   Complex_t v13 = wm1[0]*wm2[0] - wm1[1]*wm2[1] - wm1[2]*wm2[2] - wm1[3]*wm2[3];
   Complex_t v14 = wm1[0]*wp2[0] - wm1[1]*wp2[1] - wm1[2]*wp2[2] - wm1[3]*wp2[3];
   Complex_t v23 = wp1[0]*wm2[0] - wp1[1]*wm2[1] - wp1[2]*wm2[2] - wp1[3]*wm2[3];
   Complex_t v24 = wp1[0]*wp2[0] - wp1[1]*wp2[1] - wp1[2]*wp2[2] - wp1[3]*wp2[3];
   Complex_t v34 = wm2[0]*wp2[0] - wm2[1]*wp2[1] - wm2[2]*wp2[2] - wm2[3]*wp2[3];

   ANL4DVector q = wm1.fP + wp1.fP;
   ANL4DVector k = wm1.fP + wp2.fP;

   Double_t s    = q*q;
   Double_t t    = k*k;
   Double_t mz2  = mz*mz;
   Complex_t das = -1./s;
   Complex_t dat = -1./t;
   Complex_t dzs = -1./Complex_t(s-mz2, TMath::Max(TMath::Sign(mz*gamz,s), 0.));
   Complex_t dzt = -1./Complex_t(t-mz2, TMath::Max(TMath::Sign(mz*gamz,t), 0.));

   Complex_t sv1 =   (wp1.fP(0)+q(0))*wm1[0] - (wp1.fP(1)+q(1))*wm1[1]
                   - (wp1.fP(2)+q(2))*wm1[2] - (wp1.fP(3)+q(3))*wm1[3];
   Complex_t sv2 = - (wm1.fP(0)+q(0))*wp1[0] + (wm1.fP(1)+q(1))*wp1[1]
                   + (wm1.fP(2)+q(2))*wp1[2] + (wm1.fP(3)+q(3))*wp1[3];
   Complex_t sv3 =   (wp2.fP(0)-q(0))*wm2[0] - (wp2.fP(1)-q(1))*wm2[1]
                   - (wp2.fP(2)-q(2))*wm2[2] - (wp2.fP(3)-q(3))*wm2[3];
   Complex_t sv4 = - (wm2.fP(0)-q(0))*wp2[0] + (wm2.fP(1)-q(1))*wp2[1]
                   + (wm2.fP(2)-q(2))*wp2[2] + (wm2.fP(3)-q(3))*wp2[3];

   Complex_t tv1 =   (wp2.fP(0)+k(0))*wm1[0] - (wp2.fP(1)+k(1))*wm1[1]
                   - (wp2.fP(2)+k(2))*wm1[2] - (wp2.fP(3)+k(3))*wm1[3];
   Complex_t tv2 = - (wm2.fP(0)-k(0))*wp1[0] + (wm2.fP(1)-k(1))*wp1[1]
                   + (wm2.fP(2)-k(2))*wp1[2] + (wm2.fP(3)-k(3))*wp1[3];
   Complex_t tv3 =   (wp1.fP(0)-k(0))*wm2[0] - (wp1.fP(1)-k(1))*wm2[1]
                   - (wp1.fP(2)-k(2))*wm2[2] - (wp1.fP(3)-k(3))*wm2[3];
   Complex_t tv4 = - (wm1.fP(0)+k(0))*wp2[0] + (wm1.fP(1)+k(1))*wp2[1]
                   + (wm1.fP(2)+k(2))*wp2[2] + (wm1.fP(3)+k(3))*wp2[3];

   TVectorC j12, j34, j14, j32;
   j12[0] = (wm1.fP(0)-wp1.fP(0))*v12 + sv1*wp1[0] + sv2*wm1[0];
   j12[1] = (wm1.fP(1)-wp1.fP(1))*v12 + sv1*wp1[1] + sv2*wm1[1];
   j12[2] = (wm1.fP(2)-wp1.fP(2))*v12 + sv1*wp1[2] + sv2*wm1[2];
   j12[3] = (wm1.fP(3)-wp1.fP(3))*v12 + sv1*wp1[3] + sv2*wm1[3];

   j34[0] = (wm2.fP(0)-wp2.fP(0))*v34 + sv3*wp2[0] + sv4*wm2[0];
   j34[1] = (wm2.fP(1)-wp2.fP(1))*v34 + sv3*wp2[1] + sv4*wm2[1];
   j34[2] = (wm2.fP(2)-wp2.fP(2))*v34 + sv3*wp2[2] + sv4*wm2[2];
   j34[3] = (wm2.fP(3)-wp2.fP(3))*v34 + sv3*wp2[3] + sv4*wm2[3];

   j14[0] = (wm1.fP(0)-wp2.fP(0))*v14 + tv1*wp2[0] + tv4*wm1[0];
   j14[1] = (wm1.fP(1)-wp2.fP(1))*v14 + tv1*wp2[1] + tv4*wm1[1];
   j14[2] = (wm1.fP(2)-wp2.fP(2))*v14 + tv1*wp2[2] + tv4*wm1[2];
   j14[3] = (wm1.fP(3)-wp2.fP(3))*v14 + tv1*wp2[3] + tv4*wm1[3];

   j32[0] = (wm2.fP(0)-wp1.fP(0))*v23 + tv3*wp1[0] + tv2*wm2[0];
   j32[1] = (wm2.fP(1)-wp1.fP(1))*v23 + tv3*wp1[1] + tv2*wm2[1];
   j32[2] = (wm2.fP(2)-wp1.fP(2))*v23 + tv3*wp1[2] + tv2*wm2[2];
   j32[3] = (wm2.fP(3)-wp1.fP(3))*v23 + tv3*wp1[3] + tv2*wm2[3];

   Complex_t js12 = q(0)*j12[0] - q(1)*j12[1] - q(2)*j12[2] - q(3)*j12[3];
   Complex_t js34 = q(0)*j34[0] - q(1)*j34[1] - q(2)*j34[2] - q(3)*j34[3];
   Complex_t js14 = k(0)*j14[0] - k(1)*j14[1] - k(2)*j14[2] - k(3)*j14[3];
   Complex_t js32 = k(0)*j32[0] - k(1)*j32[1] - k(2)*j32[2] - k(3)*j32[3];

   Complex_t js = j12[0]*j34[0] - j12[1]*j34[1] - j12[2]*j34[2] - j12[3]*j34[3];
   Complex_t jt = j14[0]*j32[0] - j14[1]*j32[1] - j14[2]*j32[2] - j14[3]*j32[3];

   Double_t gwwa2 = gwwa*gwwa;
   Double_t gwwz2 = gwwz*gwwz;
   Double_t dgw2  = gwwa2 + gwwz2;
   Complex_t vertex = (v12*v34 + v14*v23 - 2.*v13*v24)*dgw2
                    + (dzs*gwwz2 + das*gwwa2)*js - dzs*gwwz2*js12*js34/mz2
                    + (dzt*gwwz2 + dat*gwwa2)*jt - dzt*gwwz2*js14*js32/mz2;

   *this = - vertex;
}


//----------
// GGGGXX()
//----------
HELVertex::HELVertex(const HELVector  &wm,
                     const HELVector  &w31,
                     const HELVector  &wp,
                     const HELVector  &w32,
                           Double_t     g)
{
   Complex_t v12 = wm [0]*w31[0] - wm [1]*w31[1] - wm [2]*w31[2] - wm [3]*w31[3];
   Complex_t v13 = wm [0]*wp [0] - wm [1]*wp [1] - wm [2]*wp [2] - wm [3]*wp [3];
   Complex_t v14 = wm [0]*w32[0] - wm [1]*w32[1] - wm [2]*w32[2] - wm [3]*w32[3];
   Complex_t v23 = w31[0]*wp [0] - w31[1]*wp [1] - w31[2]*wp [2] - w31[3]*wp [3];
   Complex_t v24 = w31[0]*w32[0] - w31[1]*w32[1] - w31[2]*w32[2] - w31[3]*w32[3];
   Complex_t v34 = wp [0]*w32[0] - wp [1]*w32[1] - wp [2]*w32[2] - wp [3]*w32[3];

   ANL4DVector q = wm.fP + w31.fP;

   Double_t s    = q*q;
   Complex_t dws = -1./Complex_t(s, 0.);

   Complex_t sv1 =   (w31.fP(0)+q(0))*wm [0] - (w31.fP(1)+q(1))*wm [1]
                   - (w31.fP(2)+q(2))*wm [2] - (w31.fP(3)+q(3))*wm [3];
   Complex_t sv2 = - (wm .fP(0)+q(0))*w31[0] + (wm .fP(1)+q(1))*w31[1]
                   + (wm .fP(2)+q(2))*w31[2] + (wm .fP(3)+q(3))*w31[3];
   Complex_t sv3 =   (w32.fP(0)-q(0))*wp [0] - (w32.fP(1)-q(1))*wp [1]
                   - (w32.fP(2)-q(2))*wp [2] - (w32.fP(3)-q(3))*wp [3];
   Complex_t sv4 = - (wp .fP(0)-q(0))*w32[0] + (wp .fP(1)-q(1))*w32[1]
                   + (wp .fP(2)-q(2))*w32[2] + (wp .fP(3)-q(3))*w32[3];

   TVectorC j12, j34;
   j12[0] = (wm .fP(0)-w31.fP(0))*v12 + sv1*w31[0] + sv2*wm [0];
   j12[1] = (wm .fP(1)-w31.fP(1))*v12 + sv1*w31[1] + sv2*wm [1];
   j12[2] = (wm .fP(2)-w31.fP(2))*v12 + sv1*w31[2] + sv2*wm [2];
   j12[3] = (wm .fP(3)-w31.fP(3))*v12 + sv1*w31[3] + sv2*wm [3];

   j34[0] = (wp .fP(0)-w32.fP(0))*v34 + sv3*w32[0] + sv4*wp [0];
   j34[1] = (wp .fP(1)-w32.fP(1))*v34 + sv3*w32[1] + sv4*wp [1];
   j34[2] = (wp .fP(2)-w32.fP(2))*v34 + sv3*w32[2] + sv4*wp [2];
   j34[3] = (wp .fP(3)-w32.fP(3))*v34 + sv3*w32[3] + sv4*wp [3];

   Complex_t js = j12[0]*j34[0] - j12[1]*j34[1] - j12[2]*j34[2] - j12[3]*j34[3];

   Complex_t vertex = v14*v23 - v13*v24 + dws*js;
   *this = vertex * (g*g);
}


//----------
// VSSXXX()
//----------
HELVertex::HELVertex(const HELVector &vc,
                     const HELScalar &s1,
                     const HELScalar &s2,
                           Double_t    g)
{
   ANL4DVector p = s1.fP - s2.fP;
   *this = g*s1*s2*(vc[0]*p(0) - vc[1]*p(1) - vc[2]*p(2) - vc[3]*p(3));
}

//----------
// VVSSXX()
//----------
HELVertex::HELVertex(const HELVector &v1,
                     const HELVector &v2,
                     const HELScalar &s1,
                     const HELScalar &s2,
                           Double_t    g)
{
   *this = g*s1*s2*(v1[0]*v2[0] - v1[1]*v2[1] - v1[2]*v2[2] - v1[3]*v2[3]);
}

//----------
// IOSXXX()
//----------
HELVertex::HELVertex(const HELFermion &fin,
                     const HELFermion &fout,
                     const HELScalar  &sc,
                           Complex_t   gl,
                           Complex_t   gr)
{
   *this = sc * (gl*(fout[0]*fin[0]+fout[1]*fin[1]) 
               + gr*(fout[2]*fin[2]+fout[3]*fin[3])) ;
}
