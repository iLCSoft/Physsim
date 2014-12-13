// 2013/04/25 Alpha and Sin2W changed to the values used in Whizard
//      --- Junping Tian & Keisuke Fujii
// 2014/11/25 AlphaS changed to value in Whizard (0.11780)
//            

#ifndef GENNUMCON_H
#define GENNUMCON_H
#include "TMath.h"
//*----------------------------------------------------------------------
//*     Numerical and Natural Constants
//*----------------------------------------------------------------------
static const Double_t   kSqh     = TMath::Sqrt(0.5); // sqrt(1/2)
static const Double_t   kPi      = TMath::Pi();      // Pi
static const Double_t   k2Pi     = 2*TMath::Pi();    // 2*Pi
static const Double_t   k4Pi     = 4*TMath::Pi();    // 4*Pi
static const Double_t   k8Pi     = 8*TMath::Pi();    // 8*Pi
static const Double_t   kGeV2fb  = 0.389379292e12;   // GeV to fb

//static const Double_t   kAlpha   = 1./132.50495;     // alpha(mz)
static const Double_t   kAlpha   = 1./128.;          // alpha(mz)
static const Double_t   kAlpha0  = 1./137.0359895;   // alpha(q=0) = 1/137.
static const Double_t   kAlphaS  = 0.11780;             // alphaS(mz) = 0.12
//static const Double_t   kAlphaS  = 0.12;             // alphaS(mz) = 0.12
//static const Double_t   kSin2W   = 0.222249945;      // sin^2(theta_W)
static const Double_t   kSin2W   = 0.230;            // sin^2(theta_W)
static const Double_t   kSinW    = TMath::Sqrt(kSin2W);       // sin(theta_W)
static const Double_t   kCos2W   = (1. - kSinW)*(1. + kSinW); // cos^2(theta_W)
static const Double_t   kCosW    = TMath::Sqrt(kCos2W);       // cos(theta_W)
static const Double_t   kSinCosW = kSinW*kCosW;               // sin(2theta_W)/2
static const Double_t   kGe      = TMath::Sqrt(k4Pi*kAlpha);  // e
static const Double_t   kGw      = kGe/kSinW;                 // gw
static const Double_t   kGz      = kGw/kCosW;                 // gz
  
static const Double_t   kM_e     = 0.510998902e-3;   // electron mass [GeV]
static const Double_t   kM_z     = 91.188;           // Z mass [GeV]
static const Double_t   kM_w     = kM_z*kCosW;       // W mass [GeV]
static const Char_t    *kName[2][2][3] = {{{"nu_e", "nu_mu"  , "nu_tau"},
                                           {"e"   , "mu"     , "tau"   }},
                                          {{"up"  , "charm"  , "top"   },
                                           {"down", "strange", "bottom"}}};
static const Int_t      kPID [2][2][3] = {{{    12,        14,       16},
                                           {    11,        13,       15}},
                                          {{     2,         4,        6},
                                           {     1,         3,        5}}};
static const Double_t   kChrg[2][2][3] = {{{    0.,        0.,       0.},
                                           {   -1.,       -1.,      -1.}},
                                          {{  2/3.,      2/3.,     2/3.},
                                           { -1/3.,     -1/3.,    -1/3.}}};

static const Double_t   kMass[2][2][3] = {{{0.000000, 0.00000,   0.0000},
                                           {0.511e-3, 0.10566,   1.7770}},
                                          {{0.04    , 1.5    , 175.    },
                                           {0.04    , 0.1    ,   4.7   }}};
static const Double_t   kVkm [3][3]    = {{0.975, 0.222, 0.010},
                                          {0.222, 0.974, 0.043},
                                          {0.010, 0.043, 0.999}};
#endif
