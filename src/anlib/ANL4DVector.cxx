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
//* $Id: ANL4DVector.cxx,v 1.1.1.1 2005/02/15 02:14:56 fujiik Exp $
//*************************************************************************
//
#include "ANL4DVector.h"
//_____________________________________________________________________
//  ----------------
//  Lockable LVector
//  ----------------
//

ClassImp(ANL4DVector)
