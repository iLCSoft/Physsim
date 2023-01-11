// $Id: LinkDef.h,v 1.3 2007/01/18 15:51:30 fujiik Exp $

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ANL4DVector;

#pragma link C++ class GENBranch+;
#pragma link C++ class GENFrame+;
#pragma link C++ class GENPhase2+;
#pragma link C++ class GENDecayMode+;
#pragma link C++ class GENModePicker+;
#pragma link C++ class GENPDTEntry+;
#pragma link C++ class GENPDTWBoson+;
#pragma link C++ class GENPDTZBoson+;
#pragma link C++ class GENPDTPhoton+;
#pragma link C++ class GENPDTGluon+;

#pragma link C++ class TVectorC+;
#pragma link C++ class HELFermion+;
#pragma link C++ class HELVector+;
#pragma link C++ class HELScalar+;
#pragma link C++ class HELVertex+;

#pragma link C++ class lcme::LCMEBase+;
#pragma link C++ class lcme::LCMEZHH+;
#pragma link C++ class lcme::LCMEZH+;
#pragma link C++ class lcme::LCMENNHH+;
#pragma link C++ class lcme::LCMENNH+;
#pragma link C++ class lcme::LCMEEEH+;
#pragma link C++ class lcme::LCMEEEZ+;
#pragma link C++ class lcme::LCMEZZ+;
#pragma link C++ class lcme::LCMEZZH+;
#pragma link C++ class lcme::LCMEZZZ+;
#pragma link C++ class lcme::LCMEWW+;

#pragma link C++ class lcme::TAttLockable+;

#endif
