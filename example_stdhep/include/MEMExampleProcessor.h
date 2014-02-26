#ifndef MEMExampleProcessor_h
#define MEMExampleProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "physsim/LCMEZHH.h"

/**  Example processor to use Matrix Element Method.
 * 
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author Junping Tian, KEK
 * @version $Id: MEMExampleProcessor.h,v 1.0 2013-10-2$ 
 */

class MEMExampleProcessor : public marlin::Processor {
  
 public:
  
  virtual marlin::Processor*  newProcessor() { return new MEMExampleProcessor ; }
  
  
  MEMExampleProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colMCP ;
  std::string _colMCTL ;
  std::string _colPFOs ;

  int _nRun ;
  int _nEvt ;

  lcme::LCMEZHH *_zhh;

} ;

#endif



