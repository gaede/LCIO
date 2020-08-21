#ifndef __CINT__ 
#include "lcio.h"

#include "MT/Types.h"
#include "MT/LCReader.h"
#include "MT/LCReaderListener.h"
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCEventImpl.h"
#include "UTIL/LCIterator.h"
#endif

#include <cstdlib>
#include <mutex>
#include <future>
#include <functional>
#include <thread>
#include <unistd.h>

static std::vector<std::string> FILEN ; 

using namespace std ;
using namespace lcio ;
using LCEventPtr = MT::LCEventPtr;
using LCRunHeaderPtr = MT::LCRunHeaderPtr;

std::mutex printMutex;
#define SAFE_PRINT( message ) { std::lock_guard<std::mutex> lock(printMutex); std::cout << message << std::endl; }
#define SAFE_CODE( code ) { std::lock_guard<std::mutex> lock(printMutex); code; }

template<class T>
double inv_mass(T* p1, T* p2){
  double e = p1->getEnergy()+p2->getEnergy() ;
  double px = p1->getMomentum()[0]+p2->getMomentum()[0];
  double py = p1->getMomentum()[1]+p2->getMomentum()[1];
  double pz = p1->getMomentum()[2]+p2->getMomentum()[2];
  return( sqrt( e*e - px*px - py*py - pz*pz  ) );
}

template<class T>
TLorentzVector v4(T* p){
  return TLorentzVector( p->getMomentum()[0],p->getMomentum()[1], p->getMomentum()[2],p->getEnergy());
}


//======================================================
enum histoindex{
  muonmass, jetmass, recoilmass,
  Size
};
//======================================================

template<unsigned SIZE>
class Histograms{
public:
  Histograms() = default ;

  void create(int idx, const char* n, int nBin=100, double min=0., double max=0. ){
    create( idx , n , n , nBin, min , max ) ;
  }

  void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0. ){

    _h[idx] = new TH1D( n, t , nBin , min, max ) ;

    std::cout << " create histo " <<  n << " at index " << idx << std::endl ;
  }

  void create(int idx, const char* n, const char* t,  int nBin , double* bins ){

    _h[idx] = new TH1D( n, t , nBin , bins ) ;

    std::cout  << " create histo " <<  n << " at index " << idx << std::endl ;
  }

  void fill( int idx , double val, double weight=1.0 ){
    std::lock_guard<std::mutex> lock( _m[ idx ] );
    _h[idx]->Fill( val , weight ) ;
  }

  TH1* operator[](unsigned idx) { return _h[idx] ; }

protected:

  std::array<TH1*,SIZE> _h ;
  std::array<std::mutex,SIZE> _m ;
};
//======================================================



//======================================================
class Scheduler final : public MT::LCReaderListener {
private:
  typedef double                             task_return_type ;
  typedef std::future<task_return_type>      future_type ;
  typedef std::vector<future_type>           future_list ;
  
public:
  Scheduler( unsigned int maxNTasks, Histograms<Size>* h ) :
    _maxNTasks(maxNTasks),
    _h(h) {
    SAFE_PRINT( "Scheduler created with maxNTasks = " << _maxNTasks ) ;
  }
  
  ~Scheduler() {
    waitForAll();
  }
  
  void startTask( LCEventPtr event ) {
    while ( not canStartNewTask() ) {
      processFinishedTasks() ;
      usleep(1) ;
    }
//    SAFE_PRINT( "Starting new task ..." ) ;

    _futures.push_back( std::async( std::launch::async,
				    [&](LCEventPtr evt) {
      

// ===============================================							  


   LCIterator<ReconstructedParticle> jets( evt.get(), "Jets" ) ;
   LCIterator<ReconstructedParticle> muons( evt.get(), "IsolatedMuons" ) ;
   
   if( jets.size() != 2)
     return 0.;
   
   if( muons.size() != 2)
     return 0.;
   
   auto mu1 = muons.next(); 
   auto mu2 = muons.next(); 
   _h->fill( muonmass, inv_mass( mu1, mu2) ) ;
   
   auto j1 = jets.next(); 
   auto j2 = jets.next(); 
   _h->fill( jetmass, inv_mass( j1, j2) ) ;
   
   
   // the recoil mass
   const auto& vm1 = v4(mu1) ;
   const auto& vm2 = v4(mu2) ;
   TLorentzVector ecms(0.,0.,0.,250.) ;
   TLorentzVector recoil = ecms - ( vm1 + vm2 ) ;
   _h->fill(recoilmass, recoil.M() ) ;
   

   return 0. ;
// ===============================================							  
				    }, event)) ;
  }
  
  void waitForAll() {
    SAFE_PRINT( "waitForAll()" ) ;
    for ( unsigned int i=0 ; i<_futures.size() ; ++i ) {
      ++_n ;
      _e += _futures.at(i).get();
    }
    _futures.clear();
  }
  void printResult(){

    std::cout << " Finished processing of " << _n << " events with everage energy = " << _e/_n << std::endl ;
  }
  
private:
  
  void processFinishedTasks() {
    auto iter = _futures.begin() ;
    while ( iter != _futures.end() ) {
      auto status = iter->wait_for( std::chrono::seconds(0) ) ;
      if ( status == std::future_status::ready ) {
        // process finished task and remove it for the pending task list
        _e += iter->get();
	++_n ;
        iter = _futures.erase( iter );
      }
      else {
        iter++ ;        
      }
    }
  }

  
  bool canStartNewTask() const {
    return ( _futures.size() < _maxNTasks ) ;
  }
  
  void processEvent( LCEventPtr event ) override {
    startTask( event ) ;
  }
  
  void processRunHeader( LCRunHeaderPtr hdr ) override {
    // Wait for all event task to finish
    // and then process the new run header
    waitForAll();
    UTIL::LCTOOLS::dumpRunHeader( hdr.get() );
  }
  
private:
  unsigned int                   _maxNTasks {} ;
  future_list                    _futures {} ;
  double _e = 0.;
  int _n=0;
  Histograms<Size>* _h ;
};

/** Small utility to dump events in parallel from different files
 */
int parallel_processing(char* FILEN ){


    std::vector<std::string> inputFiles ;
    inputFiles.push_back( FILEN )  ;
    
    // The LCReader to read events and runs
    MT::LCReader reader( MT::LCReader::directAccess  | MT::LCReader::lazyUnpack ) ;
    reader.setReadCollectionNames( { "PFOs","Jets","IsolatedMuons" } );
    reader.open( inputFiles );
    
    unsigned int nThreads = std::thread::hardware_concurrency() ;
    char *nthreadsenv = getenv( "LCIO_MAX_THREADS" ) ;
    if ( nthreadsenv ) {
      nThreads = atoi( nthreadsenv );
      if ( nThreads <= 0 ) {
        nThreads = 1;
      }
    }
    //--------------------------------------------------------------

    Histograms<Size> h ;
    h.create( muonmass ,"inv. mass - muons", 100,  60. , 120. );
    h.create( jetmass ,"inv. mass - jets", 100, 0. , 150. );
    h.create( recoilmass ,"recoil mass", 100, 110. , 170. );
    
    //------------------------------------------------------------
    // The task scheduler
    Scheduler scheduler(nThreads,&h);
    MT::LCReaderListenerList listeners ;
    listeners.insert( &scheduler ) ;
    
    // Read stream and process events and run headers
    while (1) {
      try {
        reader.readNextRecord( listeners );        
      }
      catch ( const IO::EndOfDataException& ) {
        break;
      }
    }

    scheduler.printResult() ; 

  //===================================================================================================
    TCanvas* c1 = new TCanvas("recoil plots");
    
    c1->Divide(2,2);
    c1->cd(1) ;
    h[muonmass]->Draw() ;
    c1->cd(2) ;
    h[jetmass]->Draw();
    c1->cd(3) ;
    h[recoilmass]->Draw();
    
    c1->Print("recoil_plots_MT.pdf") ;
 //===================================================================================================

 return 0 ;
}
