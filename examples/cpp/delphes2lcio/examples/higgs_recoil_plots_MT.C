#ifndef __CINT__ 
#include "lcio.h"
#include "MT/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCEvent.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/LCIterator.h"
#endif

#include "TH1F.h"
#include "TLorentzVector.h"
#include "tsqueue.h"

#include <thread>
#include <future>

/*
put this into your .rootlogon.C file

{
 gInterpreter->AddIncludePath("$LCIO");
 gSystem->Load("${LCIO}/lib/liblcio.so");
 gSystem->Load("${LCIO}/lib/liblcioDict.so");
}

for the LCIO API documentation see:
http://lcio.desy.de/v02-09/doc/doxygen_api/html/index.html
*/


using namespace lcio ;

const UInt_t poolSize = 3U;

std::mutex printMutex;
#define SAFE_PRINT( message ) { std::lock_guard<std::mutex> lock(printMutex); std::cout << message << std::endl; }


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

/** Example script for creating invariant mass plots from a ee -> HZ-> X mu mu sample
 *
 */
 
void higgs_recoil_plots_MT(const char* FILEN) {

 int nEvents = 0  ;
 int maxEvt = 100000 ;  // change as needed

//----
 TH1::AddDirectory(false);
//---------------------------------
 

 MT::LCReader* lcReader = new MT::LCReader( MT::LCReader::directAccess
					     | MT::LCReader::lazyUnpack ) ;
 lcReader->setReadCollectionNames( { "PFOs","Jets","IsolatedMuons" } );
 //lcReader->setReadCollectionNames( {"MCParticle", "PFOs", "RecoMCTruthLink" ,"MCTruthRecoLink","Jets","IsolatedMuons" } ) ;
 lcReader->open( FILEN ) ;
  

 Tsqueue< EVENT::LCEvent* > qu(10) ;

 // helper function to populate event queue: 
 auto fillqueue = [&](void) {
		       std::unique_ptr<EVENT::LCEvent> evt ;
		       while( ( evt = lcReader->readNextEvent()) != 0 && nEvents++ < maxEvt  ){
			 qu.push( evt.release() ) ;
		       }
		       qu.close() ;		       
		  } ;
 
 
 //==================== the event loop function ============================================================

   auto fillHistos =
     [&](int seed=0) {
       
       auto* histos = new TObjArray ;
       
       histos->Add( new TH1F("hmuonmass","inv. mass - muons", 100,  60. , 120. ) );
       histos->Add( new TH1F("hjetmass","inv. mass - jets", 100, 0. , 150. ) );
       histos->Add( new TH1F("hrecoilm","recoil mass", 100, 110. , 170. ) );

       
       SAFE_PRINT( " started new thread w/ seed: " << seed  )

       TH1F* hmuonmass = (TH1F* )histos->At(0) ;
       TH1F* hjetmass =  (TH1F* )histos->At(1) ;
       TH1F* hrecoilm =  (TH1F* )histos->At(2) ;

       while( auto eo = qu.pop() ) {
	 
	 std::unique_ptr<EVENT::LCEvent> evt( *eo ) ;
   
//	 SAFE_PRINT( " --- " << seed << " processing evt " << evt->getEventNumber()  )


	 LCIterator<ReconstructedParticle> jets( evt.get(), "Jets" ) ;
	 LCIterator<ReconstructedParticle> muons( evt.get(), "IsolatedMuons" ) ;

	 if( jets.size() != 2)
	   continue;

	 if( muons.size() != 2)
	   continue;
	 
	 auto mu1 = muons.next(); 
	 auto mu2 = muons.next(); 
	 hmuonmass->Fill( inv_mass( mu1, mu2) ) ;
	 
	 auto j1 = jets.next(); 
	 auto j2 = jets.next(); 
	 hjetmass->Fill( inv_mass( j1, j2) ) ;
	 
	 
	 // the recoil mass
	 const auto& vm1 = v4(mu1) ;
	 const auto& vm2 = v4(mu2) ;
	 TLorentzVector ecms(0.,0.,0.,250.) ;
	 TLorentzVector recoil = ecms - ( vm1 + vm2 ) ;
	 hrecoilm->Fill( recoil.M() ) ;
       }

       return histos;
     };
			 
   
  //===============================================================================================
   typedef TObjArray*                         task_return_type ;
   typedef std::future<task_return_type>      future_type ;
   typedef std::vector<future_type>           future_list ;
   typedef std::vector<task_return_type>      return_vec ;

   future_list fuhist ;

   // start the reader thread
   std::thread readerthread( fillqueue ) ; 

   // start the histogramming threads
   for(unsigned i=0 ; i < poolSize ; ++i) {
     fuhist.push_back(  std::async( std::launch::async, fillHistos , i + 1 ) );
   }

   readerthread.join() ;
   
   return_vec  histvec ;
   
   // collect the resulting histograms:
   for(auto& fu : fuhist )
      histvec.push_back( fu.get() ) ;
   
   fuhist.clear() ;
   
   
   //---
   std::cout << " done processing - now merge histograms " << std::endl ;
   //---
   
   TObjArray* histos0 = histvec[0] ;
   
   for(int j=0;j<histos0->GetEntries() ;++j){
     
     TH1F* h0 = (TH1F* )histos0->At(j) ;
     
     TList col ;
      for(unsigned i=1 ; i < histvec.size() ; ++i){
	
        TObjArray* histos2 = histvec[i] ;
        TH1F* h1 = (TH1F* ) histos2->At(j) ;
        col.Add( h1 ) ;
      }
      h0->Merge( &col ) ;
      //     col.Delete() ;
   }
 //===================================================================================================
   
   TH1F* hmuonmass = (TH1F* )histos0->At(0) ;
   TH1F* hjetmass =  (TH1F* )histos0->At(1) ;
   TH1F* hrecoilm =  (TH1F* )histos0->At(2) ;
   
 //===================================================================================================
   TCanvas* c1 = new TCanvas("recoil plots");
   
   c1->Divide(2,2);
   c1->cd(1) ;
   hmuonmass->Draw() ;
   c1->cd(2) ;
   hjetmass->Draw();
   c1->cd(3) ;
   hrecoilm->Draw();
   
   c1->Print("recoil_plots.pdf") ;
 //===================================================================================================
   
  std::cout << std::endl 
	    <<  "  " <<  nEvents 
	    << " events read from file: " 
	    << FILEN << std::endl  ;
  
  
  lcReader->close() ;
  delete lcReader ;
}
