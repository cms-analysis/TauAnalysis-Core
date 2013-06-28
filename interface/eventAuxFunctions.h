#ifndef TauAnalysis_Core_eventAuxFunctions_h
#define TauAnalysis_Core_eventAuxFunctions_h

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

template<typename T>
bool isValidRef(const edm::Ref<T>& ref)
{
  return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );  
}

template<typename T>
bool isValidRef(const edm::Ptr<T>& ptr)
{
  return ( (ptr.isAvailable() || ptr.isTransient()) && ptr.isNonnull() );  
}

#endif
