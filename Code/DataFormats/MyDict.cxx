// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MyDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "SVFitObject.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *SVFitObject_Dictionary();
   static void SVFitObject_TClassManip(TClass*);
   static void *new_SVFitObject(void *p = 0);
   static void *newArray_SVFitObject(Long_t size, void *p);
   static void delete_SVFitObject(void *p);
   static void deleteArray_SVFitObject(void *p);
   static void destruct_SVFitObject(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SVFitObject*)
   {
      ::SVFitObject *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SVFitObject));
      static ::ROOT::TGenericClassInfo 
         instance("SVFitObject", "SVFitObject.h", 13,
                  typeid(::SVFitObject), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SVFitObject_Dictionary, isa_proxy, 4,
                  sizeof(::SVFitObject) );
      instance.SetNew(&new_SVFitObject);
      instance.SetNewArray(&newArray_SVFitObject);
      instance.SetDelete(&delete_SVFitObject);
      instance.SetDeleteArray(&deleteArray_SVFitObject);
      instance.SetDestructor(&destruct_SVFitObject);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SVFitObject*)
   {
      return GenerateInitInstanceLocal((::SVFitObject*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SVFitObject*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SVFitObject_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SVFitObject*)0x0)->GetClass();
      SVFitObject_TClassManip(theClass);
   return theClass;
   }

   static void SVFitObject_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_SVFitObject(void *p) {
      return  p ? new(p) ::SVFitObject : new ::SVFitObject;
   }
   static void *newArray_SVFitObject(Long_t nElements, void *p) {
      return p ? new(p) ::SVFitObject[nElements] : new ::SVFitObject[nElements];
   }
   // Wrapper around operator delete
   static void delete_SVFitObject(void *p) {
      delete ((::SVFitObject*)p);
   }
   static void deleteArray_SVFitObject(void *p) {
      delete [] ((::SVFitObject*)p);
   }
   static void destruct_SVFitObject(void *p) {
      typedef ::SVFitObject current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SVFitObject

namespace ROOT {
   static TClass *vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_Dictionary();
   static void vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p = 0);
   static void *newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p);
   static void deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p);
   static void destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >*)
   {
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", -2, "vector", 210,
                  typeid(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >) );
      instance.SetNew(&new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >()));

      ::ROOT::AddClassAlternate("vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >, std::allocator<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >*)0x0)->GetClass();
      vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > : new vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >;
   }
   static void *newArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >[nElements] : new vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p) {
      delete ((vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >*)p);
   }
   static void deleteArray_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p) {
      delete [] ((vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >*)p);
   }
   static void destruct_vectorlEROOTcLcLMathcLcLLorentzVectorlEROOTcLcLMathcLcLPxPyPzE4DlEdoublegRsPgRsPgR(void *p) {
      typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >

namespace {
  void TriggerDictionaryInitialization_MyDict_Impl() {
    static const char* headers[] = {
"SVFitObject.h",
0
    };
    static const char* includePaths[] = {
"/grid_mnt/libcern/root/6.24.06/centos7.6-x86_64/include/",
"/grid_mnt/home-pbs/msessini/Test/workdirNewFormat18_oldGEF_Aug_17_2022/Code/DataFormats/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MyDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$SVFitObject.h")))  SVFitObject;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MyDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "SVFitObject.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"SVFitObject", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MyDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MyDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MyDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MyDict() {
  TriggerDictionaryInitialization_MyDict_Impl();
}
