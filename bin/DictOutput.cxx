// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME bindIDictOutput
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/DataReader.h"
#include "include/ExperimentalInput.h"
#include "include/GOSIAReader.h"
#include "include/GOSIASimFitter.h"
#include "include/GOSIASimMinFCN.h"
#include "include/Literature.h"
#include "include/MatrixElement.h"
#include "include/MiscFunctions.h"
#include "include/NormalizationSampler.h"
#include "include/Nucleus.h"
#include "include/NucleusReader.h"
#include "include/ScalingFitFCN.h"
#include "include/ScalingParameter.h"
#include "include/TransitionRates.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *Nucleus_Dictionary();
   static void Nucleus_TClassManip(TClass*);
   static void *new_Nucleus(void *p = 0);
   static void *newArray_Nucleus(Long_t size, void *p);
   static void delete_Nucleus(void *p);
   static void deleteArray_Nucleus(void *p);
   static void destruct_Nucleus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Nucleus*)
   {
      ::Nucleus *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Nucleus));
      static ::ROOT::TGenericClassInfo 
         instance("Nucleus", "include/Nucleus.h", 22,
                  typeid(::Nucleus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Nucleus_Dictionary, isa_proxy, 4,
                  sizeof(::Nucleus) );
      instance.SetNew(&new_Nucleus);
      instance.SetNewArray(&newArray_Nucleus);
      instance.SetDelete(&delete_Nucleus);
      instance.SetDeleteArray(&deleteArray_Nucleus);
      instance.SetDestructor(&destruct_Nucleus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Nucleus*)
   {
      return GenerateInitInstanceLocal((::Nucleus*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Nucleus*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Nucleus_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Nucleus*)0x0)->GetClass();
      Nucleus_TClassManip(theClass);
   return theClass;
   }

   static void Nucleus_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/Nucleus.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExptData_Dictionary();
   static void ExptData_TClassManip(TClass*);
   static void *new_ExptData(void *p = 0);
   static void *newArray_ExptData(Long_t size, void *p);
   static void delete_ExptData(void *p);
   static void deleteArray_ExptData(void *p);
   static void destruct_ExptData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExptData*)
   {
      ::ExptData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExptData));
      static ::ROOT::TGenericClassInfo 
         instance("ExptData", "include/ExperimentalInput.h", 20,
                  typeid(::ExptData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExptData_Dictionary, isa_proxy, 4,
                  sizeof(::ExptData) );
      instance.SetNew(&new_ExptData);
      instance.SetNewArray(&newArray_ExptData);
      instance.SetDelete(&delete_ExptData);
      instance.SetDeleteArray(&deleteArray_ExptData);
      instance.SetDestructor(&destruct_ExptData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExptData*)
   {
      return GenerateInitInstanceLocal((::ExptData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ExptData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExptData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExptData*)0x0)->GetClass();
      ExptData_TClassManip(theClass);
   return theClass;
   }

   static void ExptData_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/ExperimentalInput.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ExperimentData_Dictionary();
   static void ExperimentData_TClassManip(TClass*);
   static void *new_ExperimentData(void *p = 0);
   static void *newArray_ExperimentData(Long_t size, void *p);
   static void delete_ExperimentData(void *p);
   static void deleteArray_ExperimentData(void *p);
   static void destruct_ExperimentData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExperimentData*)
   {
      ::ExperimentData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ExperimentData));
      static ::ROOT::TGenericClassInfo 
         instance("ExperimentData", "include/ExperimentalInput.h", 94,
                  typeid(::ExperimentData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ExperimentData_Dictionary, isa_proxy, 4,
                  sizeof(::ExperimentData) );
      instance.SetNew(&new_ExperimentData);
      instance.SetNewArray(&newArray_ExperimentData);
      instance.SetDelete(&delete_ExperimentData);
      instance.SetDeleteArray(&deleteArray_ExperimentData);
      instance.SetDestructor(&destruct_ExperimentData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExperimentData*)
   {
      return GenerateInitInstanceLocal((::ExperimentData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ExperimentData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ExperimentData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ExperimentData*)0x0)->GetClass();
      ExperimentData_TClassManip(theClass);
   return theClass;
   }

   static void ExperimentData_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/ExperimentalInput.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DataReader_Dictionary();
   static void DataReader_TClassManip(TClass*);
   static void delete_DataReader(void *p);
   static void deleteArray_DataReader(void *p);
   static void destruct_DataReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DataReader*)
   {
      ::DataReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DataReader));
      static ::ROOT::TGenericClassInfo 
         instance("DataReader", "include/DataReader.h", 19,
                  typeid(::DataReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DataReader_Dictionary, isa_proxy, 4,
                  sizeof(::DataReader) );
      instance.SetDelete(&delete_DataReader);
      instance.SetDeleteArray(&deleteArray_DataReader);
      instance.SetDestructor(&destruct_DataReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DataReader*)
   {
      return GenerateInitInstanceLocal((::DataReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DataReader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DataReader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DataReader*)0x0)->GetClass();
      DataReader_TClassManip(theClass);
   return theClass;
   }

   static void DataReader_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/DataReader.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *GOSIAReader_Dictionary();
   static void GOSIAReader_TClassManip(TClass*);
   static void delete_GOSIAReader(void *p);
   static void deleteArray_GOSIAReader(void *p);
   static void destruct_GOSIAReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GOSIAReader*)
   {
      ::GOSIAReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::GOSIAReader));
      static ::ROOT::TGenericClassInfo 
         instance("GOSIAReader", "include/GOSIAReader.h", 13,
                  typeid(::GOSIAReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &GOSIAReader_Dictionary, isa_proxy, 4,
                  sizeof(::GOSIAReader) );
      instance.SetDelete(&delete_GOSIAReader);
      instance.SetDeleteArray(&deleteArray_GOSIAReader);
      instance.SetDestructor(&destruct_GOSIAReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GOSIAReader*)
   {
      return GenerateInitInstanceLocal((::GOSIAReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GOSIAReader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *GOSIAReader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::GOSIAReader*)0x0)->GetClass();
      GOSIAReader_TClassManip(theClass);
   return theClass;
   }

   static void GOSIAReader_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/GOSIAReader.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ScalingParameter_Dictionary();
   static void ScalingParameter_TClassManip(TClass*);
   static void *new_ScalingParameter(void *p = 0);
   static void *newArray_ScalingParameter(Long_t size, void *p);
   static void delete_ScalingParameter(void *p);
   static void deleteArray_ScalingParameter(void *p);
   static void destruct_ScalingParameter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ScalingParameter*)
   {
      ::ScalingParameter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ScalingParameter));
      static ::ROOT::TGenericClassInfo 
         instance("ScalingParameter", "include/ScalingParameter.h", 25,
                  typeid(::ScalingParameter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ScalingParameter_Dictionary, isa_proxy, 4,
                  sizeof(::ScalingParameter) );
      instance.SetNew(&new_ScalingParameter);
      instance.SetNewArray(&newArray_ScalingParameter);
      instance.SetDelete(&delete_ScalingParameter);
      instance.SetDeleteArray(&deleteArray_ScalingParameter);
      instance.SetDestructor(&destruct_ScalingParameter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ScalingParameter*)
   {
      return GenerateInitInstanceLocal((::ScalingParameter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ScalingParameter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ScalingParameter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ScalingParameter*)0x0)->GetClass();
      ScalingParameter_TClassManip(theClass);
   return theClass;
   }

   static void ScalingParameter_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/ScalingParameter.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LitLifetime_Dictionary();
   static void LitLifetime_TClassManip(TClass*);
   static void *new_LitLifetime(void *p = 0);
   static void *newArray_LitLifetime(Long_t size, void *p);
   static void delete_LitLifetime(void *p);
   static void deleteArray_LitLifetime(void *p);
   static void destruct_LitLifetime(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LitLifetime*)
   {
      ::LitLifetime *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LitLifetime));
      static ::ROOT::TGenericClassInfo 
         instance("LitLifetime", "include/Literature.h", 10,
                  typeid(::LitLifetime), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LitLifetime_Dictionary, isa_proxy, 4,
                  sizeof(::LitLifetime) );
      instance.SetNew(&new_LitLifetime);
      instance.SetNewArray(&newArray_LitLifetime);
      instance.SetDelete(&delete_LitLifetime);
      instance.SetDeleteArray(&deleteArray_LitLifetime);
      instance.SetDestructor(&destruct_LitLifetime);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LitLifetime*)
   {
      return GenerateInitInstanceLocal((::LitLifetime*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LitLifetime*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LitLifetime_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LitLifetime*)0x0)->GetClass();
      LitLifetime_TClassManip(theClass);
   return theClass;
   }

   static void LitLifetime_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/Literature.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LitBranchingRatio_Dictionary();
   static void LitBranchingRatio_TClassManip(TClass*);
   static void *new_LitBranchingRatio(void *p = 0);
   static void *newArray_LitBranchingRatio(Long_t size, void *p);
   static void delete_LitBranchingRatio(void *p);
   static void deleteArray_LitBranchingRatio(void *p);
   static void destruct_LitBranchingRatio(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LitBranchingRatio*)
   {
      ::LitBranchingRatio *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LitBranchingRatio));
      static ::ROOT::TGenericClassInfo 
         instance("LitBranchingRatio", "include/Literature.h", 51,
                  typeid(::LitBranchingRatio), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LitBranchingRatio_Dictionary, isa_proxy, 4,
                  sizeof(::LitBranchingRatio) );
      instance.SetNew(&new_LitBranchingRatio);
      instance.SetNewArray(&newArray_LitBranchingRatio);
      instance.SetDelete(&delete_LitBranchingRatio);
      instance.SetDeleteArray(&deleteArray_LitBranchingRatio);
      instance.SetDestructor(&destruct_LitBranchingRatio);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LitBranchingRatio*)
   {
      return GenerateInitInstanceLocal((::LitBranchingRatio*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LitBranchingRatio*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LitBranchingRatio_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LitBranchingRatio*)0x0)->GetClass();
      LitBranchingRatio_TClassManip(theClass);
   return theClass;
   }

   static void LitBranchingRatio_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/Literature.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LitMixingRatio_Dictionary();
   static void LitMixingRatio_TClassManip(TClass*);
   static void *new_LitMixingRatio(void *p = 0);
   static void *newArray_LitMixingRatio(Long_t size, void *p);
   static void delete_LitMixingRatio(void *p);
   static void deleteArray_LitMixingRatio(void *p);
   static void destruct_LitMixingRatio(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LitMixingRatio*)
   {
      ::LitMixingRatio *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LitMixingRatio));
      static ::ROOT::TGenericClassInfo 
         instance("LitMixingRatio", "include/Literature.h", 100,
                  typeid(::LitMixingRatio), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LitMixingRatio_Dictionary, isa_proxy, 4,
                  sizeof(::LitMixingRatio) );
      instance.SetNew(&new_LitMixingRatio);
      instance.SetNewArray(&newArray_LitMixingRatio);
      instance.SetDelete(&delete_LitMixingRatio);
      instance.SetDeleteArray(&deleteArray_LitMixingRatio);
      instance.SetDestructor(&destruct_LitMixingRatio);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LitMixingRatio*)
   {
      return GenerateInitInstanceLocal((::LitMixingRatio*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LitMixingRatio*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LitMixingRatio_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LitMixingRatio*)0x0)->GetClass();
      LitMixingRatio_TClassManip(theClass);
   return theClass;
   }

   static void LitMixingRatio_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/Literature.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LitMatrixElement_Dictionary();
   static void LitMatrixElement_TClassManip(TClass*);
   static void *new_LitMatrixElement(void *p = 0);
   static void *newArray_LitMatrixElement(Long_t size, void *p);
   static void delete_LitMatrixElement(void *p);
   static void deleteArray_LitMatrixElement(void *p);
   static void destruct_LitMatrixElement(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LitMatrixElement*)
   {
      ::LitMatrixElement *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LitMatrixElement));
      static ::ROOT::TGenericClassInfo 
         instance("LitMatrixElement", "include/Literature.h", 146,
                  typeid(::LitMatrixElement), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LitMatrixElement_Dictionary, isa_proxy, 4,
                  sizeof(::LitMatrixElement) );
      instance.SetNew(&new_LitMatrixElement);
      instance.SetNewArray(&newArray_LitMatrixElement);
      instance.SetDelete(&delete_LitMatrixElement);
      instance.SetDeleteArray(&deleteArray_LitMatrixElement);
      instance.SetDestructor(&destruct_LitMatrixElement);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LitMatrixElement*)
   {
      return GenerateInitInstanceLocal((::LitMatrixElement*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LitMatrixElement*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LitMatrixElement_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LitMatrixElement*)0x0)->GetClass();
      LitMatrixElement_TClassManip(theClass);
   return theClass;
   }

   static void LitMatrixElement_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/Literature.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *MiscFunctions_Dictionary();
   static void MiscFunctions_TClassManip(TClass*);
   static void *new_MiscFunctions(void *p = 0);
   static void *newArray_MiscFunctions(Long_t size, void *p);
   static void delete_MiscFunctions(void *p);
   static void deleteArray_MiscFunctions(void *p);
   static void destruct_MiscFunctions(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MiscFunctions*)
   {
      ::MiscFunctions *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MiscFunctions));
      static ::ROOT::TGenericClassInfo 
         instance("MiscFunctions", "include/MiscFunctions.h", 17,
                  typeid(::MiscFunctions), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MiscFunctions_Dictionary, isa_proxy, 4,
                  sizeof(::MiscFunctions) );
      instance.SetNew(&new_MiscFunctions);
      instance.SetNewArray(&newArray_MiscFunctions);
      instance.SetDelete(&delete_MiscFunctions);
      instance.SetDeleteArray(&deleteArray_MiscFunctions);
      instance.SetDestructor(&destruct_MiscFunctions);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MiscFunctions*)
   {
      return GenerateInitInstanceLocal((::MiscFunctions*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MiscFunctions*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MiscFunctions_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::MiscFunctions*)0x0)->GetClass();
      MiscFunctions_TClassManip(theClass);
   return theClass;
   }

   static void MiscFunctions_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/MiscFunctions.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *TransitionRates_Dictionary();
   static void TransitionRates_TClassManip(TClass*);
   static void *new_TransitionRates(void *p = 0);
   static void *newArray_TransitionRates(Long_t size, void *p);
   static void delete_TransitionRates(void *p);
   static void deleteArray_TransitionRates(void *p);
   static void destruct_TransitionRates(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TransitionRates*)
   {
      ::TransitionRates *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TransitionRates));
      static ::ROOT::TGenericClassInfo 
         instance("TransitionRates", "include/TransitionRates.h", 24,
                  typeid(::TransitionRates), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TransitionRates_Dictionary, isa_proxy, 4,
                  sizeof(::TransitionRates) );
      instance.SetNew(&new_TransitionRates);
      instance.SetNewArray(&newArray_TransitionRates);
      instance.SetDelete(&delete_TransitionRates);
      instance.SetDeleteArray(&deleteArray_TransitionRates);
      instance.SetDestructor(&destruct_TransitionRates);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TransitionRates*)
   {
      return GenerateInitInstanceLocal((::TransitionRates*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TransitionRates*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TransitionRates_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TransitionRates*)0x0)->GetClass();
      TransitionRates_TClassManip(theClass);
   return theClass;
   }

   static void TransitionRates_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/TransitionRates.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ScalingFitFCN_Dictionary();
   static void ScalingFitFCN_TClassManip(TClass*);
   static void *new_ScalingFitFCN(void *p = 0);
   static void *newArray_ScalingFitFCN(Long_t size, void *p);
   static void delete_ScalingFitFCN(void *p);
   static void deleteArray_ScalingFitFCN(void *p);
   static void destruct_ScalingFitFCN(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ScalingFitFCN*)
   {
      ::ScalingFitFCN *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ScalingFitFCN));
      static ::ROOT::TGenericClassInfo 
         instance("ScalingFitFCN", "include/ScalingFitFCN.h", 8,
                  typeid(::ScalingFitFCN), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ScalingFitFCN_Dictionary, isa_proxy, 4,
                  sizeof(::ScalingFitFCN) );
      instance.SetNew(&new_ScalingFitFCN);
      instance.SetNewArray(&newArray_ScalingFitFCN);
      instance.SetDelete(&delete_ScalingFitFCN);
      instance.SetDeleteArray(&deleteArray_ScalingFitFCN);
      instance.SetDestructor(&destruct_ScalingFitFCN);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ScalingFitFCN*)
   {
      return GenerateInitInstanceLocal((::ScalingFitFCN*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ScalingFitFCN*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ScalingFitFCN_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ScalingFitFCN*)0x0)->GetClass();
      ScalingFitFCN_TClassManip(theClass);
   return theClass;
   }

   static void ScalingFitFCN_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/ScalingFitFCN.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *GOSIASimMinFCN_Dictionary();
   static void GOSIASimMinFCN_TClassManip(TClass*);
   static void delete_GOSIASimMinFCN(void *p);
   static void deleteArray_GOSIASimMinFCN(void *p);
   static void destruct_GOSIASimMinFCN(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GOSIASimMinFCN*)
   {
      ::GOSIASimMinFCN *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::GOSIASimMinFCN));
      static ::ROOT::TGenericClassInfo 
         instance("GOSIASimMinFCN", "include/GOSIASimMinFCN.h", 53,
                  typeid(::GOSIASimMinFCN), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &GOSIASimMinFCN_Dictionary, isa_proxy, 4,
                  sizeof(::GOSIASimMinFCN) );
      instance.SetDelete(&delete_GOSIASimMinFCN);
      instance.SetDeleteArray(&deleteArray_GOSIASimMinFCN);
      instance.SetDestructor(&destruct_GOSIASimMinFCN);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GOSIASimMinFCN*)
   {
      return GenerateInitInstanceLocal((::GOSIASimMinFCN*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GOSIASimMinFCN*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *GOSIASimMinFCN_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::GOSIASimMinFCN*)0x0)->GetClass();
      GOSIASimMinFCN_TClassManip(theClass);
   return theClass;
   }

   static void GOSIASimMinFCN_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/GOSIASimMinFCN.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *GOSIASimFitter_Dictionary();
   static void GOSIASimFitter_TClassManip(TClass*);
   static void *new_GOSIASimFitter(void *p = 0);
   static void *newArray_GOSIASimFitter(Long_t size, void *p);
   static void delete_GOSIASimFitter(void *p);
   static void deleteArray_GOSIASimFitter(void *p);
   static void destruct_GOSIASimFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GOSIASimFitter*)
   {
      ::GOSIASimFitter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::GOSIASimFitter));
      static ::ROOT::TGenericClassInfo 
         instance("GOSIASimFitter", "include/GOSIASimFitter.h", 65,
                  typeid(::GOSIASimFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &GOSIASimFitter_Dictionary, isa_proxy, 4,
                  sizeof(::GOSIASimFitter) );
      instance.SetNew(&new_GOSIASimFitter);
      instance.SetNewArray(&newArray_GOSIASimFitter);
      instance.SetDelete(&delete_GOSIASimFitter);
      instance.SetDeleteArray(&deleteArray_GOSIASimFitter);
      instance.SetDestructor(&destruct_GOSIASimFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GOSIASimFitter*)
   {
      return GenerateInitInstanceLocal((::GOSIASimFitter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GOSIASimFitter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *GOSIASimFitter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::GOSIASimFitter*)0x0)->GetClass();
      GOSIASimFitter_TClassManip(theClass);
   return theClass;
   }

   static void GOSIASimFitter_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/GOSIASimFitter.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *NucleusReader_Dictionary();
   static void NucleusReader_TClassManip(TClass*);
   static void *new_NucleusReader(void *p = 0);
   static void *newArray_NucleusReader(Long_t size, void *p);
   static void delete_NucleusReader(void *p);
   static void deleteArray_NucleusReader(void *p);
   static void destruct_NucleusReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NucleusReader*)
   {
      ::NucleusReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::NucleusReader));
      static ::ROOT::TGenericClassInfo 
         instance("NucleusReader", "include/NucleusReader.h", 23,
                  typeid(::NucleusReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &NucleusReader_Dictionary, isa_proxy, 4,
                  sizeof(::NucleusReader) );
      instance.SetNew(&new_NucleusReader);
      instance.SetNewArray(&newArray_NucleusReader);
      instance.SetDelete(&delete_NucleusReader);
      instance.SetDeleteArray(&deleteArray_NucleusReader);
      instance.SetDestructor(&destruct_NucleusReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NucleusReader*)
   {
      return GenerateInitInstanceLocal((::NucleusReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NucleusReader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *NucleusReader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::NucleusReader*)0x0)->GetClass();
      NucleusReader_TClassManip(theClass);
   return theClass;
   }

   static void NucleusReader_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/NucleusReader.h");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_Nucleus(void *p) {
      return  p ? new(p) ::Nucleus : new ::Nucleus;
   }
   static void *newArray_Nucleus(Long_t nElements, void *p) {
      return p ? new(p) ::Nucleus[nElements] : new ::Nucleus[nElements];
   }
   // Wrapper around operator delete
   static void delete_Nucleus(void *p) {
      delete ((::Nucleus*)p);
   }
   static void deleteArray_Nucleus(void *p) {
      delete [] ((::Nucleus*)p);
   }
   static void destruct_Nucleus(void *p) {
      typedef ::Nucleus current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Nucleus

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExptData(void *p) {
      return  p ? new(p) ::ExptData : new ::ExptData;
   }
   static void *newArray_ExptData(Long_t nElements, void *p) {
      return p ? new(p) ::ExptData[nElements] : new ::ExptData[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExptData(void *p) {
      delete ((::ExptData*)p);
   }
   static void deleteArray_ExptData(void *p) {
      delete [] ((::ExptData*)p);
   }
   static void destruct_ExptData(void *p) {
      typedef ::ExptData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExptData

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExperimentData(void *p) {
      return  p ? new(p) ::ExperimentData : new ::ExperimentData;
   }
   static void *newArray_ExperimentData(Long_t nElements, void *p) {
      return p ? new(p) ::ExperimentData[nElements] : new ::ExperimentData[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExperimentData(void *p) {
      delete ((::ExperimentData*)p);
   }
   static void deleteArray_ExperimentData(void *p) {
      delete [] ((::ExperimentData*)p);
   }
   static void destruct_ExperimentData(void *p) {
      typedef ::ExperimentData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ExperimentData

namespace ROOT {
   // Wrapper around operator delete
   static void delete_DataReader(void *p) {
      delete ((::DataReader*)p);
   }
   static void deleteArray_DataReader(void *p) {
      delete [] ((::DataReader*)p);
   }
   static void destruct_DataReader(void *p) {
      typedef ::DataReader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DataReader

namespace ROOT {
   // Wrapper around operator delete
   static void delete_GOSIAReader(void *p) {
      delete ((::GOSIAReader*)p);
   }
   static void deleteArray_GOSIAReader(void *p) {
      delete [] ((::GOSIAReader*)p);
   }
   static void destruct_GOSIAReader(void *p) {
      typedef ::GOSIAReader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GOSIAReader

namespace ROOT {
   // Wrappers around operator new
   static void *new_ScalingParameter(void *p) {
      return  p ? new(p) ::ScalingParameter : new ::ScalingParameter;
   }
   static void *newArray_ScalingParameter(Long_t nElements, void *p) {
      return p ? new(p) ::ScalingParameter[nElements] : new ::ScalingParameter[nElements];
   }
   // Wrapper around operator delete
   static void delete_ScalingParameter(void *p) {
      delete ((::ScalingParameter*)p);
   }
   static void deleteArray_ScalingParameter(void *p) {
      delete [] ((::ScalingParameter*)p);
   }
   static void destruct_ScalingParameter(void *p) {
      typedef ::ScalingParameter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ScalingParameter

namespace ROOT {
   // Wrappers around operator new
   static void *new_LitLifetime(void *p) {
      return  p ? new(p) ::LitLifetime : new ::LitLifetime;
   }
   static void *newArray_LitLifetime(Long_t nElements, void *p) {
      return p ? new(p) ::LitLifetime[nElements] : new ::LitLifetime[nElements];
   }
   // Wrapper around operator delete
   static void delete_LitLifetime(void *p) {
      delete ((::LitLifetime*)p);
   }
   static void deleteArray_LitLifetime(void *p) {
      delete [] ((::LitLifetime*)p);
   }
   static void destruct_LitLifetime(void *p) {
      typedef ::LitLifetime current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LitLifetime

namespace ROOT {
   // Wrappers around operator new
   static void *new_LitBranchingRatio(void *p) {
      return  p ? new(p) ::LitBranchingRatio : new ::LitBranchingRatio;
   }
   static void *newArray_LitBranchingRatio(Long_t nElements, void *p) {
      return p ? new(p) ::LitBranchingRatio[nElements] : new ::LitBranchingRatio[nElements];
   }
   // Wrapper around operator delete
   static void delete_LitBranchingRatio(void *p) {
      delete ((::LitBranchingRatio*)p);
   }
   static void deleteArray_LitBranchingRatio(void *p) {
      delete [] ((::LitBranchingRatio*)p);
   }
   static void destruct_LitBranchingRatio(void *p) {
      typedef ::LitBranchingRatio current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LitBranchingRatio

namespace ROOT {
   // Wrappers around operator new
   static void *new_LitMixingRatio(void *p) {
      return  p ? new(p) ::LitMixingRatio : new ::LitMixingRatio;
   }
   static void *newArray_LitMixingRatio(Long_t nElements, void *p) {
      return p ? new(p) ::LitMixingRatio[nElements] : new ::LitMixingRatio[nElements];
   }
   // Wrapper around operator delete
   static void delete_LitMixingRatio(void *p) {
      delete ((::LitMixingRatio*)p);
   }
   static void deleteArray_LitMixingRatio(void *p) {
      delete [] ((::LitMixingRatio*)p);
   }
   static void destruct_LitMixingRatio(void *p) {
      typedef ::LitMixingRatio current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LitMixingRatio

namespace ROOT {
   // Wrappers around operator new
   static void *new_LitMatrixElement(void *p) {
      return  p ? new(p) ::LitMatrixElement : new ::LitMatrixElement;
   }
   static void *newArray_LitMatrixElement(Long_t nElements, void *p) {
      return p ? new(p) ::LitMatrixElement[nElements] : new ::LitMatrixElement[nElements];
   }
   // Wrapper around operator delete
   static void delete_LitMatrixElement(void *p) {
      delete ((::LitMatrixElement*)p);
   }
   static void deleteArray_LitMatrixElement(void *p) {
      delete [] ((::LitMatrixElement*)p);
   }
   static void destruct_LitMatrixElement(void *p) {
      typedef ::LitMatrixElement current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LitMatrixElement

namespace ROOT {
   // Wrappers around operator new
   static void *new_MiscFunctions(void *p) {
      return  p ? new(p) ::MiscFunctions : new ::MiscFunctions;
   }
   static void *newArray_MiscFunctions(Long_t nElements, void *p) {
      return p ? new(p) ::MiscFunctions[nElements] : new ::MiscFunctions[nElements];
   }
   // Wrapper around operator delete
   static void delete_MiscFunctions(void *p) {
      delete ((::MiscFunctions*)p);
   }
   static void deleteArray_MiscFunctions(void *p) {
      delete [] ((::MiscFunctions*)p);
   }
   static void destruct_MiscFunctions(void *p) {
      typedef ::MiscFunctions current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MiscFunctions

namespace ROOT {
   // Wrappers around operator new
   static void *new_TransitionRates(void *p) {
      return  p ? new(p) ::TransitionRates : new ::TransitionRates;
   }
   static void *newArray_TransitionRates(Long_t nElements, void *p) {
      return p ? new(p) ::TransitionRates[nElements] : new ::TransitionRates[nElements];
   }
   // Wrapper around operator delete
   static void delete_TransitionRates(void *p) {
      delete ((::TransitionRates*)p);
   }
   static void deleteArray_TransitionRates(void *p) {
      delete [] ((::TransitionRates*)p);
   }
   static void destruct_TransitionRates(void *p) {
      typedef ::TransitionRates current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TransitionRates

namespace ROOT {
   // Wrappers around operator new
   static void *new_ScalingFitFCN(void *p) {
      return  p ? new(p) ::ScalingFitFCN : new ::ScalingFitFCN;
   }
   static void *newArray_ScalingFitFCN(Long_t nElements, void *p) {
      return p ? new(p) ::ScalingFitFCN[nElements] : new ::ScalingFitFCN[nElements];
   }
   // Wrapper around operator delete
   static void delete_ScalingFitFCN(void *p) {
      delete ((::ScalingFitFCN*)p);
   }
   static void deleteArray_ScalingFitFCN(void *p) {
      delete [] ((::ScalingFitFCN*)p);
   }
   static void destruct_ScalingFitFCN(void *p) {
      typedef ::ScalingFitFCN current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ScalingFitFCN

namespace ROOT {
   // Wrapper around operator delete
   static void delete_GOSIASimMinFCN(void *p) {
      delete ((::GOSIASimMinFCN*)p);
   }
   static void deleteArray_GOSIASimMinFCN(void *p) {
      delete [] ((::GOSIASimMinFCN*)p);
   }
   static void destruct_GOSIASimMinFCN(void *p) {
      typedef ::GOSIASimMinFCN current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GOSIASimMinFCN

namespace ROOT {
   // Wrappers around operator new
   static void *new_GOSIASimFitter(void *p) {
      return  p ? new(p) ::GOSIASimFitter : new ::GOSIASimFitter;
   }
   static void *newArray_GOSIASimFitter(Long_t nElements, void *p) {
      return p ? new(p) ::GOSIASimFitter[nElements] : new ::GOSIASimFitter[nElements];
   }
   // Wrapper around operator delete
   static void delete_GOSIASimFitter(void *p) {
      delete ((::GOSIASimFitter*)p);
   }
   static void deleteArray_GOSIASimFitter(void *p) {
      delete [] ((::GOSIASimFitter*)p);
   }
   static void destruct_GOSIASimFitter(void *p) {
      typedef ::GOSIASimFitter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GOSIASimFitter

namespace ROOT {
   // Wrappers around operator new
   static void *new_NucleusReader(void *p) {
      return  p ? new(p) ::NucleusReader : new ::NucleusReader;
   }
   static void *newArray_NucleusReader(Long_t nElements, void *p) {
      return p ? new(p) ::NucleusReader[nElements] : new ::NucleusReader[nElements];
   }
   // Wrapper around operator delete
   static void delete_NucleusReader(void *p) {
      delete ((::NucleusReader*)p);
   }
   static void deleteArray_NucleusReader(void *p) {
      delete [] ((::NucleusReader*)p);
   }
   static void destruct_NucleusReader(void *p) {
      typedef ::NucleusReader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NucleusReader

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 469,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::__1::vector<int, std::__1::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 469,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::__1::vector<float, std::__1::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 469,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETMatrixTlEdoublegRsPgR_Dictionary();
   static void vectorlETMatrixTlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlETMatrixTlEdoublegRsPgR(void *p = 0);
   static void *newArray_vectorlETMatrixTlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlETMatrixTlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlETMatrixTlEdoublegRsPgR(void *p);
   static void destruct_vectorlETMatrixTlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TMatrixT<double> >*)
   {
      vector<TMatrixT<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TMatrixT<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TMatrixT<double> >", -2, "vector", 469,
                  typeid(vector<TMatrixT<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETMatrixTlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TMatrixT<double> >) );
      instance.SetNew(&new_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlETMatrixTlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlETMatrixTlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TMatrixT<double> > >()));

      ::ROOT::AddClassAlternate("vector<TMatrixT<double> >","std::__1::vector<TMatrixT<double>, std::__1::allocator<TMatrixT<double> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TMatrixT<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETMatrixTlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TMatrixT<double> >*)0x0)->GetClass();
      vectorlETMatrixTlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETMatrixTlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TMatrixT<double> > : new vector<TMatrixT<double> >;
   }
   static void *newArray_vectorlETMatrixTlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TMatrixT<double> >[nElements] : new vector<TMatrixT<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      delete ((vector<TMatrixT<double> >*)p);
   }
   static void deleteArray_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      delete [] ((vector<TMatrixT<double> >*)p);
   }
   static void destruct_vectorlETMatrixTlEdoublegRsPgR(void *p) {
      typedef vector<TMatrixT<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TMatrixT<double> >

namespace ROOT {
   static TClass *vectorlEScalingParametergR_Dictionary();
   static void vectorlEScalingParametergR_TClassManip(TClass*);
   static void *new_vectorlEScalingParametergR(void *p = 0);
   static void *newArray_vectorlEScalingParametergR(Long_t size, void *p);
   static void delete_vectorlEScalingParametergR(void *p);
   static void deleteArray_vectorlEScalingParametergR(void *p);
   static void destruct_vectorlEScalingParametergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ScalingParameter>*)
   {
      vector<ScalingParameter> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ScalingParameter>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ScalingParameter>", -2, "vector", 469,
                  typeid(vector<ScalingParameter>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEScalingParametergR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ScalingParameter>) );
      instance.SetNew(&new_vectorlEScalingParametergR);
      instance.SetNewArray(&newArray_vectorlEScalingParametergR);
      instance.SetDelete(&delete_vectorlEScalingParametergR);
      instance.SetDeleteArray(&deleteArray_vectorlEScalingParametergR);
      instance.SetDestructor(&destruct_vectorlEScalingParametergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ScalingParameter> >()));

      ::ROOT::AddClassAlternate("vector<ScalingParameter>","std::__1::vector<ScalingParameter, std::__1::allocator<ScalingParameter> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ScalingParameter>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEScalingParametergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ScalingParameter>*)0x0)->GetClass();
      vectorlEScalingParametergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEScalingParametergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEScalingParametergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ScalingParameter> : new vector<ScalingParameter>;
   }
   static void *newArray_vectorlEScalingParametergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ScalingParameter>[nElements] : new vector<ScalingParameter>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEScalingParametergR(void *p) {
      delete ((vector<ScalingParameter>*)p);
   }
   static void deleteArray_vectorlEScalingParametergR(void *p) {
      delete [] ((vector<ScalingParameter>*)p);
   }
   static void destruct_vectorlEScalingParametergR(void *p) {
      typedef vector<ScalingParameter> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ScalingParameter>

namespace ROOT {
   static TClass *vectorlEMatrixElementgR_Dictionary();
   static void vectorlEMatrixElementgR_TClassManip(TClass*);
   static void *new_vectorlEMatrixElementgR(void *p = 0);
   static void *newArray_vectorlEMatrixElementgR(Long_t size, void *p);
   static void delete_vectorlEMatrixElementgR(void *p);
   static void deleteArray_vectorlEMatrixElementgR(void *p);
   static void destruct_vectorlEMatrixElementgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MatrixElement>*)
   {
      vector<MatrixElement> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MatrixElement>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MatrixElement>", -2, "vector", 469,
                  typeid(vector<MatrixElement>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMatrixElementgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<MatrixElement>) );
      instance.SetNew(&new_vectorlEMatrixElementgR);
      instance.SetNewArray(&newArray_vectorlEMatrixElementgR);
      instance.SetDelete(&delete_vectorlEMatrixElementgR);
      instance.SetDeleteArray(&deleteArray_vectorlEMatrixElementgR);
      instance.SetDestructor(&destruct_vectorlEMatrixElementgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MatrixElement> >()));

      ::ROOT::AddClassAlternate("vector<MatrixElement>","std::__1::vector<MatrixElement, std::__1::allocator<MatrixElement> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MatrixElement>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMatrixElementgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MatrixElement>*)0x0)->GetClass();
      vectorlEMatrixElementgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMatrixElementgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMatrixElementgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatrixElement> : new vector<MatrixElement>;
   }
   static void *newArray_vectorlEMatrixElementgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MatrixElement>[nElements] : new vector<MatrixElement>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMatrixElementgR(void *p) {
      delete ((vector<MatrixElement>*)p);
   }
   static void deleteArray_vectorlEMatrixElementgR(void *p) {
      delete [] ((vector<MatrixElement>*)p);
   }
   static void destruct_vectorlEMatrixElementgR(void *p) {
      typedef vector<MatrixElement> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MatrixElement>

namespace ROOT {
   static TClass *vectorlELitMixingRatiogR_Dictionary();
   static void vectorlELitMixingRatiogR_TClassManip(TClass*);
   static void *new_vectorlELitMixingRatiogR(void *p = 0);
   static void *newArray_vectorlELitMixingRatiogR(Long_t size, void *p);
   static void delete_vectorlELitMixingRatiogR(void *p);
   static void deleteArray_vectorlELitMixingRatiogR(void *p);
   static void destruct_vectorlELitMixingRatiogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LitMixingRatio>*)
   {
      vector<LitMixingRatio> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LitMixingRatio>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LitMixingRatio>", -2, "vector", 469,
                  typeid(vector<LitMixingRatio>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELitMixingRatiogR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LitMixingRatio>) );
      instance.SetNew(&new_vectorlELitMixingRatiogR);
      instance.SetNewArray(&newArray_vectorlELitMixingRatiogR);
      instance.SetDelete(&delete_vectorlELitMixingRatiogR);
      instance.SetDeleteArray(&deleteArray_vectorlELitMixingRatiogR);
      instance.SetDestructor(&destruct_vectorlELitMixingRatiogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LitMixingRatio> >()));

      ::ROOT::AddClassAlternate("vector<LitMixingRatio>","std::__1::vector<LitMixingRatio, std::__1::allocator<LitMixingRatio> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LitMixingRatio>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELitMixingRatiogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LitMixingRatio>*)0x0)->GetClass();
      vectorlELitMixingRatiogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELitMixingRatiogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELitMixingRatiogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitMixingRatio> : new vector<LitMixingRatio>;
   }
   static void *newArray_vectorlELitMixingRatiogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitMixingRatio>[nElements] : new vector<LitMixingRatio>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELitMixingRatiogR(void *p) {
      delete ((vector<LitMixingRatio>*)p);
   }
   static void deleteArray_vectorlELitMixingRatiogR(void *p) {
      delete [] ((vector<LitMixingRatio>*)p);
   }
   static void destruct_vectorlELitMixingRatiogR(void *p) {
      typedef vector<LitMixingRatio> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LitMixingRatio>

namespace ROOT {
   static TClass *vectorlELitMatrixElementgR_Dictionary();
   static void vectorlELitMatrixElementgR_TClassManip(TClass*);
   static void *new_vectorlELitMatrixElementgR(void *p = 0);
   static void *newArray_vectorlELitMatrixElementgR(Long_t size, void *p);
   static void delete_vectorlELitMatrixElementgR(void *p);
   static void deleteArray_vectorlELitMatrixElementgR(void *p);
   static void destruct_vectorlELitMatrixElementgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LitMatrixElement>*)
   {
      vector<LitMatrixElement> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LitMatrixElement>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LitMatrixElement>", -2, "vector", 469,
                  typeid(vector<LitMatrixElement>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELitMatrixElementgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LitMatrixElement>) );
      instance.SetNew(&new_vectorlELitMatrixElementgR);
      instance.SetNewArray(&newArray_vectorlELitMatrixElementgR);
      instance.SetDelete(&delete_vectorlELitMatrixElementgR);
      instance.SetDeleteArray(&deleteArray_vectorlELitMatrixElementgR);
      instance.SetDestructor(&destruct_vectorlELitMatrixElementgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LitMatrixElement> >()));

      ::ROOT::AddClassAlternate("vector<LitMatrixElement>","std::__1::vector<LitMatrixElement, std::__1::allocator<LitMatrixElement> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LitMatrixElement>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELitMatrixElementgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LitMatrixElement>*)0x0)->GetClass();
      vectorlELitMatrixElementgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELitMatrixElementgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELitMatrixElementgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitMatrixElement> : new vector<LitMatrixElement>;
   }
   static void *newArray_vectorlELitMatrixElementgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitMatrixElement>[nElements] : new vector<LitMatrixElement>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELitMatrixElementgR(void *p) {
      delete ((vector<LitMatrixElement>*)p);
   }
   static void deleteArray_vectorlELitMatrixElementgR(void *p) {
      delete [] ((vector<LitMatrixElement>*)p);
   }
   static void destruct_vectorlELitMatrixElementgR(void *p) {
      typedef vector<LitMatrixElement> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LitMatrixElement>

namespace ROOT {
   static TClass *vectorlELitLifetimegR_Dictionary();
   static void vectorlELitLifetimegR_TClassManip(TClass*);
   static void *new_vectorlELitLifetimegR(void *p = 0);
   static void *newArray_vectorlELitLifetimegR(Long_t size, void *p);
   static void delete_vectorlELitLifetimegR(void *p);
   static void deleteArray_vectorlELitLifetimegR(void *p);
   static void destruct_vectorlELitLifetimegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LitLifetime>*)
   {
      vector<LitLifetime> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LitLifetime>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LitLifetime>", -2, "vector", 469,
                  typeid(vector<LitLifetime>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELitLifetimegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LitLifetime>) );
      instance.SetNew(&new_vectorlELitLifetimegR);
      instance.SetNewArray(&newArray_vectorlELitLifetimegR);
      instance.SetDelete(&delete_vectorlELitLifetimegR);
      instance.SetDeleteArray(&deleteArray_vectorlELitLifetimegR);
      instance.SetDestructor(&destruct_vectorlELitLifetimegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LitLifetime> >()));

      ::ROOT::AddClassAlternate("vector<LitLifetime>","std::__1::vector<LitLifetime, std::__1::allocator<LitLifetime> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LitLifetime>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELitLifetimegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LitLifetime>*)0x0)->GetClass();
      vectorlELitLifetimegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELitLifetimegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELitLifetimegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitLifetime> : new vector<LitLifetime>;
   }
   static void *newArray_vectorlELitLifetimegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitLifetime>[nElements] : new vector<LitLifetime>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELitLifetimegR(void *p) {
      delete ((vector<LitLifetime>*)p);
   }
   static void deleteArray_vectorlELitLifetimegR(void *p) {
      delete [] ((vector<LitLifetime>*)p);
   }
   static void destruct_vectorlELitLifetimegR(void *p) {
      typedef vector<LitLifetime> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LitLifetime>

namespace ROOT {
   static TClass *vectorlELitBranchingRatiogR_Dictionary();
   static void vectorlELitBranchingRatiogR_TClassManip(TClass*);
   static void *new_vectorlELitBranchingRatiogR(void *p = 0);
   static void *newArray_vectorlELitBranchingRatiogR(Long_t size, void *p);
   static void delete_vectorlELitBranchingRatiogR(void *p);
   static void deleteArray_vectorlELitBranchingRatiogR(void *p);
   static void destruct_vectorlELitBranchingRatiogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LitBranchingRatio>*)
   {
      vector<LitBranchingRatio> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LitBranchingRatio>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LitBranchingRatio>", -2, "vector", 469,
                  typeid(vector<LitBranchingRatio>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELitBranchingRatiogR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LitBranchingRatio>) );
      instance.SetNew(&new_vectorlELitBranchingRatiogR);
      instance.SetNewArray(&newArray_vectorlELitBranchingRatiogR);
      instance.SetDelete(&delete_vectorlELitBranchingRatiogR);
      instance.SetDeleteArray(&deleteArray_vectorlELitBranchingRatiogR);
      instance.SetDestructor(&destruct_vectorlELitBranchingRatiogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LitBranchingRatio> >()));

      ::ROOT::AddClassAlternate("vector<LitBranchingRatio>","std::__1::vector<LitBranchingRatio, std::__1::allocator<LitBranchingRatio> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<LitBranchingRatio>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELitBranchingRatiogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LitBranchingRatio>*)0x0)->GetClass();
      vectorlELitBranchingRatiogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELitBranchingRatiogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELitBranchingRatiogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitBranchingRatio> : new vector<LitBranchingRatio>;
   }
   static void *newArray_vectorlELitBranchingRatiogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LitBranchingRatio>[nElements] : new vector<LitBranchingRatio>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELitBranchingRatiogR(void *p) {
      delete ((vector<LitBranchingRatio>*)p);
   }
   static void deleteArray_vectorlELitBranchingRatiogR(void *p) {
      delete [] ((vector<LitBranchingRatio>*)p);
   }
   static void destruct_vectorlELitBranchingRatiogR(void *p) {
      typedef vector<LitBranchingRatio> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LitBranchingRatio>

namespace ROOT {
   static TClass *vectorlEExptDatagR_Dictionary();
   static void vectorlEExptDatagR_TClassManip(TClass*);
   static void *new_vectorlEExptDatagR(void *p = 0);
   static void *newArray_vectorlEExptDatagR(Long_t size, void *p);
   static void delete_vectorlEExptDatagR(void *p);
   static void deleteArray_vectorlEExptDatagR(void *p);
   static void destruct_vectorlEExptDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ExptData>*)
   {
      vector<ExptData> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ExptData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ExptData>", -2, "vector", 469,
                  typeid(vector<ExptData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEExptDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ExptData>) );
      instance.SetNew(&new_vectorlEExptDatagR);
      instance.SetNewArray(&newArray_vectorlEExptDatagR);
      instance.SetDelete(&delete_vectorlEExptDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEExptDatagR);
      instance.SetDestructor(&destruct_vectorlEExptDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ExptData> >()));

      ::ROOT::AddClassAlternate("vector<ExptData>","std::__1::vector<ExptData, std::__1::allocator<ExptData> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ExptData>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEExptDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ExptData>*)0x0)->GetClass();
      vectorlEExptDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEExptDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEExptDatagR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ExptData> : new vector<ExptData>;
   }
   static void *newArray_vectorlEExptDatagR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ExptData>[nElements] : new vector<ExptData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEExptDatagR(void *p) {
      delete ((vector<ExptData>*)p);
   }
   static void deleteArray_vectorlEExptDatagR(void *p) {
      delete [] ((vector<ExptData>*)p);
   }
   static void destruct_vectorlEExptDatagR(void *p) {
      typedef vector<ExptData> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ExptData>

namespace ROOT {
   static TClass *vectorlEExperimentDatagR_Dictionary();
   static void vectorlEExperimentDatagR_TClassManip(TClass*);
   static void *new_vectorlEExperimentDatagR(void *p = 0);
   static void *newArray_vectorlEExperimentDatagR(Long_t size, void *p);
   static void delete_vectorlEExperimentDatagR(void *p);
   static void deleteArray_vectorlEExperimentDatagR(void *p);
   static void destruct_vectorlEExperimentDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ExperimentData>*)
   {
      vector<ExperimentData> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ExperimentData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ExperimentData>", -2, "vector", 469,
                  typeid(vector<ExperimentData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEExperimentDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ExperimentData>) );
      instance.SetNew(&new_vectorlEExperimentDatagR);
      instance.SetNewArray(&newArray_vectorlEExperimentDatagR);
      instance.SetDelete(&delete_vectorlEExperimentDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEExperimentDatagR);
      instance.SetDestructor(&destruct_vectorlEExperimentDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ExperimentData> >()));

      ::ROOT::AddClassAlternate("vector<ExperimentData>","std::__1::vector<ExperimentData, std::__1::allocator<ExperimentData> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<ExperimentData>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEExperimentDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ExperimentData>*)0x0)->GetClass();
      vectorlEExperimentDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEExperimentDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEExperimentDatagR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ExperimentData> : new vector<ExperimentData>;
   }
   static void *newArray_vectorlEExperimentDatagR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ExperimentData>[nElements] : new vector<ExperimentData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEExperimentDatagR(void *p) {
      delete ((vector<ExperimentData>*)p);
   }
   static void deleteArray_vectorlEExperimentDatagR(void *p) {
      delete [] ((vector<ExperimentData>*)p);
   }
   static void destruct_vectorlEExperimentDatagR(void *p) {
      typedef vector<ExperimentData> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ExperimentData>

namespace {
  void TriggerDictionaryInitialization_DictOutput_Impl() {
    static const char* headers[] = {
"include/DataReader.h",
"include/ExperimentalInput.h",
"include/GOSIAReader.h",
"include/GOSIASimFitter.h",
"include/GOSIASimMinFCN.h",
"include/Literature.h",
"include/MatrixElement.h",
"include/MiscFunctions.h",
"include/NormalizationSampler.h",
"include/Nucleus.h",
"include/NucleusReader.h",
"include/ScalingFitFCN.h",
"include/ScalingParameter.h",
"include/TransitionRates.h",
0
    };
    static const char* includePaths[] = {
"/Users/jh0110/Programs/root/root-install/include/",
"/Volumes/Seagate_1/GOSIAFitter/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DictOutput dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/Nucleus.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/DataReader.h")))  Nucleus;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/ExperimentalInput.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/DataReader.h")))  ExptData;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/ExperimentalInput.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/DataReader.h")))  ExperimentData;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/DataReader.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/DataReader.h")))  DataReader;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/GOSIAReader.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIAReader.h")))  GOSIAReader;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/ScalingParameter.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  ScalingParameter;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/Literature.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  LitLifetime;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/Literature.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  LitBranchingRatio;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/Literature.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  LitMixingRatio;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/Literature.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  LitMatrixElement;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/MiscFunctions.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  MiscFunctions;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/TransitionRates.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  TransitionRates;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/ScalingFitFCN.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  ScalingFitFCN;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/GOSIASimMinFCN.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  GOSIASimMinFCN;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/GOSIASimFitter.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/GOSIASimFitter.h")))  GOSIASimFitter;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@include/NucleusReader.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/NucleusReader.h")))  NucleusReader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DictOutput dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/DataReader.h"
#include "include/ExperimentalInput.h"
#include "include/GOSIAReader.h"
#include "include/GOSIASimFitter.h"
#include "include/GOSIASimMinFCN.h"
#include "include/Literature.h"
#include "include/MatrixElement.h"
#include "include/MiscFunctions.h"
#include "include/NormalizationSampler.h"
#include "include/Nucleus.h"
#include "include/NucleusReader.h"
#include "include/ScalingFitFCN.h"
#include "include/ScalingParameter.h"
#include "include/TransitionRates.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"DataReader", payloadCode, "@",
"ExperimentData", payloadCode, "@",
"ExptData", payloadCode, "@",
"GOSIAReader", payloadCode, "@",
"GOSIASimFitter", payloadCode, "@",
"GOSIASimMinFCN", payloadCode, "@",
"LitBranchingRatio", payloadCode, "@",
"LitLifetime", payloadCode, "@",
"LitMatrixElement", payloadCode, "@",
"LitMixingRatio", payloadCode, "@",
"MiscFunctions", payloadCode, "@",
"MiscFunctions::c", payloadCode, "@",
"MiscFunctions::hbar", payloadCode, "@",
"Nucleus", payloadCode, "@",
"NucleusReader", payloadCode, "@",
"ScalingFitFCN", payloadCode, "@",
"ScalingParameter", payloadCode, "@",
"TransitionRates", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DictOutput",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DictOutput_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DictOutput_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DictOutput() {
  TriggerDictionaryInitialization_DictOutput_Impl();
}
