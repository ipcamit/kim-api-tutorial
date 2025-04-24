#ifndef LJ_HPP_ 
#define LJ_HPP_ 

#include "./kim-api/KIM_ModelDriverHeaders.hpp"

extern "C" { // treat this function like a C function
int model_driver_create(KIM::ModelDriverCreate * const modelDriverCreate,
                        KIM::LengthUnit const requestedLengthUnit,
                        KIM::EnergyUnit const requestedEnergyUnit,
                        KIM::ChargeUnit const requestedChargeUnit,
                        KIM::TemperatureUnit const requestedTemperatureUnit,
                        KIM::TimeUnit const requestedTimeUnit);
}

class LennardJones612
{
 public:
  LennardJones612(KIM::ModelDriverCreate * const modelDriverCreate,
                  KIM::LengthUnit const requestedLengthUnit,
                  KIM::EnergyUnit const requestedEnergyUnit,
                  KIM::ChargeUnit const requestedChargeUnit,
                  KIM::TemperatureUnit const requestedTemperatureUnit,
                  KIM::TimeUnit const requestedTimeUnit,
                  int * const ier);
  ~LennardJones612();

   // verbatim
  // no need to make these "extern" since KIM will only access them
  // via function pointers.  "static" is required so that there is not
  // an implicit this pointer added to the prototype by the C++ compiler
  // static = @classmethod, remove "self" for C compatibility
  static int Destroy(KIM::ModelDestroy * const modelDestroy);
  static int Refresh(KIM::ModelRefresh * const modelRefresh);
  static int
  Compute(KIM::ModelCompute const * const modelCompute,
          KIM::ModelComputeArguments const * const modelComputeArguments); // <- main code goes here
  static int ComputeArgumentsCreate(
      KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate); // <- what this model can provide
  static int ComputeArgumentsDestroy(
      KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy);
  private:
  std::string species;
  double cutoff, sigma, epsilon;
  int compute_arguments_create(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate);
};

#endif  // LJ_HPP_
