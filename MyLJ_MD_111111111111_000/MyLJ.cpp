#include "MyLJ.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

extern "C" {
int model_driver_create(KIM::ModelDriverCreate * const modelDriverCreate,
                        KIM::LengthUnit const requestedLengthUnit,
                        KIM::EnergyUnit const requestedEnergyUnit,
                        KIM::ChargeUnit const requestedChargeUnit,
                        KIM::TemperatureUnit const requestedTemperatureUnit,
                        KIM::TimeUnit const requestedTimeUnit)
{
    int ier; // error code

    LennardJones612 * modelObject;
    modelObject = new LennardJones612(modelDriverCreate,
                                      requestedLengthUnit,
                                      requestedEnergyUnit,
                                      requestedChargeUnit,
                                      requestedTemperatureUnit,
                                      requestedTimeUnit,
                                      & ier); // allocate memory 
                                              // hopefully only new you will even need
    if (ier != 0){
        delete modelObject;
        return ier;
    }

    modelDriverCreate->SetModelBufferPointer(static_cast<void *>(modelObject));
    return 0;

}
} // extern C

LennardJones612::LennardJones612(KIM::ModelDriverCreate * const modelDriverCreate,
                  KIM::LengthUnit const requestedLengthUnit,
                  KIM::EnergyUnit const requestedEnergyUnit,
                  KIM::ChargeUnit const requestedChargeUnit,
                  KIM::TemperatureUnit const requestedTemperatureUnit,
                  KIM::TimeUnit const requestedTimeUnit,
                  int * const ier){
    // load the parameters
    int numberOfParamFiles;
    modelDriverCreate->GetNumberOfParameterFiles(&numberOfParamFiles);
    if (numberOfParamFiles != 1){
        std::cerr << "Wrong number of parameter files, 1 expected, got " << numberOfParamFiles << "\n";
        *ier = 1;
        return;
    }

    *ier = modelDriverCreate->SetModelNumbering(KIM::NUMBERING::zeroBased);
    // Convert units
    KIM::LengthUnit fromLength = KIM::LENGTH_UNIT::A;
    KIM::EnergyUnit fromEnergy = KIM::ENERGY_UNIT::eV;
    KIM::ChargeUnit fromCharge = KIM::CHARGE_UNIT::e;
    KIM::TemperatureUnit fromTemperature = KIM::TEMPERATURE_UNIT::K;
    KIM::TimeUnit fromTime = KIM::TIME_UNIT::ps;

    double convertLength = 1.0, convertEnergy = 1.0;

    *ier = KIM::ModelDriverCreate::ConvertUnit(fromLength,
                                               fromEnergy,
                                               fromCharge,
                                               fromTemperature,
                                               fromTime,
                                               requestedLengthUnit,
                                               requestedEnergyUnit,
                                               requestedChargeUnit,
                                               requestedTemperatureUnit,
                                               requestedTimeUnit,
                                               1.0,
                                               0.0,
                                               0.0,
                                               0.0,
                                               0.0,
                                               &convertLength);

     *ier = KIM::ModelDriverCreate::ConvertUnit(fromLength,
                                               fromEnergy,
                                               fromCharge,
                                               fromTemperature,
                                               fromTime,
                                               requestedLengthUnit,
                                               requestedEnergyUnit,
                                               requestedChargeUnit,
                                               requestedTemperatureUnit,
                                               requestedTimeUnit,
                                               0.0,
                                               1.0,
                                               0.0,
                                               0.0,
                                               0.0,
                                               &convertEnergy);


    *ier = modelDriverCreate->SetUnits(requestedLengthUnit,
                                       requestedEnergyUnit,
                                       KIM::CHARGE_UNIT::unused,
                                       KIM::TEMPERATURE_UNIT::unused,
                                       KIM::TIME_UNIT::unused);


    modelDriverCreate->SetInfluenceDistancePointer(&cutoff);

    const int modelWillNotRequestNeighborsOfNoncontributingParticles_ = static_cast<int>(true);
    modelDriverCreate->SetNeighborListPointers(1, &cutoff, &modelWillNotRequestNeighborsOfNoncontributingParticles_);


    // check if all the functions have correct signature
    KIM::ModelDestroyFunction * destroy = LennardJones612::Destroy; // can be any function with same signature
    KIM::ModelRefreshFunction * refresh = LennardJones612::Refresh;
    KIM::ModelComputeFunction * compute = LennardJones612::Compute;
    KIM::ModelComputeArgumentsCreateFunction * computeArgumentsCreate = LennardJones612::ComputeArgumentsCreate;
    KIM::ModelComputeArgumentsDestroyFunction * computeArgumentsDestroy = LennardJones612::ComputeArgumentsDestroy;


    KIM::SpeciesName KIMSpeciesCode("Si");
    *ier = modelDriverCreate->SetSpeciesCode(KIMSpeciesCode, 0); // only one species so 0


    // register the functions so KIM-API knows what to call
    *ier =   modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Compute,
                                                KIM::LANGUAGE_NAME::cpp,
                                                true,
                                                reinterpret_cast<KIM::Function *>(compute)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate,
                                                KIM::LANGUAGE_NAME::cpp,
                                                true,
                                                reinterpret_cast<KIM::Function *>(computeArgumentsCreate)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy,
                                                KIM::LANGUAGE_NAME::cpp,
                                                true,
                                                reinterpret_cast<KIM::Function *>(computeArgumentsDestroy)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Destroy, 
                                                KIM::LANGUAGE_NAME::cpp,
                                                true,
                                                reinterpret_cast<KIM::Function *>(destroy));


    const std::string * parmFileDir, * paramFileName;
    modelDriverCreate->GetParameterFileDirectoryName(&parmFileDir);
    *ier = modelDriverCreate->GetParameterFileBasename(0, &paramFileName);
    if (ier !=0){
        return;
    }


    std::string filename = *parmFileDir + "/" + *paramFileName;

    std::ifstream paramFile(filename.c_str()); // old C++ compatibility
    std::stringstream buffer;
    // check if file is open
    if (!paramFile.is_open()){
        std::cerr << "Error opening parameter file: " << filename << "\n";
        *ier = 1;
        return;
    }
    buffer << paramFile.rdbuf();
    paramFile.close();

    buffer >> species;
    buffer >> cutoff;
    buffer >> epsilon;
    buffer >> sigma;

    cutoff *= convertLength;
    sigma *= convertLength;
    epsilon *= convertEnergy;



    if (ier != 0){
        std::cerr << "Error setting species code\n";
        return;
    }
    // Species mixing
    *ier = modelDriverCreate->SetParameterPointer(
            1,
            &cutoff,
            "cutoff",
            "Cutoff of the LJ model");
    *ier = modelDriverCreate->SetParameterPointer(
            1,
            &sigma,
            "sigma",
            "Sigma for the LJ model");
    *ier = modelDriverCreate->SetParameterPointer(
            1,
            &epsilon,
            "eps",
            "Epsilon for LJ");

    if (ier != 0){
        std::cerr << "Error setting routine pointers\n";
        return;
    }

    *ier = modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Refresh,
                                                KIM::LANGUAGE_NAME::cpp,
                                                true,
                                                reinterpret_cast<KIM::Function *>(refresh));
}
LennardJones612::~LennardJones612(){};

int LennardJones612::Destroy(KIM::ModelDestroy * const modelDestroy){
    // nothing to destroy
    return 0;
};
int LennardJones612::Refresh(KIM::ModelRefresh * const modelRefresh){
    // TODO: implement refresh if needed
    // Right now it is not
    return 0;
};

int LennardJones612::ComputeArgumentsCreate(
      KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate){


    LennardJones612 * modelObject = NULL;
    modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));
    modelObject->compute_arguments_create(modelComputeArgumentsCreate);

}; 

int LennardJones612::ComputeArgumentsDestroy(
      KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy){
    // nothing to destroy
    return 0;
};

int
  LennardJones612::Compute(KIM::ModelCompute const * const modelCompute,
          KIM::ModelComputeArguments const * const modelComputeArguments){
  int const * numberOfParticles = NULL,
            * particleSpeciesCodes = NULL, 
            * particleContributing= NULL; // nullptr in modern C++
  double const * coordinates = NULL;
  double       * forces = NULL,
               * energy = NULL,
               * particleEnergy = NULL;
  int ier;
  ier = modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &numberOfParticles) + 
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
            &particleSpeciesCodes) +
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
            &particleContributing) + 
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::coordinates,
            (double const **) &coordinates) +
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &energy) +
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, &particleEnergy) +
        modelComputeArguments->GetArgumentPointer(
            KIM::COMPUTE_ARGUMENT_NAME::partialForces,
            (double const **) &forces);
    

    // parameters

    LennardJones612 * modelObject = NULL;
    modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));
    if (modelObject == NULL){
        std::cerr << "Error getting model object\n";
        return 1;
    }
    double sigma = modelObject->sigma; // static parameters issue
    double epsilon = modelObject->epsilon;
    double cutoff = modelObject->cutoff;

    if (ier != 0){
        std::cerr << "Error getting argument pointers\n";
        return ier;
    }

    int numnei = 0;
    int const * n1atom = NULL;

    for (int i = 0; i < *numberOfParticles; i++){
        if (particleContributing[i] == 0){
            continue;
        }

        modelComputeArguments->GetNeighborList(0, //# neighbor list
                                               i, // particle index
                                               &numnei, // number of neighbors
                                               &n1atom); // neighbor list
        for (int j = 0; j < numnei; j++){
            int jindex = n1atom[j];
            double dx = coordinates[3 * i] - coordinates[3 * jindex];
            double dy = coordinates[3 * i + 1] - coordinates[3 * jindex + 1];
            double dz = coordinates[3 * i + 2] - coordinates[3 * jindex + 2];

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > cutoff * cutoff){
                continue;
            }
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double r6r12 = sigma * sigma / (r6 * r12);
            double fpair = 24.0 * epsilon * (2.0 * r6r12 - r6) / (r2 * r2);
            if (forces != NULL){
                forces[3*i] += fpair*dx;
                forces[3*i+1] += fpair*dy;
                forces[3*i+2] += fpair*dz;

                forces[3*jindex] -= fpair*dx;
                forces[3*jindex+1] -= fpair*dy;
                forces[3*jindex+2] -= fpair*dz;
            }

            double pair_energy = 2.0 * epsilon * (r6r12 - r6); // pair energy/2
            if (particleContributing[jindex] == 1){
                if (energy != NULL){
                    *energy += pair_energy;
                }
                if (particleEnergy != NULL){
                    particleEnergy[i] += pair_energy;
                    particleEnergy[jindex] += pair_energy;
                }
            } else {
                if (particleEnergy != NULL){
                    particleEnergy[i] += pair_energy;
                }
            }
        }
    }


  };



int LennardJones612::compute_arguments_create(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate){

    int error = modelComputeArgumentsCreate->SetArgumentSupportStatus(
                            KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, 
                            KIM::SUPPORT_STATUS::optional) ||
                modelComputeArgumentsCreate->SetArgumentSupportStatus(
                            KIM::COMPUTE_ARGUMENT_NAME::partialForces, 
                            KIM::SUPPORT_STATUS::optional) ||
                modelComputeArgumentsCreate->SetArgumentSupportStatus(
                            KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
                            KIM::SUPPORT_STATUS::optional);

    std::cout << error << std::endl;
    if (error != 0){
        std::cerr << "Error setting argument support status\n";
        return error;
    }
    error = modelComputeArgumentsCreate->SetCallbackSupportStatus(
                            KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
                            KIM::SUPPORT_STATUS::notSupported) ||
            modelComputeArgumentsCreate->SetCallbackSupportStatus(
                            KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
                            KIM::SUPPORT_STATUS::notSupported);
    
    
    std::cout << error << "<<<<" << std::endl;

    // Discuss
    //   modelComputeArguments->IsCallbackPresent(
    //   KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm, &compProcess_dEdr);
    //   modelComputeArguments->IsCallbackPresent(
    //   KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term, &compProcess_d2Edr2);
 
    if (error != 0){
        std::cout << error << "<<<" << std::endl;
        std::cerr << "Error setting callback support status\n";
        return error;
    }
};