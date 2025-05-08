#include "MyLJ.hpp"
#include <iostream>  // For standard input/output operations (like std::cout, std::cerr)
#include <sstream>   // For string stream operations (like std::stringstream)
#include <fstream>   // For file input/output operations (like std::ifstream)
#include <string>    // For string manipulation (like std::string)

// The extern "C" block is used to prevent C++ name mangling for the
// enclosed functions. This is crucial when creating a library that needs
// to be callable from C or other languages that expect C-style linkage.
// In the context of KIM, the API often expects these C-style entry points.
extern "C" {

// model_driver_create is a mandatory function for KIM models.
// It's the entry point called by the KIM API when a model is instantiated.
// Its purpose is to create an instance of the model and set it up.
int model_driver_create(KIM::ModelDriverCreate * const modelDriverCreate, // Pointer to KIM object for model creation
                        KIM::LengthUnit const requestedLengthUnit,        // Requested unit for length
                        KIM::EnergyUnit const requestedEnergyUnit,        // Requested unit for energy
                        KIM::ChargeUnit const requestedChargeUnit,        // Requested unit for charge
                        KIM::TemperatureUnit const requestedTemperatureUnit, // Requested unit for temperature
                        KIM::TimeUnit const requestedTimeUnit)            // Requested unit for time
{
    int ier; // Variable to store error codes returned by KIM API calls or internal checks.

    // Declare a pointer to our LennardJones612 class.
    LennardJones612 * modelObject;
    // Allocate memory for a new LennardJones612 object.
    // The constructor of LennardJones612 is called here, passing along the
    // KIM creation object and requested units.
    // The 'new' keyword dynamically allocates memory on the heap.
    // The address of 'ier' is passed to the constructor so it can report errors.
    modelObject = new LennardJones612(modelDriverCreate,
                                      requestedLengthUnit,
                                      requestedEnergyUnit,
                                      requestedChargeUnit,
                                      requestedTemperatureUnit,
                                      requestedTimeUnit,
                                      &ier); // Pass address of ier for error reporting

    // Check if an error occurred during the model object's construction.
    if (ier != 0) {
        // If an error occurred, it's important to clean up allocated memory
        // to prevent memory leaks. 'delete' calls the destructor and frees memory.
        delete modelObject;
        // Return the error code to the KIM API.
        return ier;
    }

    // If construction was successful, store a pointer to the created model object
    // within the KIM framework. This allows KIM to access this model instance later.
    // static_cast<void *> is used to cast the modelObject pointer to a generic void pointer,
    // static_cast is a safer way to convert between types in C++, e.g. 
    // double to int, or in this case, LennardJones612* to void*.
    // Because KIM stores it in a generic way.
    modelDriverCreate->SetModelBufferPointer(static_cast<void *>(modelObject));
    // save the pointer to the model object in the KIM API

    // Return 0 to indicate successful creation.
    return 0;
}
} // end extern "C"

// Constructor for the LennardJones612 class.
// This is called when 'new LennardJones612(...)' is executed in model_driver_create.
LennardJones612::LennardJones612(KIM::ModelDriverCreate * const modelDriverCreate,
                                 KIM::LengthUnit const requestedLengthUnit,
                                 KIM::EnergyUnit const requestedEnergyUnit,
                                 KIM::ChargeUnit const requestedChargeUnit,
                                 KIM::TemperatureUnit const requestedTemperatureUnit,
                                 KIM::TimeUnit const requestedTimeUnit,
                                 int * const ier) { // Pointer to an integer for error reporting
    // Initialize the error code to 0 (no error).
    *ier = 0;

    // --- Load Parameters from File ---
    int numberOfParamFiles;
    const std::string *parmFileDir, *paramFileName; // Pointers to strings for directory and filename

    // Get the directory where parameter files are located.
    modelDriverCreate->GetParameterFileDirectoryName(&parmFileDir);
    // Get the base name of the first parameter file (index 0).
    // only one parameter file is expected for this model.
    *ier = modelDriverCreate->GetParameterFileBasename(0, &paramFileName);
    if (*ier != 0) {
        std::cerr << "Error getting parameter file basename.\n";
        return; // Exit constructor if there's an error
    }

    // Get the total number of parameter files provided.
    modelDriverCreate->GetNumberOfParameterFiles(&numberOfParamFiles);
    if (numberOfParamFiles != 1) {
        std::cerr << "Wrong number of parameter files, 1 expected, got " << numberOfParamFiles << "\n";
        *ier = 1; // Set error code
        return;   // Exit constructor
    }

    // Construct the full path to the parameter file.
    std::string filename = *parmFileDir + "/" + *paramFileName;

    // Open the parameter file for reading.
    // filename.c_str() is used for compatibility with older C++ versions
    // or C APIs that expect a const char*.
    std::ifstream paramFile(filename.c_str());
    std::stringstream buffer; // Use a stringstream to read the file content easily.
    // It provides a convenient way to convert content to different datatypes

    // Check if the parameter file was successfully opened.
    if (!paramFile.is_open()) {
        std::cerr << "Error opening parameter file: " << filename << "\n";
        *ier = 1; // Set error code
        return;   // Exit constructor
    }

    // Set the model's internal numbering scheme (e.g., zero-based indexing for arrays).
    *ier = modelDriverCreate->SetModelNumbering(KIM::NUMBERING::zeroBased);
    if (*ier != 0) {
        std::cerr << "Error setting model numbering.\n";
        return;
    }

    // Define the units in which the parameters in the parameter file are specified.
    // These are the "from" units for conversion. (default units for this model)
    KIM::LengthUnit fromLength = KIM::LENGTH_UNIT::A;          // Angstroms
    KIM::EnergyUnit fromEnergy = KIM::ENERGY_UNIT::eV;          // electronVolts
    KIM::ChargeUnit fromCharge = KIM::CHARGE_UNIT::e;          // elementary charge (not used in LJ)
    KIM::TemperatureUnit fromTemperature = KIM::TEMPERATURE_UNIT::K; // Kelvin (not used in this LJ)
    KIM::TimeUnit fromTime = KIM::TIME_UNIT::ps;                // picoseconds (not used in this LJ)

    // Initialize conversion factors.
    double convertLength = 1.0;
    double convertEnergy = 1.0;

    // Calculate the conversion factor for length.
    // This converts a value of 1.0 from 'fromLength' to 'requestedLengthUnit'.
    *ier = KIM::ModelDriverCreate::ConvertUnit(fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
                                              requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
                                              requestedTemperatureUnit, requestedTimeUnit,
                                              1.0, 0.0, 0.0, 0.0, 0.0, // Value to convert: 1.0 length unit
                                              &convertLength);         // Output conversion factor
    if (*ier != 0) {
        std::cerr << "Error during length unit conversion factor calculation.\n";
        return;
    }

    // Calculate the conversion factor for energy.
    // This converts a value of 1.0 from 'fromEnergy' to 'requestedEnergyUnit'.
    *ier = KIM::ModelDriverCreate::ConvertUnit(fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
                                              requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
                                              requestedTemperatureUnit, requestedTimeUnit,
                                              0.0, 1.0, 0.0, 0.0, 0.0, // Value to convert: 1.0 energy unit
                                              &convertEnergy);         // Output conversion factor
    if (*ier != 0) {
        std::cerr << "Error during energy unit conversion factor calculation.\n";
        return;
    }

    // Inform the KIM API about the units this model will use internally for calculations
    // after conversion. Charge, Temperature, and Time are marked as unused for this LJ model.
    *ier = modelDriverCreate->SetUnits(requestedLengthUnit,
                                      requestedEnergyUnit,
                                      KIM::CHARGE_UNIT::unused,
                                      KIM::TEMPERATURE_UNIT::unused,
                                      KIM::TIME_UNIT::unused);
    if (*ier != 0) {
        std::cerr << "Error setting model units.\n";
        return;
    }

    // --- Read Parameters from File Buffer ---
    // Read the entire content of the parameter file into the stringstream buffer.
    buffer << paramFile.rdbuf();
    paramFile.close(); // Close the file as it's no longer needed.

    // Read the model parameters from the buffer.
    // It's assumed the parameter file contains these values in order:
    // species_name cutoff_distance epsilon_value sigma_value
    buffer >> species; // Reads the species name (e.g., "Ar"), single species
    buffer >> cutoff;  // Reads the cutoff distance for the potential
    buffer >> epsilon; // Reads the energy well depth (epsilon)
    buffer >> sigma;   // Reads the finite distance at which inter-particle potential is zero (sigma)

    // --- KIM API Setup for Species and Parameters ---
    // Convert the C++ string species name to a KIM::SpeciesName object.
    KIM::SpeciesName KIMSpeciesCode(species.c_str());
    // Set the species code for the model. This model supports only one species (index 0).
    *ier = modelDriverCreate->SetSpeciesCode(KIMSpeciesCode, 0);
    if (*ier != 0) {
        std::cerr << "Error setting species code.\n";
        return;
    }

    // Apply the calculated unit conversion factors to the parameters.
    cutoff *= convertLength;
    sigma *= convertLength;
    epsilon *= convertEnergy;

    // Expose model parameters to the KIM API.
    // This allows other tools or the KIM framework to query these parameter values.
    // The '1' indicates that the parameter is a scalar (not an array).
    *ier = modelDriverCreate->SetParameterPointer(
        1,                  // Number of elements (1 for scalar)
        &cutoff,            // Pointer to the cutoff variable
        "cutoff",           // Name of the parameter (as string)
        "Cutoff of the LJ model"); // Description of the parameter
    if (*ier != 0) { std::cerr << "Error setting cutoff parameter pointer.\n"; return; }

    *ier = modelDriverCreate->SetParameterPointer(
        1,
        &sigma,
        "sigma",
        "Sigma for the LJ model");
    if (*ier != 0) { std::cerr << "Error setting sigma parameter pointer.\n"; return; }

    *ier = modelDriverCreate->SetParameterPointer(
        1,
        &epsilon,
        "epsilon", // KIM standard name for epsilon is often "eps" or "epsilon"
        "Epsilon for LJ");
    if (*ier != 0) { std::cerr << "Error setting epsilon parameter pointer.\n"; return; }


    // Set the influence distance for the model. 
    modelDriverCreate->SetInfluenceDistancePointer(&cutoff);

    // Configure neighbor list settings.
    // 'modelWillNotRequestNeighborsOfNoncontributingParticles_' indicates whether the model
    // needs neighbor lists for particles that don't contribute to energy (ghost atoms).
    // For many pair potentials, this can be true.
    const int modelWillNotRequestNeighborsOfNoncontributingParticles_ = static_cast<int>(true);
    // C++ 0 = false, non-zero = true
    // C < 99 did not have bool, so this is a common pattern.
    // static_cast<int>(true) is a way to ensure the value is treated as an integer.
    // The '1' indicates we are setting this for the first (and only) neighbor list type.
    // Typically, this is for the standard neighbor list based on the cutoff.
    modelDriverCreate->SetNeighborListPointers(1, &cutoff, &modelWillNotRequestNeighborsOfNoncontributingParticles_);
    // for multiple neighbor lists, (first param > 1), cutoff will be an array of cutoffs


    // --- Register Model Routines (Function Pointers) ---
    // Somewhat complicated pattern to ensure the correct function pointers are set up.
    // These lines get pointers to the static member functions of the LennardJones612 class.
    // For nw think of this as a copy-paste template
    KIM::ModelDestroyFunction * destroy = LennardJones612::Destroy;
    KIM::ModelRefreshFunction * refresh = LennardJones612::Refresh;
    KIM::ModelComputeFunction * compute = LennardJones612::Compute;
    KIM::ModelComputeArgumentsCreateFunction * computeArgumentsCreate = LennardJones612::ComputeArgumentsCreate;
    KIM::ModelComputeArgumentsDestroyFunction * computeArgumentsDestroy = LennardJones612::ComputeArgumentsDestroy;
    // KIM will call these functions as KIM::ModelDestroyFunction(), so abstract away the class.

    // Register the core functions of the model with the KIM API.
    // This tells KIM which C++ functions to call for specific operations.
    // KIM::LANGUAGE_NAME::cpp indicates these are C++ functions.
    // 'true' indicates that these functions are required by the KIM API.
    // reinterpret_cast<KIM::Function *> casts the specific function pointer to a generic KIM::Function pointer.
    *ier =  modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Compute,
                                                KIM::LANGUAGE_NAME::cpp, true,
                                                reinterpret_cast<KIM::Function *>(compute)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate,
                                                KIM::LANGUAGE_NAME::cpp, true,
                                                reinterpret_cast<KIM::Function *>(computeArgumentsCreate)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy,
                                                KIM::LANGUAGE_NAME::cpp, true,
                                                reinterpret_cast<KIM::Function *>(computeArgumentsDestroy)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Destroy,
                                                KIM::LANGUAGE_NAME::cpp, true,
                                                reinterpret_cast<KIM::Function *>(destroy)) ||
            modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Refresh,
                                                KIM::LANGUAGE_NAME::cpp, true,
                                                reinterpret_cast<KIM::Function *>(refresh));
    if (*ier != 0) {
        std::cerr << "Error setting routine pointers.\n";
        // No return here, as the original code didn't have one, but it might be good practice to add.
    }
}

// Destructor for the LennardJones612 class.
// This is called when 'delete modelObject' is executed.
// Its purpose is to release any resources acquired by the class instance.
LennardJones612::~LennardJones612() {
    // In this simple LJ model, no dynamic memory was allocated directly by the
    // LennardJones612 object itself (beyond what 'new' did for the object itself,
    // which 'delete' handles). If there were, this is where it would be freed.
}

// Static member function: Destroy
// This function is called by the KIM API when the model instance is being destroyed.
// It's responsible for cleaning up anything set up by model_driver_create,
// primarily deleting the model object itself.
int LennardJones612::Destroy(KIM::ModelDestroy * const modelDestroy) {
    // Retrieve the pointer to the model object that was stored during creation.
    LennardJones612 * modelObject = NULL;
    modelDestroy->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

    // If the modelObject pointer is valid, delete the object.
    // This calls the LennardJones612 destructor.
    if (modelObject != NULL) {
        delete modelObject;
    }
    // Return 0 to indicate success.
    return 0;
}

// Static member function: Refresh
// This function is called by the KIM API if the model's parameters might have changed
// externally (e.g., via KIM API calls to modify parameters). The model should
// re-read or update its internal state if necessary.
int LennardJones612::Refresh(KIM::ModelRefresh * const modelRefresh) {
    // TODO: Implement refresh functionality if parameters can be changed
    // dynamically and the model needs to react to those changes.
    // For this simple model, parameters are read once at creation and assumed constant.
    // Thus, Refresh currently does nothing.
    return 0;
}

// Static member function: ComputeArgumentsCreate
// This function is called by the KIM API before a computation.
// Its purpose is to inform the KIM API about which computational arguments
// (like energy, forces, particleEnergy) this model supports and can compute.
// It also specifies support for callbacks.
int LennardJones612::ComputeArgumentsCreate(
    KIM::ModelCompute const * const modelCompute,                       // KIM object for compute context
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) { // KIM object for specifying arg support

    // Retrieve the model object pointer (though not strictly needed in this function
    // as we are just declaring support, not using model parameters).
    LennardJones612 * modelObject = NULL;
    modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

    int error = 0;
    // Declare support status for various compute arguments.
    // KIM::SUPPORT_STATUS::optional means the model can compute these if requested,
    // but the calling simulation might not always ask for them.
    error = modelComputeArgumentsCreate->SetArgumentSupportStatus(
                KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,      // Total potential energy
                KIM::SUPPORT_STATUS::optional) ||
            modelComputeArgumentsCreate->SetArgumentSupportStatus(
                KIM::COMPUTE_ARGUMENT_NAME::partialForces,      // Forces on particles
                KIM::SUPPORT_STATUS::optional) ||
            modelComputeArgumentsCreate->SetArgumentSupportStatus(
                KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, // Per-particle energy
                KIM::SUPPORT_STATUS::optional);

    if (error != 0) {
        std::cerr << "Error setting argument support status in ComputeArgumentsCreate.\n";
        return error;
    }

    // Declare support status for compute callbacks.
    // Callbacks are more advanced KIM features where the model can ask the
    // simulation environment to perform certain actions during computation.
    // This LJ model does not use these specific callbacks.
    error = modelComputeArgumentsCreate->SetCallbackSupportStatus(
                KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,    // Callback for first derivative of energy
                KIM::SUPPORT_STATUS::notSupported) ||
            modelComputeArgumentsCreate->SetCallbackSupportStatus(
                KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,  // Callback for second derivative of energy
                KIM::SUPPORT_STATUS::notSupported);

    if (error != 0) {
        std::cerr << "Error setting callback support status in ComputeArgumentsCreate.\n";
        return error;
    }
    return 0; // Success
}

// Static member function: ComputeArgumentsDestroy
// This function is called by the KIM API after a computation.
// If ComputeArgumentsCreate allocated any resources, this function
// would be responsible for freeing them.
int LennardJones612::ComputeArgumentsDestroy(
    KIM::ModelCompute const * const modelCompute,                                 // KIM object for compute context
    KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy) {     // KIM object for arg destruction
    // In this simple model, ComputeArgumentsCreate did not allocate any dynamic
    // resources, so there's nothing to destroy here.
    // I am bit unsure where this will be used at all, as Compute Arguments are usually handled by the simulator
    return 0; // Success
}

// Static member function: Compute
// This is the core computational routine of the model.
// It's called by the KIM API to calculate properties like energy and forces.
int LennardJones612::Compute(
    KIM::ModelCompute const * const modelCompute,                         // KIM object for compute context
    KIM::ModelComputeArguments const * const modelComputeArguments) {   // KIM object for accessing compute arguments

    // Retrieve the pointer to our LennardJones612 model object.
    // This gives access to the model's parameters (sigma, epsilon, cutoff).
    // In most common cases you will use the PIMPL design pattern and will call the
    // implementation compute function here.
    LennardJones612 * modelObject = NULL; // NULL as nullptr only came in C++11
    modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject)); 
    // Cast the void pointer back to the specific model type.
    // Done by casting the LennardJones612 pointer to a void pointer and passing it
    // to the KIM API for storing the address (all address are same size). But the underlying
    // memory pattern is now of type LennardJones612*. And we can use it as such.

    // If modelObject is null, something went wrong (e.g., it wasn't set correctly).
    if (modelObject == NULL) {
        std::cerr << "Error: Model object is null in Compute.\n";
        return 1; // Indicate an error
    }

    // --- Get Pointers to Simulation Data from KIM ---
    // These pointers will give us access to the atomic configuration and arrays
    // where we need to store calculated results.
    int const * numberOfParticles = NULL;    // Total number of particles
    int const * particleSpeciesCodes = NULL; // Array of species codes for each particle
    int const * particleContributing = NULL; // Array indicating if a particle contributes to the energy
    double const * coordinates = NULL;       // Array of particle coordinates (x1,y1,z1, x2,y2,z2, ...)

    // Pointers to arrays where results should be stored.
    // These will be non-null only if the simulation requested them.
    double * forces = NULL;            // Array for forces (fx1,fy1,fz1, fx2,fy2,fz2, ...)
    double * energy = NULL;            // Pointer to a single double for total energy
    double * particleEnergy = NULL;    // Array for per-particle energies

    int ier; // Error code

    // Request pointers to the necessary arguments from the KIM API.
    // The KIM API manages the memory for these arrays.
    // Note: The original code sums the return codes; KIM typically returns non-zero on error for each call.
    // It's often better to check each call individually or use a more robust error aggregation.
    ier = modelComputeArguments->GetArgumentPointer(
              KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &numberOfParticles);
    if (ier != 0) { std::cerr << "Error getting numberOfParticles pointer.\n"; return ier; }

    ier = modelComputeArguments->GetArgumentPointer(
              KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, &particleSpeciesCodes);
    if (ier != 0) { std::cerr << "Error getting particleSpeciesCodes pointer.\n"; return ier; }

    ier = modelComputeArguments->GetArgumentPointer(
              KIM::COMPUTE_ARGUMENT_NAME::particleContributing, &particleContributing);
    if (ier != 0) { std::cerr << "Error getting particleContributing pointer.\n"; return ier; }

    ier = modelComputeArguments->GetArgumentPointer(
              KIM::COMPUTE_ARGUMENT_NAME::coordinates, (double const **) &coordinates);
    if (ier != 0) { std::cerr << "Error getting coordinates pointer.\n"; return ier; }

    // These are optional arguments. 
    // The pointers (energy, particleEnergy, forces) will remain NULL if not available.
    // The subsequent checks (if (forces != NULL), etc.) handle this.
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &energy);
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, &particleEnergy);
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialForces, (double const **) &forces); // Cast to (double const **) for KIM, then used as double*

    // initialize energy to zero if it is not null
    if (energy != NULL) {
        *energy = 0.0; // Initialize total energy to zero
    }
    // initialize particleEnergy to zero if it is not null
    if (particleEnergy != NULL) {
        for (int i = 0; i < *numberOfParticles; ++i) {
            particleEnergy[i] = 0.0; // Initialize per-particle energy to zero
        }
    }
    // initialize forces to zero if it is not null
    if (forces != NULL) {
        for (int i = 0; i < *numberOfParticles; ++i) {
            forces[3 * i + 0] = 0.0; // Initialize forces to zero
            forces[3 * i + 1] = 0.0;
            forces[3 * i + 2] = 0.0;
        }
    }

    // --- Retrieve Model Parameters ---
    // Get the Lennard-Jones parameters (sigma, epsilon, cutoff) from our modelObject.
    // These were read from the parameter file and converted to consistent units during construction.
    double sigma_val = modelObject->sigma;
    double epsilon_val = modelObject->epsilon;
    double cutoff_val = modelObject->cutoff;
    double cutoff_sq = cutoff_val * cutoff_val;

    // --- Main Computation Loop ---
    // Initialize total energy if requested.
    // The KIM API expects the model to *add* to existing values if "partial" is in the name
    // so total energy and forces are summed globally. But each domain is computed in distributed
    // fashion.

    int numnei = 0;          // Number of neighbors for a particle
    int const * n1atom = NULL; // Pointer to the neighbor list for a particle

    // Loop over each particle 'i' in the system.
    for (int i = 0; i < *numberOfParticles; ++i) {
        // Skip particle 'i' if it's not a "contributing" particle.
        // not needed for energy
        if (particleContributing[i] == 0) {
            continue;
        }

        // Get the neighbor list for particle 'i'.
        // The '0' refers to the first (and in this case, only) type of neighbor list
        // requested/configured by the model (based on the cutoff).
        ier = modelComputeArguments->GetNeighborList(0,       // Neighbor list index (0 for the primary list)
                                                     i,       // Index of the central particle
                                                     &numnei, // Output: number of neighbors found
                                                     &n1atom);// Output: pointer to the list of neighbor indices
        if (ier != 0) {
            std::cerr << "Error getting neighbor list for particle " << i << "\n";
            return ier;
        }

        // Loop over each neighbor 'j' of particle 'i'.
        for (int j = 0; j < numnei; ++j) {
            int jindex = n1atom[j]; // Get the actual index of the j-th neighbor

            if (jindex < i && particleContributing[jindex] == 1) {
                // if j < i (particle already encountered) and jindex is not contributing,
                 continue;
            }

            // Calculate the distance vector components between particle i and particle jindex.
            double dx = coordinates[3 * i + 0] - coordinates[3 * jindex + 0];
            double dy = coordinates[3 * i + 1] - coordinates[3 * jindex + 1];
            double dz = coordinates[3 * i + 2] - coordinates[3 * jindex + 2];

            // Calculate the squared distance.
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > cutoff_sq) {
                continue; // Skip if outside cutoff or self
            }

            double sigma_sq_div_r2 = (sigma_val * sigma_val) / r2; // (sigma/r)^2
            double sigma6_div_r6 = sigma_sq_div_r2 * sigma_sq_div_r2 * sigma_sq_div_r2; // (sigma/r)^6
            double sigma12_div_r12 = sigma6_div_r6 * sigma6_div_r6; // (sigma/r)^12

            double pair_energy = 4.0 * epsilon_val * (sigma12_div_r12 - sigma6_div_r6);
            double f_over_r = (24.0 * epsilon_val / r2) * (2.0 * sigma12_div_r12 - sigma6_div_r6);

            if (particleContributing[jindex] == 0){
                f_over_r *= 0.5; // Only add half of the force to the non-contributing particle
            }

            // If force calculation is requested, add the force contributions.
            if (forces != NULL) {
                double force_x_component = f_over_r * dx;
                double force_y_component = f_over_r * dy;
                double force_z_component = f_over_r * dz;

                // Add force to particle i
                forces[3 * i + 0] += force_x_component;
                forces[3 * i + 1] += force_y_component;
                forces[3 * i + 2] += force_z_component;

                forces[3 * jindex + 0] -= force_x_component;
                forces[3 * jindex + 1] -= force_y_component;
                forces[3 * jindex + 2] -= force_z_component;
            }

            double half_pair_interaction_energy = pair_energy * 0.5;

            if (energy != NULL) {
                if (particleContributing[jindex] == 1) { // Only add to j if it's a real particle
                    *energy += pair_energy;
                } else {
                    *energy += half_pair_interaction_energy;
                }
            }

            // If per-particle energy calculation is requested.
            if (particleEnergy != NULL) {
                // Each particle in an interacting pair gets half of the pair's interaction energy.
                particleEnergy[i] += half_pair_interaction_energy;
                if (particleContributing[jindex] == 1) { // Only add to j if it's a real particle
                    particleEnergy[jindex] += half_pair_interaction_energy;
                }
            }
        } // End of neighbor loop (j)
    }     // End of particle loop (i)

    return 0; // Indicate successful computation
}