cmake_minimum_required(VERSION 3.10)
list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API 2.0 REQUIRED CONFIG)
if(NOT TARGET kim-api)
  enable_testing()
  project("${KIM_API_PROJECT_NAME}" VERSION "${KIM_API_VERSION}"
    LANGUAGES CXX C Fortran)
endif()

# End preamble


add_kim_api_model_library(
  NAME            "LJSi_MO_111111111110_000"
  DRIVER_NAME     "MyLJ_MD_111111111111_000"
  PARAMETER_FILES "si.param"
  )
