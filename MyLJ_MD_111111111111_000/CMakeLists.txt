cmake_minimum_required(VERSION 3.10)
list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API 2.0 REQUIRED CONFIG)
if(NOT TARGET kim-api)
  enable_testing()
  project("${KIM_API_PROJECT_NAME}" VERSION "${KIM_API_VERSION}"
    LANGUAGES CXX C Fortran)
endif()

# End preamble


set(MODEL_DRIVER_NAME "MyLJ_MD_111111111111_000")

add_kim_api_model_driver_library(
  NAME                    ${MODEL_DRIVER_NAME}
  CREATE_ROUTINE_NAME     "model_driver_create"
  CREATE_ROUTINE_LANGUAGE "cpp"
  )

target_sources(${MODEL_DRIVER_NAME} PRIVATE
  MyLJ.cpp
  )
