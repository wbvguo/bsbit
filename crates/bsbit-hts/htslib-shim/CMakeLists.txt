cmake_minimum_required(VERSION 3.20)
project(bsbit_htslib_shim LANGUAGES C)

include(CTest)
find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)
pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib=1.24)

add_library(bsbit_htslib_shim bsbit_hts.c)
target_compile_features(bsbit_htslib_shim PUBLIC c_std_11)
target_include_directories(
  bsbit_htslib_shim
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(bsbit_htslib_shim PRIVATE PkgConfig::HTSLIB)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(
    bsbit_htslib_shim PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror
  )
endif()

if(BUILD_TESTING)
  add_executable(bsbit_htslib_shim_smoke tests/shim_smoke.c)
  target_compile_features(bsbit_htslib_shim_smoke PRIVATE c_std_11)
  target_link_libraries(
    bsbit_htslib_shim_smoke PRIVATE bsbit_htslib_shim PkgConfig::HTSLIB ZLIB::ZLIB
  )
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
      bsbit_htslib_shim_smoke PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror
    )
  endif()
  add_test(
    NAME bsbit_htslib_shim_smoke
    COMMAND bsbit_htslib_shim_smoke ${CMAKE_CURRENT_BINARY_DIR}/smoke-artifacts
  )

  add_executable(bsbit_htslib_shim_mutation_smoke tests/mutation_smoke.c)
  target_compile_features(bsbit_htslib_shim_mutation_smoke PRIVATE c_std_11)
  target_link_libraries(
    bsbit_htslib_shim_mutation_smoke
    PRIVATE bsbit_htslib_shim PkgConfig::HTSLIB ZLIB::ZLIB
  )
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
      bsbit_htslib_shim_mutation_smoke
      PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror
    )
  endif()
  add_test(
    NAME bsbit_htslib_shim_mutation_smoke
    COMMAND
      bsbit_htslib_shim_mutation_smoke
      ${CMAKE_CURRENT_BINARY_DIR}/mutation-artifacts
  )

  if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(bsbit_htslib_shim_fault_smoke tests/fault_smoke.c)
    target_compile_features(bsbit_htslib_shim_fault_smoke PRIVATE c_std_11)
    target_link_libraries(
      bsbit_htslib_shim_fault_smoke
      PRIVATE bsbit_htslib_shim PkgConfig::HTSLIB ZLIB::ZLIB
    )
    target_compile_options(
      bsbit_htslib_shim_fault_smoke
      PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror
    )
    target_link_options(
      bsbit_htslib_shim_fault_smoke
      PRIVATE
        -Wl,--wrap=calloc
        -Wl,--wrap=malloc
        -Wl,--wrap=bgzf_open
        -Wl,--wrap=bgzf_compression
        -Wl,--wrap=bgzf_read
        -Wl,--wrap=bgzf_write
        -Wl,--wrap=bgzf_flush
        -Wl,--wrap=bgzf_close
        -Wl,--wrap=sam_hdr_parse
        -Wl,--wrap=bam_init1
        -Wl,--wrap=hts_open
        -Wl,--wrap=sam_hdr_write
        -Wl,--wrap=sam_parse1
        -Wl,--wrap=sam_write1
        -Wl,--wrap=hts_close
    )
    add_test(
      NAME bsbit_htslib_shim_fault_smoke
      COMMAND
        bsbit_htslib_shim_fault_smoke
        ${CMAKE_CURRENT_BINARY_DIR}/fault-artifacts
    )
  endif()
endif()
