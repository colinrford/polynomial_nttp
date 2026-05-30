# LamFeatureDetect.cmake — shared build-feature detection for lam submodules.
#
# Master copy at <workspace>/cmake/. Vendored per-submodule (same convention as
# LamModule.cmake / GenerateLamConfig.cmake) so standalone builds stay
# independent. Provides:
#   lam_check_int128(<out_bool> [TYPE_ALIAS_VAR <v> TYPE_ALIAS_NAME <n>])
#   lam_setup_tbb(<out_bool> [DEFAULT ON|OFF] [DEFINE <macro>] [QUIET])

if(_LAM_FEATURE_DETECT_INCLUDED)
  return()
endif()
set(_LAM_FEATURE_DETECT_INCLUDED TRUE)

# -----------------------------------------------------------------------------
# lam_check_int128(<out_bool> [TYPE_ALIAS_VAR <var> TYPE_ALIAS_NAME <name>])
#
# Sets <out_bool> (in PARENT_SCOPE) to the string "true"/"false" for whether
# `unsigned __int128` is available, suitable for substitution into a *_config
# template's `inline constexpr bool ... = @...@;`.
#
# Uses a compiler/platform branch rather than check_cxx_source_compiles: the
# latter inherits CMAKE_CXX_MODULE_STD ON and fails in the TryCompile scratch
# dir (learned in ctbignum). __int128 is a well-known Clang/GCC extension on
# 64-bit targets; MSVC never has it. (32-bit Clang/GCC is out of scope here.)
#
# If TYPE_ALIAS_VAR + TYPE_ALIAS_NAME are given, also sets <var> to a C++ alias
# string: "using <name> = unsigned __int128;" or "using <name> = void;".
# -----------------------------------------------------------------------------
function(lam_check_int128 out_bool)
  cmake_parse_arguments(PARSE_ARGV 1 _a "" "TYPE_ALIAS_VAR;TYPE_ALIAS_NAME" "")

  if(MSVC)
    set(_has "false")
    set(_alias_t "void")
  else()
    set(_has "true")
    set(_alias_t "unsigned __int128")
  endif()

  set(${out_bool} "${_has}" PARENT_SCOPE)
  if(_a_TYPE_ALIAS_VAR AND _a_TYPE_ALIAS_NAME)
    set(${_a_TYPE_ALIAS_VAR} "using ${_a_TYPE_ALIAS_NAME} = ${_alias_t};" PARENT_SCOPE)
  endif()
endfunction()

# -----------------------------------------------------------------------------
# lam_setup_tbb(<out_bool> [DEFAULT ON|OFF] [DEFINE <macro>] [QUIET])
#
# Declares option(LAM_USE_TBB) (DEFAULT defaults to OFF), and if enabled runs
# find_package(TBB). Sets <out_bool> (PARENT_SCOPE) to "true"/"false" for the
# config template. If found and DEFINE is given, adds that compile definition to
# the calling directory. If requested-but-not-found, flips LAM_USE_TBB OFF in
# the caller's scope so a subsequent `if(LAM_USE_TBB) target_link_libraries(...
# TBB::tbb)` is skipped. QUIET passes through to find_package.
# -----------------------------------------------------------------------------
function(lam_setup_tbb out_bool)
  cmake_parse_arguments(PARSE_ARGV 1 _a "QUIET" "DEFAULT;DEFINE" "")
  if(NOT DEFINED _a_DEFAULT)
    set(_a_DEFAULT OFF)
  endif()
  option(LAM_USE_TBB "Enable TBB for parallel/bulk operations" ${_a_DEFAULT})

  set(_bool "false")
  if(LAM_USE_TBB)
    if(_a_QUIET)
      find_package(TBB QUIET)
    else()
      find_package(TBB)
    endif()
    if(TBB_FOUND)
      message(STATUS "TBB found (${TBB_VERSION}): enabling TBB support")
      set(_bool "true")
      if(_a_DEFINE)
        add_compile_definitions(${_a_DEFINE})
      endif()
    else()
      message(STATUS "TBB not found: disabling TBB support")
      set(LAM_USE_TBB OFF PARENT_SCOPE)
    endif()
  endif()
  set(${out_bool} "${_bool}" PARENT_SCOPE)
endfunction()
