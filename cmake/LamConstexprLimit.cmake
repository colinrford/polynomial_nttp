# LamConstexprLimit.cmake
# Portable constexpr evaluation step limit:
#   Clang / clang-cl : -fconstexpr-steps=N
#   GCC              : -fconstexpr-ops-limit=N
#   MSVC (cl.exe)    : /constexpr:stepsN   (default is only 1048576, so the
#                                           heavy compile-time tests need this)
#
# LAM_CONSTEXPR_STEPS_CAP (optional cache var): clamp every requested limit to
# this ceiling. Lets a weak machine (e.g. Surface Pro 3) cap heavy compile-time
# work at configure time -- -DLAM_CONSTEXPR_STEPS_CAP=20000000 -- without
# touching call sites. Empty (default) means no cap.

set(LAM_CONSTEXPR_STEPS_CAP "" CACHE STRING
    "Optional ceiling for constexpr step limits (empty = no cap)")

function(lam_target_constexpr_limit target limit)
  if(LAM_CONSTEXPR_STEPS_CAP AND limit GREATER LAM_CONSTEXPR_STEPS_CAP)
    set(limit ${LAM_CONSTEXPR_STEPS_CAP})
  endif()
  target_compile_options(${target} PRIVATE
    "$<$<CXX_COMPILER_ID:Clang,AppleClang>:-fconstexpr-steps=${limit}>"
    "$<$<CXX_COMPILER_ID:GNU>:-fconstexpr-ops-limit=${limit}>"
    "$<$<CXX_COMPILER_ID:MSVC>:/constexpr:steps${limit}>"
  )
endfunction()
