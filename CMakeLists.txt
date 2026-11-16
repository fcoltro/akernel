cmake_minimum_required(VERSION 3.20)
project(akernel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# akernel is header-only. Rigour of the interval and expansion arithmetic
# depends on strict IEEE-754 semantics: no contraction, no reordering.
add_library(akernel INTERFACE)
target_include_directories(akernel INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(MSVC)
  target_compile_options(akernel INTERFACE /W4 /WX /EHsc /fp:strict /permissive-)
else()
  # Not verified on this toolchain; flags chosen to match /fp:strict.
  target_compile_options(akernel INTERFACE -Wall -Wextra -Werror
                         -ffp-contract=off -fno-fast-math)
endif()

enable_testing()

function(ak_add_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE akernel)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

ak_add_test(escalation_policy_test tests/num/escalation_policy_test.cpp)
ak_add_test(interval_test          tests/num/interval_test.cpp)
ak_add_test(expansion_test         tests/num/expansion_test.cpp)
ak_add_test(predicates_test        tests/num/predicates_test.cpp)
ak_add_test(bernstein_test         tests/geom/bernstein_test.cpp)
ak_add_test(nurbs_test             tests/geom/nurbs_test.cpp)
ak_add_test(continuity_test        tests/geom/continuity_test.cpp)
ak_add_test(seam_certificate_test  tests/geom/seam_certificate_test.cpp)
ak_add_test(ssi_test               tests/solve/ssi_test.cpp)
ak_add_test(fill_test              tests/var/fill_test.cpp)
ak_add_test(complex_test           tests/topo/complex_test.cpp)
ak_add_test(rate_test              tests/diag/rate_test.cpp)
ak_add_test(csi_test               tests/solve/csi_test.cpp)
ak_add_test(boolean_test           tests/op/boolean_test.cpp)
ak_add_test(tessellate_test        tests/op/tessellate_test.cpp)
ak_add_test(identity_test          tests/op/identity_test.cpp)
ak_add_test(nonmanifold_test       tests/op/nonmanifold_test.cpp)
ak_add_test(quadric_test           tests/op/quadric_test.cpp)
ak_add_test(contact_test           tests/op/contact_test.cpp)
ak_add_test(flat_test              tests/op/flat_test.cpp)
ak_add_test(exact_test             tests/op/exact_test.cpp)

add_executable(ak_bench bench/bench_main.cpp)
target_link_libraries(ak_bench PRIVATE akernel)
target_include_directories(ak_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)

add_executable(ak_demo examples/intersect_demo.cpp)
target_link_libraries(ak_demo PRIVATE akernel)

add_executable(ak_boolean_demo examples/boolean_demo.cpp)
target_link_libraries(ak_boolean_demo PRIVATE akernel)
