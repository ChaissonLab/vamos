# CMake generated Testfile for 
# Source directory: /home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2
# Build directory: /home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_isa "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/build/test_isa")
set_tests_properties(test_isa PROPERTIES  _BACKTRACE_TRIPLES "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;1129;ADD_TEST;/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;0;")
add_test(test_basic "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/build/test_basic")
set_tests_properties(test_basic PROPERTIES  _BACKTRACE_TRIPLES "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;1130;ADD_TEST;/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;0;")
add_test(test_verify "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/build/test_verify" "-f" "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/data/test_small_2.fasta")
set_tests_properties(test_verify PROPERTIES  _BACKTRACE_TRIPLES "/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;1131;ADD_TEST;/home/mchaisso/projects/vamos-pathogenic/vamos/src/parasail-2.6.2/CMakeLists.txt;0;")
