file(REMOVE_RECURSE
  "libparasail.a"
  "libparasail.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/parasail.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
