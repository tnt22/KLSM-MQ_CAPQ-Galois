function(rgraph base numnodes numedges a b c d)
  add_custom_command(OUTPUT ${base}.dimacs
    COMMAND python ${CMAKE_BINARY_DIR}/tools/generators/rmat.py
      ${CMAKE_BINARY_DIR}/tools/bin/GTgraph-rmat -n ${numnodes} -m ${numedges} -a ${a} -b ${b} -c ${c} -d ${d} ${base}.rmat
    DEPENDS ${CMAKE_BINARY_DIR}/tools/bin/GTgraph-rmat)
  add_custom_command(OUTPUT ${base}.gr
    COMMAND graph-convert-standalone -dimacs2gr ${base}.rmat ${base}.gr
    DEPENDS ${base}.dimacs graph-convert-standalone)
endfunction(rgraph)

rgraph(rmat8-2e14 16384 131072 0.45 0.15 0.15 0.25)
rgraph(rmat8-2e18 262144 2097152 0.45 0.15 0.15 0.25)
rgraph(rmat16-2e24 16777216 268435456 0.45 0.15 0.15 0.25)
rgraph(rmat16p-2e24 16777216 268435456 0.5 0.1 0.1 0.3)
rgraph(rmat16p-2e27 134217728 2147483648 0.57 0.19 0.19 0.05)
rgraph(rmat8-2e26 67108864 536870912 0.45 0.15 0.15 0.25)

add_custom_target(more-scalefree-graphs DEPENDS rmat16-2e24.gr rmat16p-2e24.gr rmat16p-2e27.gr)
add_dependencies(more-inputs more-scalefree-graphs)
add_dependencies(more-scalefree-graphs more-tools)
