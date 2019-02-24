add_library(capq STATIC
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/capq.c
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/fat_skiplist.c
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/qdlocksrc/locks/qd_lock.c
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/qdlocksrc/locks/tatas_lock.c
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/gc/gc.c
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/gc/ptst.c
)

if(NOT DEFINED ENV{CACHE_LINE_SIZE} )
    message( FATAL_ERROR "need to set CACHE_LINE_SIZE variable (-DCACHE_LINE_SIZE for cmake) for capq (128 should be safe to assume)")
endif()

set_target_properties(capq PROPERTIES COMPILE_FLAGS
    "-w -std=gnu11 -DINTEL -DCACHE_LINE_SIZE=$ENV{CACHE_LINE_SIZE}"
)

target_include_directories(capq PUBLIC
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library
    ${CMAKE_SOURCE_DIR}/include/Galois/WorkList/capq/capq_c_library/qdlocksrc
)
