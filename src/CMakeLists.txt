cmake_minimum_required(VERSION 2.8.4)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
IF(UNIX)
 link_libraries(m)
ENDIF(UNIX)

project (smash-contigs) 
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_BUILD_TYPE Release)

add_executable (smash-map smash-map.c mem.c seq.c time.c pos.c msg.c parser.c common.c buffer.c hash.c rmodel.c defs.h param.h) 
add_executable (smash-reduce smash-reduce.c mem.c lines.c time.c msg.c common.c defs.h param.h) 
add_executable (smash-visual smash-visual.c mem.c time.c msg.c paint.c common.c defs.h param.h)
TARGET_LINK_LIBRARIES(smash-map pthread)
TARGET_LINK_LIBRARIES(smash-reduce)
TARGET_LINK_LIBRARIES(smash-visual)
