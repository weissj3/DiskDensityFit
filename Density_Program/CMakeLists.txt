project(DiskFit)

cmake_minimum_required(VERSION 2.6)

#Broke out the code into a testable library

option(BUILD_TESTS "Build test programs" 0)

include_directories(${PROJECT_SOURCE_DIR}/includes)
include_directories(${PROJECT_SOURCE_DIR}/../tao)
include_directories(${PROJECT_SOURCE_DIR}/../undvc_common)

set(INCLUDES    includes/Density_to_Star_Counts.h)

set(SOURCE     src/Density_to_StarCounts_TrimmedDown.cpp
                src/Density_to_StarCounts_main.cpp)
                
if(${TEST})
    include_directories(${PROJECT_SOURCE_DIR}/Tests/includes)
    set(INCLUDES    ${INCLUDE}
                    Tests/includes/Density_toStarCounts_Tests.h)

    set(SOURCE     ${SOURCE}
                    Tests/src/Tests_main.cpp)
                
endif(${TEST})

link_directories(${Boost_LIBRARY_DIRS})

add_executable(${PROJECT_NAME} ${SOURCE})

if(${TEST})

    target_compile_definitions(${PROJECT_NAME} PUBLIC TEST)
    
endif(${TEST})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} asynchronous_algorithms synchronous_algorithms tao_util undvc_common ${Boost_LIBRARIES})


