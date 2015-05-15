cmake_minimum_required(VERSION 3.2)
project(esmf_proj4)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(SOURCE_FILES
##    test/test_transform.cpp
##    test/test_transform_data.h
#    doc/examples/transform_example.cpp
#    src/transform.cpp
#    )

# DEPENDENCIES #########################################################################################################

find_path(PROJ_API_H proj_api.h)
find_library(PROJ_LIB proj)

include_directories(include ${PROJ_API_H})
#link_directories(${PROJ_LIB})
link_libraries(${PROJ_LIB})

# SHARE ################################################################################################################

add_library(transform SHARED src/transform.cpp)
#target_include_directories (transform PUBLIC include)
#target_link_libraries(transform PUBLIC ${PROJ_LIB})

link_libraries(transform)

# EXE ##################################################################################################################

#add_executable(esmf_proj4 ${SOURCE_FILES})
add_executable(run_transform_example doc/examples/transform_example.cpp)
#target_link_libraries(run_transform_example transform)
add_executable(run_transform_tests test/test_transform.cpp)
#target_link_libraries(run_transform_tests transform)
add_executable(esmf_proj4 test/test_transform.cpp)
#add_executable(esmf_proj4 doc/examples/transform_example.cpp)
#target_link_libraries(esmf_proj4 transform)
#target_include_directories (esmf_proj4 PUBLIC include)
#target_link_libraries(esmf_proj4 proj)

install(TARGETS transform
#        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

#install(TARGETS run_transform_example run_transform_tests
#        RUNTIME DESTINATION bin)

install(FILES include/transform.h
        DESTINATION include)

# TESTING ##############################################################################################################

enable_testing()

add_test(test_transform run_transform_tests)
add_test(test_transform_example run_transform_example)