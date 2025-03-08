file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(unZip ALL
    COMMAND tar xvfz data/mars/mars_data.tar.gz -C data/mars
    COMMAND gunzip -c data/gato/diagnostics.dat.gz > data/gato/diagnostics.dat
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(example example.cpp)
target_link_libraries(example PRIVATE fusionio::fusionio)

add_executable(example_c example_c.c)
target_link_libraries(example_c PRIVATE fusionio::fusionio)

add_executable(example_fortran example_fortran.f90)
target_link_libraries(example_fortran PRIVATE fusionio::fusionio)

add_executable(example_push example_push.f90)
target_link_libraries(example_push PRIVATE fusionio::push)

if (FUSIONIO_ENABLE_PCMS)
    add_executable(example_pcms example_pcms.cpp)
    target_link_libraries(example_pcms PRIVATE fusionio::fusionio pcms::core)
    # add_test(pcms_example_m3dc1 ./example_pcms m3dc1)
    # add_test(pcms_example_geqdsk ./example_pcms geqdsk)
    # add_test(pcms_example_gpec ./example_pcms gpec)
endif()

add_test(cpp_example_m3dc1 ./example m3dc1)
add_test(cpp_example_geqdsk ./example geqdsk)
add_test(cpp_example_gpec ./example gpec)

add_test(c_example ./example_c)
add_test(fortran_example ./example_fortran)