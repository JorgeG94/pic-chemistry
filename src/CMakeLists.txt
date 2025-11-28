add_subdirectory(basis)
add_subdirectory(cli)
add_subdirectory(core)
add_subdirectory(fragmentation)
add_subdirectory(io)
add_subdirectory(parallel)

target_sources(${main_lib} PRIVATE mqc_driver.f90)
