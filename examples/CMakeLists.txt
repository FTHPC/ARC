
find_package(LibPressio REQUIRED)
find_package(std_compat REQUIRED)
find_package(SZ REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)
pkg_search_module(ZSTD IMPORTED_TARGET GLOBAL libzstd)
pkg_search_module(FFTW3 QUIET IMPORTED_TARGET GLOBAL fftw3)

add_executable(arc_pressio_example arc_pressio_example.c)
target_link_libraries(arc_pressio_example PRIVATE arc LibPressio::libpressio SZ)

add_executable(libpressio_example_sz libpressio_example_sz.c)
target_link_libraries(libpressio_example_sz PRIVATE arc LibPressio::libpressio SZ)

add_executable(libpressio_example_zfp libpressio_example_zfp.c)
target_link_libraries(libpressio_example_zfp PRIVATE arc LibPressio::libpressio SZ)
