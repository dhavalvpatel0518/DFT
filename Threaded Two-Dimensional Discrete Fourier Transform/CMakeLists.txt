cmake_minimum_required(VERSION 2.8)
project(P2)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
        Complex.cc
        Complex.h
        InputImage.cc
        InputImage.h
        fft2d.cc
        )

find_package (Threads)
add_executable(P2
        fft2d.cc
        Complex.cc
        Complex.h
        InputImage.cc
        InputImage.h)
target_link_libraries (P2 ${CMAKE_THREAD_LIBS_INIT})

