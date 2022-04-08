cmake_minimum_required(VERSION 3.15)

project(dowser VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib/JUCE)

juce_add_console_app(dowser
    PRODUCT_NAME "dowser")

target_sources(dowser
    PRIVATE
        src/main.cpp)

target_include_directories(dowser PRIVATE src)

target_link_libraries(dowser
    PRIVATE
        juce::juce_core
        juce::juce_dsp
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_warning_flags)


# never mind this; juce/cmake ends up doing the same by default
# if(CMAKE_BUILD_TYPE MATCHES "Debug")
#     set(cflags_dowser "-O0 -g")
# elseif(CMAKE_BUILD_TYPE MATCHES "Release")
#     set(cflags_dowser "-O3")
# endif()

# set(set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${cflags_dowser} -Wall -Werror"))