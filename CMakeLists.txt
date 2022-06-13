cmake_minimum_required(VERSION 3.15)

project(dowser VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib/JUCE)

juce_add_console_app(dowser
    PRODUCT_NAME "dowser")

target_sources(dowser
    PRIVATE
        src/main.cpp)

target_include_directories(dowser PRIVATE src lib)

target_link_libraries(dowser
    PRIVATE
        juce::juce_core
        juce::juce_dsp
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_warning_flags)

add_executable(test-watershed src/test_watershed.cpp src/output.cpp)
target_include_directories(test-watershed PRIVATE src)