cmake_minimum_required(VERSION 3.15)

project(dowser VERSION 0.0.1)


set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib/JUCE)

juce_add_console_app(dowser
    PRODUCT_NAME "dowser")

target_sources(dowser
    PRIVATE
        main.cpp)

target_link_libraries(dowser
    PRIVATE
        juce::juce_core
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_warning_flags)