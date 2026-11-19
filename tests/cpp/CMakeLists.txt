include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(oehealth_tests
    test_context.cpp
    test_oehealth.cpp
    test_result.cpp
)

# Tests reach private headers under src/ (reference tables, variant
# selection, helpers) as well as the public API.
target_include_directories(oehealth_tests
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_definitions(oehealth_tests
    PRIVATE
        OEHEALTH_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data"
)

target_link_libraries(oehealth_tests
    PRIVATE
        oehealth
        nlohmann_json::nlohmann_json
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oehealth_tests)
