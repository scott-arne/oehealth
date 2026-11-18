include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(oehealth_tests
    test_oehealth.cpp
)

target_link_libraries(oehealth_tests
    PRIVATE
        oehealth
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oehealth_tests)
