include(FetchContent)

# Add Google Test
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1  # Specify the version of GoogleTest
)
FetchContent_MakeAvailable(googletest)
