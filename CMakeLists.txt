cmake_minimum_required(VERSION 3.14)
project(lbm_sim CXX)

# ---------------------------------------------------------------------------
# Always emit compile_commands.json (for clangd / Neovim)
# ---------------------------------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# C++ standard
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Default build type -> Debug if none specified
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")

# Debug -> -O0, full debug symbols, no optimization (best for gdb / lldb)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra")

# Release -> -O3, no debug symbols, NDEBUG defined (disables assert())
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---------------------------------------------------------------------------
# Dependencies via FetchContent (raylib + glm) — no system install needed
# ---------------------------------------------------------------------------
include(FetchContent)

# --- raylib -----------------------------------------------------------------
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build raylib's own examples
set(BUILD_GAMES    OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG        5.5
    GIT_SHALLOW    TRUE
)

FetchContent_MakeAvailable(raylib)

# ---------------------------------------------------------------------------
# Your executable
# ---------------------------------------------------------------------------
add_executable(lbm_sim
    src/main.cpp
    # add more .cpp files here as the project grows
)

target_include_directories(lbm_sim PRIVATE src)

target_link_libraries(lbm_sim PRIVATE
    raylib
)

# raylib on Linux needs these system libs; FetchContent's raylib target
# usually propagates them itself, but this is here in case your platform
# doesn't pull them in automatically.
if(UNIX AND NOT APPLE)
    target_link_libraries(lbm_sim PRIVATE GL m pthread dl rt X11)
endif()
