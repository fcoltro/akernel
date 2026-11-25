# ak_workbench: a minimal desktop app to try the kernel on real geometry.
#
# Dependencies are fetched at configure time (needs network once):
#   GLFW 3.4        window, input and OpenGL context
#   Dear ImGui      buttons and panels (GLFW + OpenGL 2 backends), needs git
# The viewport uses fixed-function OpenGL 1.1/2.1, so no GL loader is needed:
# opengl32 on Windows, libGL on Linux.
include(FetchContent)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
# X11 only on Linux (runs under XWayland too), so no wayland-scanner is needed.
set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "")

FetchContent_Declare(glfw
  URL https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.zip
  URL_HASH SHA256=b5ec004b2712fd08e8861dc271428f048775200a2df719ccf575143ba749a3e9)
# Dear ImGui publishes no release archives; a shallow clone of the tag (pinned
# by its commit below for reference: 97428e8ac99e339ce05eee531cf55b77b29ea709).
FetchContent_Declare(imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  GIT_TAG v1.91.9
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(glfw imgui)

find_package(OpenGL REQUIRED)

add_library(ak_imgui STATIC
  ${imgui_SOURCE_DIR}/imgui.cpp
  ${imgui_SOURCE_DIR}/imgui_draw.cpp
  ${imgui_SOURCE_DIR}/imgui_tables.cpp
  ${imgui_SOURCE_DIR}/imgui_widgets.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp)
target_include_directories(ak_imgui SYSTEM PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(ak_imgui PUBLIC glfw OpenGL::GL)

add_executable(ak_workbench main.cpp)
target_link_libraries(ak_workbench PRIVATE akernel ak_imgui)
if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(ak_workbench PRIVATE Threads::Threads)
endif()
