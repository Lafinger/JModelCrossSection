cmake_minimum_required(VERSION 3.16)
project(MeshSlicer)

# Libigl
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(LIBIGL_EMBREE           "Build target igl::embree"           ON)
option(LIBIGL_GLFW             "Build target igl::glfw"             ON)
option(LIBIGL_IMGUI            "Build target igl::imgui"            ON)
option(LIBIGL_OPENGL           "Build target igl::opengl"           ON)
option(LIBIGL_PNG              "Build target igl::png"              ON)
option(LIBIGL_PREDICATES       "Build target igl::predicates"       ON)
# option(LIBIGL_XML              "Build target igl::xml"              ON)
option(LIBIGL_COPYLEFT_CGAL    "Build target igl_copyleft::cgal"    ON)
option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" ON)
include(libigl)

include_directories(./glm)
include_directories(./utils)


add_executable(${PROJECT_NAME} main.cpp ./utils/MyQuaternion.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::embree igl::glfw igl::imgui igl::opengl igl::png igl::predicates igl_restricted::triangle igl::core)
