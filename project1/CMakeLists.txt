cmake_minimum_required(VERSION 3.0)
project(my_project)

# set(CMAKE_CXX_STANDARD 23)

# 设置OpenCV依赖项
find_package(OpenCV REQUIRED)

# 添加可执行文件
add_executable(main src/main.cpp)
add_executable(cys src/cys.cpp)

# 链接OpenCV库
target_link_libraries(main ${OpenCV_LIBS})
target_link_libraries(cys ${OpenCV_LIBS})
