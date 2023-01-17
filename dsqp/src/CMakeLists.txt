# MIT License
# Copyright (c) 2023 Goesta Stomberg, Henrik Ebel, Timm Faulwasser, Peter Eberhard
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required (VERSION 3.0)
# if(NOT CMAKE_VERSION VERSION_LESS 3.0)
#     cmake_policy(SET CMP0048 NEW)
# endif()
project (dsqp)

set(LCM_INCLUDE_DIR "/usr/local/include" CACHE FILEPATH "headers of lcm should be in the directory \${LCM_INCLUDE_DIR}/lcm")
set(LCM_LIBRARY_DIR "/usr/local/lib" CACHE FILEPATH "liblcm.so should be located in the directory \${LCM_LIBRARY_DIR}")
set(QPOASES_LIBRARY_DIR "/usr/local/lib" CACHE FILEPATH "libqpOASES.so should be located in the directory \${QPOASES_LIBRARY_DIR}")
set(QPOASES_INCLUDE_DIR "/usr/local/include" CACHE FILEPATH "headers of qpOASES should be in the directory \${QPOASES_LIBRARY_DIR} and its subfolder qpOASES")
#set(ENV{PKG_CONFIG_PATH} "/usr/local/lib/pkgconfig" CACHE FILEPATH "pkgconfig files for casadi and related packages need to be in the directory \${ENV{PKG_CONFIG_PATH}}, i.e., folder needs to contain casadi.pc, coinmetis.pc, conmumps.pc, ipopt.pc, opencv.pc") # pkgconfig files for casadi and related packages (folder needs to contain casadi.pc, coinmetis.pc, conmumps.pc, ipopt.pc, opencv.pc)

include_directories("${LCM_INCLUDE_DIR}")

# locate some libraries
find_library(LCM_LIB lcm ${LCM_LIBRARY_DIR})
find_library(QPOASES_LIB qpOASES ${QPOASES_LIBRARY_DIR} NO_DEFAULT_PATH)

# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(CASADI casadi)
if(NOT CASADI_LIBRARIES)
  message(SEND_ERROR "Could not find casadi.pc - did you set PKG_CONFIG_PATH?")
endif()

message("CASADI_LIBRARY_DIRS: ${CASADI_LIBRARY_DIRS}")
message("CASADI_LIBRARIES: ${CASADI_LIBRARIES}")
message("CASADI_INCLUDEDIR: ${CASADI_INCLUDEDIR}")
link_directories(${CASADI_LIBRARY_DIRS})


add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/lcm/itmessage/vector_t.hpp ${CMAKE_CURRENT_SOURCE_DIR}/lcm/itmessage/vector_idx_t.hpp
                   COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/lcm/build_all.sh 
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lcm/
                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lcm/vector_t.lcm ${CMAKE_CURRENT_SOURCE_DIR}/lcm/vector_idx_t.lcm)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
if ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	# set (CMAKE_CXX_FLAGS "-std=c++11 -Og -pthread -g -O0 -Wall -Wextra -pedantic -pedantic-errors -Werror=pedantic -Wcast-align -Wdisabled-optimization -Wformat=2  -Wlogical-op -Wmissing-include-dirs -Wstrict-overflow=2 -Wswitch-default -Wundef -fsanitize=address -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function")
	set (CMAKE_CXX_FLAGS "-std=c++11 -O3 -pthread -fno-math-errno")
endif()
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${QPOASES_INCLUDE_DIR}")

add_executable(dmpc_example dmpc_example.cpp handler_vector_idx_t.cpp lcm/itmessage/vector_t.hpp handler_stateVector.cpp lcm/itmessage/vector_idx_t.hpp sProb.cpp dsqpAgent.cpp admmAgent.cpp doptTimer.cpp dmpcAgent.cpp) 
target_link_libraries(dmpc_example ${LCM_LIB} ${QPOASES_LIB} casadi)
target_include_directories(dmpc_example PRIVATE ${CASADI_INCLUDEDIR})
