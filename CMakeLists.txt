# Copyright 2012-2014 Axel Huebl
#
# This file is part of MDSim.
#
# MDSim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MDSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MDSim.
# If not, see <http://www.gnu.org/licenses/>.
#

################################################################################
# Required CMake Version
################################################################################

cmake_minimum_required(VERSION 2.8.5)


################################################################################
# Pre-Project (Set Compiler Here)
################################################################################



################################################################################
# Project 
################################################################################

project(MDSim)

# additional pathes to find libraries
set(CMAKE_PREFIX_PATH
    "/usr/lib/x86_64-linux-gnu/"
    "$ENV{MPI_ROOT}"
   )


################################################################################
# Find MPI
################################################################################

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_C_LIBRARIES})

# Some MPI implementations' C libs depend on their C++ counter part
if(MPI_CXX_FOUND)
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)


################################################################################
# Compile & Link MDSim
################################################################################

file(GLOB SRCFILES "*.cpp")

add_executable(mdsim
    ${SRCFILES}
)

target_link_libraries(mdsim ${LIBS} m)


################################################################################
# Install MDSim
################################################################################

install(TARGETS mdsim
        RUNTIME DESTINATION bin)
