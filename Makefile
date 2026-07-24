##################################################################
##################################################################
# User settings: modify these for your use case
##################################################################
##################################################################
# Host compiler choice (needs support for C++20)
HOST_COMPILER ?= g++


##################################################################
##################################################################
# Non-user-settings: You probably won't need to change these.
##################################################################
##################################################################
# File names and file paths for the program
program_NAME := main
src_DIR := src
program_C_SRCS := $(wildcard $(src_DIR)/*.c)
program_CXX_SRCS := $(wildcard $(src_DIR)/*.cpp)
program_H_SRCS := $(wildcard $(src_DIR)/*.h)
program_HPP_SRCS := $(wildcard $(src_DIR)/*.hpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_CXX_ASMS := ${program_CXX_SRCS:.cpp=.s}

program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := "../../BHQuasinormalModes/QuasiNormalModes/external"
program_LIBRARY_DIRS :=
program_LIBRARIES := m dl quadmath openblas lapacke lapack


# Compiler flags
CXXFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CXXFLAGS += -std=c++20 -Wall
#-fext-numeric-literals  	#-DEIGEN_HAS_CONSTEXPR=1 #-DEIGEN_NO_DEBUG
CXXFLAGS += -march=native -pthread -fopenmp
#CXXFLAGS += -march=native -pthread
CXXFLAGS += -O3 -ffast-math
#CXXFLAGS += -g -fno-omit-frame-pointer -fext-numeric-literals
CXXFLAGS += -DNDEBUG


# Add linker flags
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir)) 
LDLIBS += $(foreach library,$(program_LIBRARIES),-l$(library))


.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDLIBS)

$(program_OBJS): $(program_H_SRCS) $(program_HPP_SRCS) $(program_CUH_SRCS) $(program_GEN_SRCS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -S -fverbose-asm $< -o $@

asm: $(program_CXX_ASMS)

clean:
	$(RM) $(program_NAME)
	$(RM) $(program_OBJS)
	$(RM) $(program_CXX_ASMS)
	$(RM) $(wildcard *~)
	$(RM) -r html latex

distclean: clean

show:
	echo $(CXX)
	echo $(GXX)
	echo $(GCC)
	echo $(LINK.cc)
	echo $(CC)
	echo $(CPP)
	echo $(RM)
	echo $(CXXFLAGS)
	echo $(NVCC)
	echo $(program_CXX_SRCS) "\n"
	echo $(program_HPP_SRCS) "\n"
	echo $(program_CXX_OBJS) "\n"
	echo $(program_OBJS) "\n"
	echo $(program_CU_SRCS) "\n"
	echo $(program_CUH_SRCS) "\n"
	echo $(program_CU_OBJS) "\n"
	echo $(device_link_OBJ) "\n"
	echo $(program_CXX_ASMS) "\n"
