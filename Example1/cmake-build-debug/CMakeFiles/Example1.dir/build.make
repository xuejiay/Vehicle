# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xuejiaoyang/Documents/github/Example1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Example1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Example1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Example1.dir/flags.make

CMakeFiles/Example1.dir/main.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Example1.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/main.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/main.cpp

CMakeFiles/Example1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/main.cpp > CMakeFiles/Example1.dir/main.cpp.i

CMakeFiles/Example1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/main.cpp -o CMakeFiles/Example1.dir/main.cpp.s

CMakeFiles/Example1.dir/DI_CST.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/DI_CST.cpp.o: ../DI_CST.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Example1.dir/DI_CST.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/DI_CST.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/DI_CST.cpp

CMakeFiles/Example1.dir/DI_CST.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/DI_CST.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/DI_CST.cpp > CMakeFiles/Example1.dir/DI_CST.cpp.i

CMakeFiles/Example1.dir/DI_CST.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/DI_CST.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/DI_CST.cpp -o CMakeFiles/Example1.dir/DI_CST.cpp.s

CMakeFiles/Example1.dir/DI_CST_advance.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/DI_CST_advance.cpp.o: ../DI_CST_advance.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Example1.dir/DI_CST_advance.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/DI_CST_advance.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/DI_CST_advance.cpp

CMakeFiles/Example1.dir/DI_CST_advance.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/DI_CST_advance.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/DI_CST_advance.cpp > CMakeFiles/Example1.dir/DI_CST_advance.cpp.i

CMakeFiles/Example1.dir/DI_CST_advance.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/DI_CST_advance.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/DI_CST_advance.cpp -o CMakeFiles/Example1.dir/DI_CST_advance.cpp.s

CMakeFiles/Example1.dir/SDI.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/SDI.cpp.o: ../SDI.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Example1.dir/SDI.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/SDI.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/SDI.cpp

CMakeFiles/Example1.dir/SDI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/SDI.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/SDI.cpp > CMakeFiles/Example1.dir/SDI.cpp.i

CMakeFiles/Example1.dir/SDI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/SDI.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/SDI.cpp -o CMakeFiles/Example1.dir/SDI.cpp.s

CMakeFiles/Example1.dir/rhs.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/rhs.cpp.o: ../rhs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Example1.dir/rhs.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/rhs.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/rhs.cpp

CMakeFiles/Example1.dir/rhs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/rhs.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/rhs.cpp > CMakeFiles/Example1.dir/rhs.cpp.i

CMakeFiles/Example1.dir/rhs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/rhs.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/rhs.cpp -o CMakeFiles/Example1.dir/rhs.cpp.s

CMakeFiles/Example1.dir/trajSampling.cpp.o: CMakeFiles/Example1.dir/flags.make
CMakeFiles/Example1.dir/trajSampling.cpp.o: ../trajSampling.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/Example1.dir/trajSampling.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example1.dir/trajSampling.cpp.o -c /Users/xuejiaoyang/Documents/github/Example1/trajSampling.cpp

CMakeFiles/Example1.dir/trajSampling.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example1.dir/trajSampling.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xuejiaoyang/Documents/github/Example1/trajSampling.cpp > CMakeFiles/Example1.dir/trajSampling.cpp.i

CMakeFiles/Example1.dir/trajSampling.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example1.dir/trajSampling.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xuejiaoyang/Documents/github/Example1/trajSampling.cpp -o CMakeFiles/Example1.dir/trajSampling.cpp.s

# Object files for target Example1
Example1_OBJECTS = \
"CMakeFiles/Example1.dir/main.cpp.o" \
"CMakeFiles/Example1.dir/DI_CST.cpp.o" \
"CMakeFiles/Example1.dir/DI_CST_advance.cpp.o" \
"CMakeFiles/Example1.dir/SDI.cpp.o" \
"CMakeFiles/Example1.dir/rhs.cpp.o" \
"CMakeFiles/Example1.dir/trajSampling.cpp.o"

# External object files for target Example1
Example1_EXTERNAL_OBJECTS =

Example1: CMakeFiles/Example1.dir/main.cpp.o
Example1: CMakeFiles/Example1.dir/DI_CST.cpp.o
Example1: CMakeFiles/Example1.dir/DI_CST_advance.cpp.o
Example1: CMakeFiles/Example1.dir/SDI.cpp.o
Example1: CMakeFiles/Example1.dir/rhs.cpp.o
Example1: CMakeFiles/Example1.dir/trajSampling.cpp.o
Example1: CMakeFiles/Example1.dir/build.make
Example1: /Users/xuejiaoyang/sundial/instdir/lib/libsundials_cvode.a
Example1: /Users/xuejiaoyang/sundial/instdir/lib/libsundials_nvecserial.a
Example1: CMakeFiles/Example1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable Example1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Example1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Example1.dir/build: Example1

.PHONY : CMakeFiles/Example1.dir/build

CMakeFiles/Example1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Example1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Example1.dir/clean

CMakeFiles/Example1.dir/depend:
	cd /Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xuejiaoyang/Documents/github/Example1 /Users/xuejiaoyang/Documents/github/Example1 /Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug /Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug /Users/xuejiaoyang/Documents/github/Example1/cmake-build-debug/CMakeFiles/Example1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Example1.dir/depend

