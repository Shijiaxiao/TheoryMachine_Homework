# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build

# Include any dependencies generated for this target.
include CMakeFiles/cys.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cys.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cys.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cys.dir/flags.make

CMakeFiles/cys.dir/src/cys.cpp.o: CMakeFiles/cys.dir/flags.make
CMakeFiles/cys.dir/src/cys.cpp.o: /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/src/cys.cpp
CMakeFiles/cys.dir/src/cys.cpp.o: CMakeFiles/cys.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cys.dir/src/cys.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/cys.dir/src/cys.cpp.o -MF CMakeFiles/cys.dir/src/cys.cpp.o.d -o CMakeFiles/cys.dir/src/cys.cpp.o -c /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/src/cys.cpp

CMakeFiles/cys.dir/src/cys.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cys.dir/src/cys.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/src/cys.cpp > CMakeFiles/cys.dir/src/cys.cpp.i

CMakeFiles/cys.dir/src/cys.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cys.dir/src/cys.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/src/cys.cpp -o CMakeFiles/cys.dir/src/cys.cpp.s

# Object files for target cys
cys_OBJECTS = \
"CMakeFiles/cys.dir/src/cys.cpp.o"

# External object files for target cys
cys_EXTERNAL_OBJECTS =

cys: CMakeFiles/cys.dir/src/cys.cpp.o
cys: CMakeFiles/cys.dir/build.make
cys: /opt/homebrew/lib/libopencv_gapi.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_stitching.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_alphamat.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_aruco.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_barcode.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_bgsegm.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_bioinspired.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_ccalib.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_dnn_objdetect.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_dnn_superres.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_dpm.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_face.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_freetype.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_fuzzy.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_hfs.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_img_hash.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_intensity_transform.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_line_descriptor.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_mcc.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_quality.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_rapid.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_reg.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_rgbd.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_saliency.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_sfm.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_stereo.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_structured_light.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_superres.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_surface_matching.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_tracking.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_videostab.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_viz.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_wechat_qrcode.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_xfeatures2d.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_xobjdetect.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_xphoto.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_shape.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_highgui.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_datasets.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_plot.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_text.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_ml.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_phase_unwrapping.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_optflow.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_ximgproc.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_video.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_videoio.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_imgcodecs.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_objdetect.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_calib3d.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_dnn.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_features2d.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_flann.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_photo.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_imgproc.4.7.0.dylib
cys: /opt/homebrew/lib/libopencv_core.4.7.0.dylib
cys: CMakeFiles/cys.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cys"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cys.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cys.dir/build: cys
.PHONY : CMakeFiles/cys.dir/build

CMakeFiles/cys.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cys.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cys.dir/clean

CMakeFiles/cys.dir/depend:
	cd /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1 /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1 /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build /Users/shijiaxiao/Desktop/cmake_cpp_projects/Theory_machine_project1/build/CMakeFiles/cys.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cys.dir/depend

