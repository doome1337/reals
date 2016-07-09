
ifndef CPPC
	CPPC=g++
endif

CPP_COMMON = CL_libs

CVFLAGS=-I/usr/local/include/opencv -I/usr/local/include -L/usr/local/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_core

CCFLAGS=-std=c++11

INC = -I $(CPP_COMMON)

LIBS = -lOpenCL -lrt

# Change this variable to specify the device type
# to the OpenCL device type of choice. You can also
# edit the variable in the source.
ifndef DEVICE
	DEVICE = CL_DEVICE_TYPE_DEFAULT
endif

# Check our platform and make sure we define the APPLE variable
# and set up the right compiler flags and libraries
PLATFORM = $(shell uname -s)
ifeq ($(PLATFORM), Darwin)
	CPPC = clang++
 	CCFLAGS += -stdlib=libc++
 	LIBS = -framework OpenCL
endif

CCFLAGS += -D DEVICE=$(DEVICE)

reals: reals.cpp
	$(CPPC) $^ $(INC) $(CCFLAGS) $(CVFLAGS) $(LIBS) -o $@


clean:
	rm -f vadd
