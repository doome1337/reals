#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

#define __CL_ENABLE_EXCEPTIONS
#include <OpenCL/opencl.h>
#include <OpenCL/cl.h>

#include "CL_libs/cl.hpp"
#include "CL_libs/util.hpp"
#include "CL_libs/err_code.h"

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <fstream>

// pick up device type from compiler command line or from the default type
#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif

#define TICKS_PER_SECOND (60)
#define SECONDS_OF_MEMORY (150)

#define MAX_INTENSITY_AT_REST (255)
#define SPEED_OF_LIGHT (299792458.0)
// ^ Make this a command line argument.
#define VISIBLE_PARALLAX (512)
#define GRAVITATIONAL_CONSTANT (0.0000000000000667408)

int main(int argc, char** argv) {
        int num_objects = 10;

        std::vector<cl_float3> h_positions(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_velocities(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_r(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_f(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_u(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float> h_local_time(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float> h_masses(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_int> h_deprecated(num_objects);

        return(0);
}
