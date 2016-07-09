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

#define I_KEY (1)
#define I_FILE (2)
#define O_WINDOW (1)
#define O_FILE (2)

void init_output(int OUTPUT_METHOD) {
        if (OUTPUT_METHOD & O_WINDOW) {
        }
        if (OUTPUT_METHOD & O_FILE) {
        }
}

int fetch_input(int method) {
        if (method & I_KEY) {
                // return cv::waitKey(1000/TICKS_PER_SECOND);
                return (-1);
                // TODO: Uncomment upon finishing of `init_output`.
        } else if (method & I_FILE) {
                return (-1);
                // TODO: Make this read a file.
        } else {
                return (-1);
        }
}

int main(int argc, char** argv) {
        int num_objects = 10;

        std::vector<float> h_positions(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY*3);
        std::vector<float> h_velocities(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY*3);
        std::vector<float> h_orientation_r(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY*3);
        std::vector<float> h_orientation_f(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY*3);
        std::vector<float> h_orientation_u(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY*3);
        std::vector<float> h_local_time(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<float> h_masses(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<int> h_deprecated(num_objects);

        int input;
        int INPUT_METHOD = I_KEY;
        int OUTPUT_METHOD = O_WINDOW | O_FILE;
        if (INPUT_METHOD & I_KEY) {
                OUTPUT_METHOD |= O_WINDOW;
        }

        init_output(OUTPUT_METHOD);

        while ((input = fetch_input(INPUT_METHOD)) != -1) {

        }
        return(0);
}
