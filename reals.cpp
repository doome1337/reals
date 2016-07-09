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

#define WIDTH (160)
#define HEIGHT (120)

#define I_KEY (1)
#define I_FILE (2)
#define O_WINDOW (1)
#define O_FILE (2)

#define WINDOW_NAME "REALS"
#define OUTPUT_FILE_NAME "test.avi"

#define EXIT_CODE (81)

void init_output(int output_method, cv::VideoWriter* video_output) {
        if (output_method & O_WINDOW) {
                cv::namedWindow(WINDOW_NAME, cv::WINDOW_AUTOSIZE);
        }
        if (output_method & O_FILE) {
                video_output->open(
                        OUTPUT_FILE_NAME,
                        CV_FOURCC('F','M','P','4'),
                        TICKS_PER_SECOND,
                        cv::Size(WIDTH, HEIGHT),
                        true);
        }
}

int fetch_input(int input_method) {
        if (input_method & I_KEY) {
                int input = cv::waitKey(1000/TICKS_PER_SECOND);
                if (input == -1) {
                        return 0;
                }
                if (input == EXIT_CODE) {
                        return -1;
                }
                return input;
        } else if (input_method & I_FILE) {
                return (-1);
                // TODO: Make this read a file.
        } else {
                return (-1);
        }
}

void output_image(int output_method, cv::Mat frame, cv::VideoWriter* video_output) {
        if (output_method & O_WINDOW) {
                cv::imshow(WINDOW_NAME, frame);
        }
        if (output_method & O_FILE) {
                (*video_output) << frame;
        }
}

int main(int argc, char** argv) {
        int num_objects = 10;

        const unsigned int history_length = TICKS_PER_SECOND * SECONDS_OF_MEMORY;
        std::vector<cl_float3> h_positions(num_objects*history_length);
        std::vector<cl_float3> h_velocities(num_objects*history_length);
        std::vector<cl_float3> h_orientations_r(num_objects*history_length);
        std::vector<cl_float3> h_orientations_f(num_objects*history_length);
        std::vector<cl_float3> h_orientations_u(num_objects*history_length);
        std::vector<cl_float> h_local_times(num_objects*history_length);
        std::vector<cl_float> h_masses(num_objects*history_length);
        std::vector<cl_float> h_optical_radii(num_objects);
        std::vector<cl_int> h_deprecated(num_objects);
        int start_tick = 0;
        int end_tick = 0;

        cl::Buffer d_positions;
        cl::Buffer d_velocities;
        cl::Buffer d_orientations_r;
        cl::Buffer d_orientations_f;
        cl::Buffer d_orientations_u;
        cl::Buffer d_local_times;
        cl::Buffer d_masses;
        cl::Buffer d_optical_radii;
        cl::Buffer d_deprecated;

        int input;
        int INPUT_METHOD = I_KEY;
        int OUTPUT_METHOD = O_WINDOW | O_FILE;
        if (INPUT_METHOD & I_KEY) {
                OUTPUT_METHOD |= O_WINDOW;
        }

        cv::VideoWriter video_output;

        init_output(OUTPUT_METHOD, &video_output);

        cv::Mat frame(HEIGHT, WIDTH, CV_8UC3, cv::Scalar(0,0,0));

        try 
        {
            // Create a context
            cl::Context context(DEVICE);

            // Load in kernel source, creating a program object for the context

            cl::Program program(context, util::loadProgram("reals.cl"), true);

            // Get the command queue
            cl::CommandQueue queue(context);

            // Create the kernel functor

            auto reals = cl::make_kernel<
                cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
                cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
                cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
                float,
                unsigned int, unsigned int, unsigned int, unsigned int,
                cl::Buffer, cl::Buffer,
                unsigned int, unsigned int,
                float>(program, "reals");

            /*d_colors  = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*3*WIDTH*HEIGHT);

            sim(
                    cl::EnqueueArgs(
                        queue,
                        cl::NDRange(HEIGHT, WIDTH)), 
                    d_colors,
                    HEIGHT,
                    WIDTH,
                    DISTANCE,
                    SCREEN_DIST2,
                    PIX_HEIGHT,
                    PIX_WIDTH);

            queue.finish();

            cl::copy(queue, d_colors, begin(h_colors), end(h_colors));*/

            while ((input = fetch_input(INPUT_METHOD)) != -1) {
                output_image(OUTPUT_METHOD, frame, &video_output);
            }
            return 0;
        }
        catch (cl::Error err) {
            std::cout << "Exception\n";
            std::cerr 
                << "ERROR: "
                << err.what()
                << "("
                << err_code(err.err())
                << ")"
                << std::endl;
        }
        return(0);
}
