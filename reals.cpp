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

#define TICKS_PER_SECOND (1)
#define SECONDS_OF_MEMORY (150)

#define WIDTH (160)
#define HEIGHT (120)
#define PIX_SIZE (0.01)

#define I_KEY (1)
#define I_FILE (2)
#define O_WINDOW (1)
#define O_FILE (2)

#define WINDOW_NAME "REALS"
#define OUTPUT_FILE_NAME "test.avi"
#define CONFIG_FILE "data.csv"
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
        std::vector<cl_float> h_top_speeds(num_objects);
        std::vector<cl_float> h_wavelengths(num_objects*history_length);
        std::vector<cl_int> h_deprecated(num_objects*history_length);
        std::vector<cl_int> h_is_black_hole(num_objects);
        std::vector<cl_int> h_is_sphere(num_objects);
        std::vector<cl_int> h_boolean_constants(7, false);
        std::vector<cl_uchar3> h_colors(WIDTH*HEIGHT, ((cl_uchar3) {0, 0, 0}));

        /*std::string line;
        std::ifstream myfile (CONFIG_FILE);
        if (myfile.is_open()) {
            myfile >> line;
            myfile.close();
        }*/

        float cur_time = 0.0;
        unsigned int start_tick = 0;
        unsigned int end_tick = 100;

        for (int i = 0; i < end_tick+10; i++) {
            h_positions[2*i+0] = (cl_float3) {0.0, 0.0, 0.0};
            h_positions[2*i+1] = (cl_float3) {3.0, 0.0, 0.0};
            h_velocities[2*i+0] = (cl_float3) {0.0, 0.0, 0.0};
            h_velocities[2*i+1] = (cl_float3) {0.0, 0.0, 0.0};
            h_orientations_r[2*i+0] = (cl_float3) {0.0, -1.0, 0.0};
            h_orientations_r[2*i+1] = (cl_float3) {1.0, 0.0, 0.0};
            h_orientations_f[2*i+0] = (cl_float3) {1.0, 0.0, 0.0};
            h_orientations_f[2*i+1] = (cl_float3) {0.0, 1.0, 0.0};
            h_orientations_u[2*i+0] = (cl_float3) {0.0, 0.0, 1.0};
            h_orientations_u[2*i+1] = (cl_float3) {0.0, 0.0, 1.0};
            h_local_times[2*i+0] = 1.0;
            h_local_times[2*i+1] = 1.0;
            h_masses[2*i+0] = 0.0;
            h_masses[2*i+1] = 100.0;
            h_wavelengths[2*i+0] = 450.0;
            h_wavelengths[2*i+1] = 500.0;
            h_deprecated[2*i+0] = 0;
            h_deprecated[2*i+1] = 0;
        }
        h_optical_radii[0] = 0.0;
        h_optical_radii[1] = 0.5;
        h_top_speeds[0] = 0.0;
        h_top_speeds[1] = 0.0;
        h_is_black_hole[0] = 0;
        h_is_black_hole[1] = 0;
        h_is_sphere[0] = 0;
        h_is_sphere[1] = 1;

        cl::Buffer d_positions;
        cl::Buffer d_velocities;
        cl::Buffer d_orientations_r;
        cl::Buffer d_orientations_f;
        cl::Buffer d_orientations_u;
        cl::Buffer d_local_times;
        cl::Buffer d_masses;
        cl::Buffer d_optical_radii;
        cl::Buffer d_top_speeds;
        cl::Buffer d_wavelengths;
        cl::Buffer d_deprecated;
        cl::Buffer d_is_black_hole;
        cl::Buffer d_is_sphere;
        cl::Buffer d_boolean_constants;
        cl::Buffer d_colors;

        int input;
        int INPUT_METHOD = I_KEY;
        int OUTPUT_METHOD = O_WINDOW;
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
                unsigned int, unsigned int,
                cl::Buffer,
                unsigned int, unsigned int, unsigned int,
                cl::Buffer
                >(program, "reals");

            d_positions = cl::Buffer(
                context,
                begin(h_positions),
                end(h_positions),
                true);
            d_velocities = cl::Buffer(
                context,
                begin(h_velocities),
                end(h_velocities),
                true);
            d_orientations_r = cl::Buffer(
                context,
                begin(h_orientations_r),
                end(h_orientations_r),
                true);
            d_orientations_f = cl::Buffer(
                context,
                begin(h_orientations_f),
                end(h_orientations_f),
                true);
            d_orientations_u = cl::Buffer(
                context,
                begin(h_orientations_u),
                end(h_orientations_u),
                true);
            d_local_times = cl::Buffer(
                context,
                begin(h_local_times),
                end(h_local_times),
                true);
            d_masses = cl::Buffer(
                context,
                begin(h_masses),
                end(h_masses),
                true);
            d_optical_radii = cl::Buffer(
                context,
                begin(h_optical_radii),
                end(h_optical_radii),
                true);
            d_top_speeds = cl::Buffer(
                context,
                begin(h_top_speeds),
                end(h_top_speeds),
                true);
            d_wavelengths = cl::Buffer(
                context,
                begin(h_wavelengths),
                end(h_wavelengths),
                true);
            d_deprecated = cl::Buffer(
                context,
                begin(h_deprecated),
                end(h_deprecated),
                true);
            d_is_black_hole = cl::Buffer(
                context,
                begin(h_is_black_hole),
                end(h_is_black_hole),
                true);
            d_is_sphere = cl::Buffer(
                context,
                begin(h_is_sphere),
                end(h_is_sphere),
                true);
            d_boolean_constants = cl::Buffer(
                context,
                begin(h_boolean_constants),
                end(h_boolean_constants),
                true);
            d_colors = cl::Buffer(
                context,
                CL_MEM_WRITE_ONLY,
                sizeof(cl_uchar3)*WIDTH*HEIGHT);

            reals(
                    cl::EnqueueArgs(
                        queue,
                        cl::NDRange(HEIGHT, WIDTH)),
                    d_positions,
                    d_velocities,
                    d_orientations_r,
                    d_orientations_f,
                    d_orientations_u,
                    d_local_times,
                    d_masses,
                    d_optical_radii,
                    d_top_speeds,
                    d_wavelengths,
                    d_deprecated,
                    d_is_black_hole,
                    d_is_sphere,
                    cur_time,
                    start_tick,
                    end_tick,
                    d_boolean_constants,
                    history_length,
                    WIDTH,
                    HEIGHT,
                    d_colors);

            queue.finish();

            cl::copy(queue, d_colors, begin(h_colors), end(h_colors));

            for (int i = 0; i < HEIGHT; i++) {
                for (int j = 0; j < WIDTH; j++) {
                    cv::Vec3b color;
                    color[0] = h_colors[i*WIDTH+j].s[2];
                    color[1] = h_colors[i*WIDTH+j].s[1];
                    color[2] = h_colors[i*WIDTH+j].s[0];
                    frame.at<cv::Vec3b>(cv::Point(j,i)) = color;
                }
            }

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
