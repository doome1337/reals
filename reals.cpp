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
#define INPUT_FILE_NAME "input.txt" 
#define OUTPUT_FILE_NAME "test.avi"

#define EXIT_CODE (81)

ifstream in_file;

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
uchar cur_com = 0;
int cur_dur = 0;

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
		if (cur_dur) {
			cur_dur--;
			return cur_com;
		}
		if (!in_file) {
			return -1;
		}
	        string line;
		in_file >> line;
		if (line.length() == 0) {
			return 0;
		}
		if (line.length() == 1) {
			return line.at(0);
		}
		cur_com = line.at(0);
		cur_dur = std::stoi(line.substr(1), NULL)-1;
		return cur_com;
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

        std::vector<cl_float3> h_positions(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_velocities(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_r(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_f(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float3> h_orientation_u(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float> h_local_time(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_float> h_masses(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        std::vector<cl_int> h_deprecated(num_objects);

        int input;
        int INPUT_METHOD = I_FILE;
        int OUTPUT_METHOD = O_FILE;
        if (INPUT_METHOD & I_KEY) {
                OUTPUT_METHOD |= O_WINDOW;
        }
	if (INPUT_METHOD & I_FILE) {
            in_file = ifstream(INPUT_FILE);
	    if (!in_file.is_open()) {
                INPUT_METHOD ^= I_FILE;
	    }
	}
	while((input = fetch_input(INPUT_METHOD)) != -1) {
		std::cout << input << std::endl;
	}

        cv::VideoWriter video_output;

        init_output(OUTPUT_METHOD, &video_output);

        cv::Mat frame(HEIGHT, WIDTH, CV_8UC3, cv::Scalar(0,0,0));

        while ((input = fetch_input(INPUT_METHOD)) != -1) {
                output_image(OUTPUT_METHOD, frame, &video_output);
        }
	if (in_file.is_open()) {
            in_file.close();
	}
        return(0);
}
