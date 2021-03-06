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

#define TICKS_PER_SECOND (30)
#define SECONDS_OF_MEMORY (9)

#define MAX_INTENSITY_AT_REST (255)
#define SPEED_OF_LIGHT (299792458.0)
// ^ Make this a command line argument.
#define VISIBLE_PARALLAX (512)
#define GRAVITATIONAL_CONSTANT (0.0000000000000667408)
#define GRAVITATIONAL_RATIO (40.0)
#define MASSIVE_BOUND (1000000000.0)
#define DECAY_PER_ORBIT (0.95)
#define DECAY_COEFFICIENT (0.2)
#define LIGHT_SLOWING_RATIO (4)

#define USE_INSTANT_GRAVITY (0)
#define APPLY_ORBIT_DECAY (0)
#define USE_ALTERNATE_FACTOR (0)

#define WIDTH (160)
#define HEIGHT (120)
#define PIX_SIZE (0.01)

#define I_KEY (1)
#define I_FILE (2)
#define O_WINDOW (1)
#define O_FILE (2)

#define WINDOW_NAME "REALS"
#define INPUT_FILE_NAME "input.txt" 
#define OUTPUT_FILE_NAME "test.avi"
#define CONFIG_FILE "data.csv"
#define EXIT_CODE (81)

std::ifstream in_file;

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

int num_objects = 2;

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
        std::string line;
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

void update();

std::vector<cl_float3> h_positions;
std::vector<cl_float3> h_velocities;
std::vector<cl_float3> h_orientations_r;
std::vector<cl_float3> h_orientations_f;
std::vector<cl_float3> h_orientations_u;
std::vector<cl_float> h_local_times;
std::vector<cl_float> h_masses;
std::vector<cl_float> h_top_speeds;
std::vector<cl_float> h_optical_radii;
std::vector<cl_int> h_deprecated;
int ticks;
int main(int argc, char** argv) {
        int num_objects = 2;
        const unsigned int history_length = TICKS_PER_SECOND * SECONDS_OF_MEMORY;

        std::vector<cl_float> h_wavelengths(num_objects*history_length);
        std::vector<cl_int> h_is_black_hole(num_objects);
        std::vector<cl_int> h_is_sphere(num_objects);
        std::vector<cl_int> h_boolean_constants(7, false);
        std::vector<cl_uchar3> h_colors(WIDTH*HEIGHT, ((cl_uchar3) {0, 0, 0}));

        ticks = 100;

        h_positions = std::vector<cl_float3>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_velocities = std::vector<cl_float3>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_orientations_r = std::vector<cl_float3>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_orientations_f = std::vector<cl_float3>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_orientations_u = std::vector<cl_float3>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_local_times = std::vector<cl_float>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_masses = std::vector<cl_float>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);
        h_top_speeds = std::vector<cl_float>(num_objects);
        h_optical_radii = std::vector<cl_float>(num_objects);
        h_deprecated = std::vector<cl_int>(num_objects*TICKS_PER_SECOND*SECONDS_OF_MEMORY);

        /*std::string line;
        std::ifstream myfile (CONFIG_FILE);
        if (myfile.is_open()) {
            myfile >> line;
            myfile.close();
        }*/

        float cur_time = 0.0;
        unsigned int start_tick = 0;
        unsigned int end_tick = ticks;

        for (int i = 0; i < history_length-1; i++) {
            h_positions.at(2*i+0) = (cl_float3) {0.0, 0.0, 0.0};
            h_positions.at(2*i+1) = (cl_float3) {2.0, 0.0, 0.0};
            h_velocities.at(2*i+0) = (cl_float3) {0.0, 0.0, 0.0};
            h_velocities.at(2*i+1) = (cl_float3) {0.0, 0.0, 0.0};
            h_orientations_r.at(2*i+0) = (cl_float3) {0.0, -1.0, 0.0};
            h_orientations_r.at(2*i+1) = (cl_float3) {1.0, 0.0, 0.0};
            h_orientations_f.at(2*i+0) = (cl_float3) {1.0, 0.0, 0.0};
            h_orientations_f.at(2*i+1) = (cl_float3) {0.0, 1.0, 0.0};
            h_orientations_u.at(2*i+0) = (cl_float3) {0.0, 0.0, 1.0};
            h_orientations_u.at(2*i+1) = (cl_float3) {0.0, 0.0, 1.0};
            h_local_times.at(2*i+0) = 1.0;
            h_local_times.at(2*i+1) = 1.0;
            h_masses.at(2*i+0) = 0.0;
            h_masses.at(2*i+1) = 100.0;
            h_wavelengths.at(2*i+0) = 450.0;
            h_wavelengths.at(2*i+1) = 500.0;
            h_deprecated.at(2*i+0) = 0;
            h_deprecated.at(2*i+1) = 0;
        }
        h_optical_radii.at(0) = 0.0;
        h_optical_radii.at(1) = 0.5;
        h_top_speeds.at(0) = 0.005;
        h_top_speeds.at(1) = 0.005;
        h_is_black_hole.at(0) = 0;
        h_is_black_hole.at(1) = 0;
        h_is_sphere.at(0) = 0;
        h_is_sphere.at(1) = 1;

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
        int INPUT_METHOD = I_FILE;
        int OUTPUT_METHOD = O_FILE;
        if (INPUT_METHOD & I_KEY) {
                OUTPUT_METHOD |= O_WINDOW;
        }
	if (INPUT_METHOD & I_FILE) {
            in_file = std::ifstream(INPUT_FILE_NAME);
	    if (!in_file.is_open()) {
                INPUT_METHOD ^= I_FILE;
	    }
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

            while ((input = fetch_input(INPUT_METHOD)) != -1) {
                end_tick = ticks;
                cl::copy(queue, begin(h_positions), end(h_positions), d_positions);
                cl::copy(queue, begin(h_velocities), end(h_velocities), d_velocities);
                cl::copy(queue, begin(h_top_speeds), end(h_top_speeds), d_top_speeds);
                //cl::copy(queue, begin(h_velocities), end(h_velocities), d_velocities);

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

                queue.flush();

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
                output_image(OUTPUT_METHOD, frame, &video_output);
                update();
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
	if (in_file.is_open()) {
            in_file.close();
	}
        return(0);
}

int worldline(int min_time, int max_time, int time, int object_index, cl_float3 object_position);

cl_float factor(cl_float ratio);

int index_time_obj(int tick, int object, int num_objects) {
    return tick*num_objects + object;
}

cl_float3 vec_add(cl_float3 vec1, cl_float3 vec2) {
    return (cl_float3) {
        vec1.s[0] + vec2.s[0],
        vec1.s[1] + vec2.s[1],
        vec1.s[2] + vec2.s[2]};
}

cl_float3 vec_minus(cl_float3 vec1, cl_float3 vec2) {
    return (cl_float3) {
        vec1.s[0] - vec2.s[0],
        vec1.s[1] - vec2.s[1],
        vec1.s[2] - vec2.s[2]};
}

cl_float dot(cl_float3 vec1, cl_float3 vec2) {
    return (vec1.s[0]*vec2.s[0] + vec1.s[1]*vec2.s[1] + vec1.s[2]*vec2.s[2]);
}

cl_float length(cl_float3 vec) {
    return sqrt(dot(vec, vec));
}

cl_float schwarzschild_radius(cl_float mass) {
    return 2 * mass * GRAVITATIONAL_CONSTANT / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
}

cl_float3 scale(cl_float scalar, cl_float3 vec) {
    return (cl_float3) {scalar*vec.s[0], scalar*vec.s[1], scalar*vec.s[2]};
}

void update() {
    for (int i = 0; i < num_objects; i++) {
        cl_float3 new_position = h_positions[index_time_obj(ticks, i, num_objects)];
        cl_float3 new_velocity = h_velocities[index_time_obj(ticks, i, num_objects)];
        cl_float new_mass = h_masses[index_time_obj(ticks, i, num_objects)];
        cl_float new_times = h_local_times[index_time_obj(ticks, i, num_objects)];
        new_position = vec_add(new_position, h_velocities[index_time_obj(ticks, i, num_objects)]);
        cl_float combined_ratio = 0;
        for (int j = 0; j < num_objects; j++) {
            if (i == j) {
                continue;
            }
            cl_float3 relative_position = vec_minus(h_positions[index_time_obj(ticks, i, num_objects)],
                h_positions[index_time_obj(ticks, j, num_objects)]);
            if (length(relative_position) >
                GRAVITATIONAL_RATIO * schwarzschild_radius(h_masses[index_time_obj(ticks, j, num_objects)])) {
                continue;
            }
            cl_float wave_ticks = USE_INSTANT_GRAVITY ? ticks :
                worldline(
                    ticks - length(relative_position) / (SPEED_OF_LIGHT - h_top_speeds[j]),
                    ticks - length(relative_position) / (SPEED_OF_LIGHT + h_top_speeds[j]),
                    ticks, j, h_positions[index_time_obj(ticks, i, num_objects)]);
            if (!USE_INSTANT_GRAVITY) {
                relative_position = vec_minus(
                    h_positions[index_time_obj(ticks, i, num_objects)],
                    h_positions[index_time_obj(ticks, j, num_objects)]);
            }
            combined_ratio += schwarzschild_radius(h_masses[index_time_obj(wave_ticks, j, num_objects)]) /
                length(relative_position);
            if (APPLY_ORBIT_DECAY && h_masses[index_time_obj(ticks, i, num_objects)] > MASSIVE_BOUND &&
                h_masses[index_time_obj(wave_ticks, j, num_objects)] > MASSIVE_BOUND) {
                new_position = vec_add(
                    new_position,
                    scale(
                        0.5 * (1 - pow(DECAY_PER_ORBIT, DECAY_COEFFICIENT *
                    sqrt(h_masses[index_time_obj(wave_ticks, j, num_objects)] / length(relative_position) /
                    length(relative_position) / length(relative_position)))),
                        relative_position));
            }
            new_velocity = vec_add(
                new_velocity,
                scale(
                    GRAVITATIONAL_CONSTANT * h_masses[index_time_obj(wave_ticks, j, num_objects)]
                    / dot(relative_position, relative_position) / length(relative_position),
                    relative_position));
        }
        new_times += factor(combined_ratio) *
            sqrt(1 - length(h_velocities[index_time_obj(ticks, i, num_objects)]) *
            length(h_velocities[index_time_obj(ticks, i, num_objects)]) / SPEED_OF_LIGHT / SPEED_OF_LIGHT);
        h_positions[index_time_obj(ticks+1, i, num_objects)] = new_position;
        h_velocities[index_time_obj(ticks+1, i, num_objects)] = new_velocity;
        h_local_times[index_time_obj(ticks+1, i, num_objects)] = new_times;
        if (length(new_velocity) > h_top_speeds[i]) {
            h_top_speeds[i] = length(new_velocity) ;
        }
    }
        std::cout << h_positions[index_time_obj(ticks, 0, num_objects)].s[0] << std::endl;
        std::cout << h_positions[index_time_obj(ticks+1, 0, num_objects)].s[0] << std::endl;
        std::cout << h_positions[index_time_obj(ticks, 1, num_objects)].s[0] << std::endl;
        std::cout << h_positions[index_time_obj(ticks+1, 1, num_objects)].s[0] << std::endl;
    ticks++;
    for (int i = 0; i < num_objects; i++) {
        for (int j = 0; j < num_objects; j++) {
            if (i == j || h_masses[index_time_obj(ticks, i, num_objects)] < h_masses[index_time_obj(ticks, j, num_objects)]) {
                continue;
            }
            cl_float3 relative_position = vec_minus(
                h_positions[index_time_obj(ticks, i, num_objects)],
                h_positions[index_time_obj(ticks, j, num_objects)]);
            if (length(relative_position) < h_optical_radii[i]) {
                h_deprecated[j] = true;
                h_positions[index_time_obj(ticks, i, num_objects)] = vec_add(
                    h_positions[index_time_obj(ticks, i, num_objects)],
                    scale(h_masses[index_time_obj(ticks, j, num_objects)] /
                    (h_masses[index_time_obj(ticks, i, num_objects)] + h_masses[index_time_obj(ticks, j, num_objects)]),
                    relative_position));
                h_masses[index_time_obj(ticks, i, num_objects)] += h_masses[index_time_obj(ticks, j, num_objects)];
            }
        }
    }
    h_velocities[index_time_obj(ticks, 1, num_objects)] = vec_add(
        h_velocities[index_time_obj(ticks, 1, num_objects)],
        (cl_float3) {(cl_float) (0.0003*cos(0.2*ticks)),(cl_float) (0.0003*sin(0.2*ticks)), 0.0});
}

int worldline(int min_time, int max_time, int time, int object_index, cl_float3 object_position) {
    int lower_bound = fmax(min_time, 0);
    int upper_bound = fmin(max_time, TICKS_PER_SECOND * SECONDS_OF_MEMORY - 1);
    int iterations = TICKS_PER_SECOND * SECONDS_OF_MEMORY / 2;
    for (int i = 0; i < iterations; i++) {
        if (upper_bound == lower_bound) {
            return lower_bound;
        }
        int half_time = (lower_bound + upper_bound) / 2;
        if (SPEED_OF_LIGHT * (half_time - time) > length(vec_minus(h_positions[index_time_obj(ticks, object_index, num_objects)], object_position))) {
            upper_bound = half_time;
        } else {
            lower_bound = half_time;
        }
    }
    return lower_bound;
}

cl_float factor(cl_float ratio) {
    if (USE_ALTERNATE_FACTOR) {
        return (1 - LIGHT_SLOWING_RATIO * ratio) * (1 - LIGHT_SLOWING_RATIO * ratio);
    } else {
        return (1 - ratio) / (1 + ratio) / (1 + ratio) / (1 + ratio);
    }
}
