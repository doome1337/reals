// TODO LIST:
// Make ray_time an int
// Tick_index wrap in index_time_obj
#define PIX_SIZE (0.005)
#define GRAVITATIONAL_CONSTANT (0.0)
#define GRAVITATIONAL_RATIO (4000.0)
#define INTENSITY_FACTOR (1.0)
#define LIGHT_SLOWING_RATIO (4.0)
#define MASSIVE_BOUND (100000000.0)
#define SPEED_OF_LIGHT (0.05)
#define VISIBLE_PARALLAX (100.0)

#define NUM_OBJECTS (2)

#pragma OPENCL EXTENSION cl_intel_printf : enable

uchar3 ray_trace(
        // Histories
        __global float3* positions,
        __global float3* velocities,
        __global float3* orientations_r,
        __global float3* orientations_f,
        __global float3* orientations_u,
        __global float* masses,
        __global float* optical_radii,
        // Values
        __global float* top_speeds,
        __global float* wavelengths,
        __global int* is_black_hole,
        __global int* is_sphere,
        __global int* deprecated,
        float3 start_point,
        float3 start_ray,
        // Constants for this frame
        const float cur_time,
        const int num_objects,
        const int history_length,
        const int start_tick,
        const int end_tick,
        // Constant Settings
        const int APPLY_GR_RS,
        const int APPLY_INTENSITY,
        const int APPLY_LENSING,
        const int APPLY_SR_RS,
        const int USE_ALTERNATE_FACTOR,
        const int USE_INSTANT_GRAVITY,
        const int USE_LONG_STEP);

__kernel void reals (
        // Histories
        __global float3* positions,
        __global float3* velocities,
        __global float3* orientations_r,
        __global float3* orientations_f,
        __global float3* orientations_u,
        __global float* local_times,
        __global float* masses,
        __global float* optical_radii,
        // Values
        __global float* top_speeds,
        __global float* wavelengths,
        __global int* deprecated,
        __global int* is_black_hole,
        __global int* is_sphere,
        // Constants for this frame
        const float cur_time,
        const unsigned int start_tick,
        const unsigned int end_tick,
        // Constants for all frames
        __global int* boolean_constants,
        const unsigned int history_length,
        const unsigned int width,
        const unsigned int height,
        __global uchar3* colors) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if (i < height && j < width) {
        const int APPLY_GR_RS = boolean_constants[0];
        const int APPLY_INTENSITY = boolean_constants[1];
        const int APPLY_LENSING = boolean_constants[2];
        const int APPLY_SR_RS = boolean_constants[3];
        const int USE_ALTERNATE_FACTOR = boolean_constants[4];
        const int USE_INSTANT_GRAVITY = boolean_constants[5];
        const int USE_LONG_STEP = boolean_constants[6];
        colors[i*width + j] = ray_trace(
            // Histories
            positions,
            velocities,
            orientations_r,
            orientations_f,
            orientations_u,
            masses,
            optical_radii,
            // Values
            top_speeds,
            wavelengths,
            is_black_hole,
            is_sphere,
            deprecated,
            positions[index_time_obj(
                end_tick,
                0,
                NUM_OBJECTS)],
            orientations_f[index_time_obj(
                end_tick,
                0,
                NUM_OBJECTS)] +
            orientations_r[index_time_obj(
                end_tick,
                0,
                NUM_OBJECTS)] * (i - height/2) * PIX_SIZE +
            orientations_u[index_time_obj(
                end_tick,
                0,
                NUM_OBJECTS)] * (j - width/2) * PIX_SIZE,
            // Constants for this frame
            cur_time,
            NUM_OBJECTS,
            history_length,
            start_tick,
            end_tick,
            // Constant Settings
            APPLY_GR_RS,
            APPLY_INTENSITY,
            APPLY_LENSING,
            APPLY_SR_RS,
            USE_ALTERNATE_FACTOR,
            USE_INSTANT_GRAVITY,
            USE_LONG_STEP);
    }
}

float schwarzschild_radius(float mass);

int worldline (
        int min_time,
        int max_time,
        int time,
        int object_index,
        __global float3* positions,
        float3 ray_position,
        const int end_tick,
        const int history_length,
        const int num_objects);

float3 acceleration (
        float3 relative_position,
        float3 velocity,
        float mass);

float factor(float ratio, int USE_ALTERNATE_FACTOR);

uchar3 perceived_colour(
    int object_index,
    float3 velocity_hit,
    int time_hit,
    float wavelength,
    float3 initial_velocity,
    float time,
    int num_objects,
    int end_tick,
    int history_length,
    int APPLY_SR_RS,
    int APPLY_GR_RS,
    int APPLY_INTENSITY,
    __global float* masses,
    __global float3* positions,
    __global float3* velocities);

uchar3 ray_trace (
        // Histories
        __global float3* positions,
        __global float3* velocities,
        __global float3* orientations_r,
        __global float3* orientations_f,
        __global float3* orientations_u,
        __global float* masses,
        __global float* optical_radii,
        // Values
        __global float* top_speeds,
        __global float* wavelengths,
        __global int* is_black_hole,
        __global int* is_sphere,
        __global int* deprecated,
        float3 start_point,
        float3 start_ray,
        // Constants for this frame
        const float cur_time,
        const int num_objects,
        const int history_length,
        const int start_tick,
        const int end_tick,
        // Constants for all frames
        const int APPLY_GR_RS,
        const int APPLY_INTENSITY,
        const int APPLY_LENSING,
        const int APPLY_SR_RS,
        const int USE_ALTERNATE_FACTOR,
        const int USE_INSTANT_GRAVITY,
        const int USE_LONG_STEP) {
    float gravitational_radii[NUM_OBJECTS];
    float relevant_radii[NUM_OBJECTS];
    int sizeable_objects[NUM_OBJECTS];
    float distance_traveled = 0.0;
    int num_sizeable = num_objects;
    int ray_time = 0;
    float3 ray_position = start_point;
    float combined_ratio = 0.0;

    for (int i = 1; i < num_objects; i++) {
        gravitational_radii[i] = GRAVITATIONAL_RATIO * schwarzschild_radius(masses[i]);
        relevant_radii[i] = max(gravitational_radii[i], optical_radii[i]);
        sizeable_objects[i] = 1;
        combined_ratio += (
            schwarzschild_radius(masses[index_time_obj(end_tick, i, num_objects)]) /
            distance(start_point, positions[index_time_obj(
                end_tick,
                i,
                num_objects)]));
    }

    float3 ray_velocity = (
        SPEED_OF_LIGHT *
        factor(combined_ratio, USE_ALTERNATE_FACTOR) *
        normalize(start_ray));
    while(num_sizeable > 0) {
        for (int i = 0; i < num_objects; i++) {
            if (sizeable_objects[i] == 1 && (deprecated[i/*ndex_time_obj(
                    time_index(
                        ray_time,
                        history_length,
                        end_tick),
                    i,
                    num_objects)*/] || distance_traveled
                >= VISIBLE_PARALLAX * relevant_radii[i])) {
                sizeable_objects[i] = 0;
                num_sizeable--;
            }
        }
        int relevant_objects[NUM_OBJECTS];
        int num_relevant = 0;
        float3 relative_positions[NUM_OBJECTS];
        for (int i = 0; i < num_objects; i++) {
            relative_positions[i] = ray_position - positions[index_time_obj(
                tick_index(
                    ray_time,
                    history_length,
                    end_tick),
                i,
                num_objects)];
        }
        float distance_squared[NUM_OBJECTS];
        for (int i = 0; i < num_objects; i++) {
            distance_squared[i] = dot(relative_positions[i],
                relative_positions[i]);
        }
        for (int i = 0; i < num_objects; i++) {
            if (sizeable_objects[i] == 1 && masses[index_time_obj(
                    tick_index(
                        ray_time,
                        history_length,
                        end_tick),
                    i,
                    num_objects)] >
                MASSIVE_BOUND && length(relative_positions[i]) <
                gravitational_radii[index_time_obj(
                    tick_index(
                        ray_time,
                        history_length,
                        end_tick),
                    i,
                    num_objects)]) {
                relevant_objects[i] = 1;
                num_relevant++;
            }
        }
        /*if (!APPLY_LENSING || num_relevant == 0) {
            float min_step = history_length;
            for (int i = 0; i < num_objects; i++) {
                if (sizeable_objects[i] == 1) {
                    float quadratic_term = (
                        SPEED_OF_LIGHT * SPEED_OF_LIGHT -
                        top_speeds[i] * top_speeds[i]);
                    float linear_term = (
                        SPEED_OF_LIGHT * dot(relative_positions[i],
                        normalize(ray_velocity)) + top_speeds[i] *
                        gravitational_radii[index_time_obj(
                            tick_index(ray_time, history_length, end_tick),
                            i,
                            num_objects)]);
                    min_step = round(min(USE_LONG_STEP ? (linear_term -
                        sqrt(linear_term * linear_term - quadratic_term *
                        (distance_squared[i] -
                        gravitational_radii[index_time_obj(
                            tick_index(ray_time, history_length, end_tick),
                            i,
                            num_objects)] *
                        gravitational_radii[index_time_obj(
                            tick_index(ray_time, history_length, end_tick),
                            i,
                            num_objects)]))) / quadratic_term : 0.5 *
                        quadratic_term / linear_term, min_step));
                }
                ray_position += min_step * ray_velocity;
                ray_time += min_step;
                continue;
            }
        }*/

        float3 gravitational_positions[NUM_OBJECTS];
        for (int i = 0; i < num_objects; i++) {
            gravitational_positions[i] = (
                ray_position - positions[(USE_INSTANT_GRAVITY ?
                    index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects) :
                    worldline(
                        ((int) floor(ray_time + (length(relative_positions[i]) /
                        (SPEED_OF_LIGHT + top_speeds[i])))),
                        ((int) ceil(ray_time + (length(relative_positions[i]) /
                        (SPEED_OF_LIGHT - top_speeds[i])))),
                        ray_time, i, positions, ray_position,
                        end_tick, history_length, num_objects))]);
        }

        float3 relative_velocities[NUM_OBJECTS];
        float3 acceleration1 = 0;
        float3 acceleration2 = 0;
        float3 acceleration3 = 0;
        float3 acceleration4 = 0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration1 += acceleration(gravitational_positions[i],
                    relative_velocities[i],
                    masses[index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects)]);
            }
        }
        float3 velocity2 = ray_velocity + 0.5 * acceleration1;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration2 += acceleration(gravitational_positions[i],
                    length(velocity2) * normalize(velocity2 -
                    velocities[index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects)]),
                    masses[index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects)]);
            }
        }
        float3 velocity3 = ray_velocity + 0.5 * acceleration2;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration3 += acceleration(gravitational_positions[i],
                    length(velocity3) * normalize(velocity3 -
                    velocities[index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects)]),
                    masses[index_time_obj(
                        tick_index(ray_time, history_length, end_tick),
                        i,
                        num_objects)]);
            }
        }
        float3 velocity4 = ray_velocity + acceleration3;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration4 += acceleration(
                    gravitational_positions[i],
                    length(velocity4) * normalize(velocity4 -
                    velocities[tick_index(ray_time, history_length, end_tick)]),
                    masses[tick_index(ray_time, history_length, end_tick)]);
            }
        }

        float3 new_position = ray_position + 0.166667 * (ray_velocity + 2.0 *
            (velocity2 + velocity3) + velocity4);
        ray_time++;
        if (ray_time == 120) {
            return(uchar3)(100*length(ray_velocity),0,0);
        }
        if (ray_time > history_length) {
            return (uchar3) (0, 0, 0);
        }
        combined_ratio = 0.0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                combined_ratio += (
                    schwarzschild_radius(
                        masses[index_time_obj(
                            tick_index(ray_time, history_length, end_tick),
                            i,
                            num_objects)])
                    / distance(
                        new_position,
                        positions[end_tick * num_objects + i]));
            }
        }

        // ray_velocity = new_position - ray_position;
        ray_velocity = (
            SPEED_OF_LIGHT * factor(
                combined_ratio,
                USE_ALTERNATE_FACTOR) *
            normalize(ray_velocity + 0.166667 * (
                    acceleration1 +
                    2 * (acceleration2 + acceleration3) +
                    acceleration4)));
        for (int i = 0; i < num_objects; i++) {
            if (sizeable_objects[i]) {
                if (is_sphere[i]) {
                    if (distance(
                            new_position,
                            positions[index_time_obj(
                                tick_index(
                                    ray_time,
                                    history_length,
                                    end_tick),
                                i,
                                num_objects)]) < (
                            is_black_hole[i] ?
                            schwarzschild_radius(
                                masses[index_time_obj(
                                    tick_index(
                                        ray_time,
                                        history_length,
                                        end_tick),
                                    i,
                                    num_objects)]) :
                            optical_radii[i])) {
                        return perceived_colour(
                            i,
                            ray_velocity,
                            ray_time,
                            wavelengths[index_time_obj(
                                tick_index(
                                    ray_time,
                                    history_length,
                                    end_tick),
                                i,
                                num_objects)],
                            start_ray,
                            cur_time,
                            num_objects,
                            end_tick,
                            history_length,
                            APPLY_SR_RS,
                            APPLY_GR_RS,
                            APPLY_INTENSITY,
                            masses,
                            positions,
                            velocities);
                    }
                }/* else {
                    for (int j = 0; j < num_faces[i]; j++) {
                        float stretch = dot(normals[i][j] - ray_position,
                            normals[i][j]) / dot(new_position - ray_position,
                            normals[i][j]);
                        if (stretch > 0 && stretch < 1) {
                            float3 u = faces[i][j][1] - faces[i][j][0];
                            float3 v = faces[i][j][2] - faces[i][j][0];
                            float3 w = (ray_position + stretch * (new_position -
                                ray_position)) - faces[i][j][0];
                            if (length(cross(u, w)) + length(cross(w, v)) <
                                length(cross(u, v)) && signbit(dot(cross(u, w),
                                cross(u, v))) == 0 && signbit(dot(cross(w, v),
                                cross(u, v))) == 0) {
                                return colours[index_time_obj(
                                    tick_index(ray_time, history_length, end_tick),
                                    i,
                                    num_objects)][j];
                            }
                        }
                    }
                }*/
            }
        }
        distance_traveled += distance(new_position, ray_position);
        ray_position = new_position;
    }
    return (uchar3) (0, 0, 0);
}

int tick_index(int ray_time, int history_length, int end_tick) {
    return (end_tick - ray_time + history_length) % history_length;
}

int index_time_obj (int tick, int object, int num_objects) {
    return tick*num_objects + object;
}

float schwarzschild_radius(float mass) {
    return 2 * GRAVITATIONAL_CONSTANT * mass / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
}

float factor(float ratio, int USE_ALTERNATE_FACTOR) {
    if (USE_ALTERNATE_FACTOR) {
        return (1 - LIGHT_SLOWING_RATIO * ratio) * (1 - LIGHT_SLOWING_RATIO *
        ratio);
    } else {
        return (1 - ratio) / ((1 + ratio) * (1 + ratio) * (1 + ratio));
    }
}

float3 acceleration (
        float3 relative_position,
        float3 velocity,
        float mass) {
    return -1.5 * relative_position * schwarzschild_radius(mass) *
        dot(
        cross(relative_position, velocity),
        cross(relative_position, velocity))
        / length(relative_position)
        / length(relative_position)
        / length(relative_position)
        / length(relative_position)
        / length(relative_position);

}

// THINGS TO ADD:
// use ITO
// convert what is needed to int
int worldline (
        int min_time,
        int max_time,
        int time,
        int object_index,
        __global float3* positions,
        float3 ray_position,
        const int end_tick,
        const int history_length,
        const int num_objects
        ) {
    int lower_bound = max(min_time, 0);
    int upper_bound = min(max_time, history_length - 1);
    int iterations = history_length / 2;
    for (int i = 0; i < iterations; i++) {
        if (upper_bound == lower_bound) {
             return lower_bound;
        }
        int half_time = (lower_bound + upper_bound) / 2;
        if (
            (SPEED_OF_LIGHT * (half_time - time)) >
            distance(
                positions[index_time_obj(
                    tick_index(
                        half_time,
                        history_length,
                        end_tick),
                    object_index,
                    num_objects)],
                ray_position)) {
            // If (ct)>d, we're looking too far in the past.
            upper_bound = half_time;
        } else {
            lower_bound = half_time+1;
        }
    }
    return lower_bound;
}


uchar3 perceived_colour(
        int object_index,
        float3 velocity_hit,
        int time_hit,
        float wavelength,
        float3 initial_velocity,
        float time,
        int num_objects,
        int end_tick,
        int history_length,
        int APPLY_SR_RS,
        int APPLY_GR_RS,
        int APPLY_INTENSITY,
        __global float* masses,
        __global float3* positions,
        __global float3* velocities) {
    if (APPLY_SR_RS) {
        float outgoing_speed = dot(
            normalize(initial_velocity),
            normalize(velocities[index_time_obj(
                tick_index(
                    time,
                    history_length,
                    end_tick),
                0,
                num_objects)]));

        float incoming_speed = dot(
            normalize(velocity_hit),
            normalize(velocities[index_time_obj(
                tick_index(
                    time_hit,
                    history_length,
                    end_tick),
                object_index,
                num_objects)]));

        float ratio = (outgoing_speed + incoming_speed) / (1 -
            fabs(outgoing_speed * incoming_speed));
        wavelength *= sqrt((1 + ratio) / (1 - ratio));
    }
    if (APPLY_GR_RS) {
        for (int i = 0; i < num_objects; i++) {
            if (masses[index_time_obj(
                    tick_index(
                        time,
                        history_length,
                        end_tick),
                    i,
                    num_objects)] > MASSIVE_BOUND) {
                wavelength *= sqrt((1 - schwarzschild_radius(
                    masses[index_time_obj(tick_index(
                    time,
                    history_length,
                    end_tick), i, num_objects)]) /
                    length(positions[index_time_obj(tick_index(
                    time_hit,
                    history_length,
                    end_tick), object_index, num_objects)] -
                    positions[index_time_obj(tick_index(
                    time_hit,
                    history_length,
                    end_tick), i, num_objects)])) /
                    (1 - schwarzschild_radius(
                    masses[index_time_obj(tick_index(
                    time,
                    history_length,
                    end_tick), i, num_objects)]) /
                    length(positions[index_time_obj(tick_index(
                    time,
                    history_length,
                    end_tick), 0, num_objects)] -
                    positions[index_time_obj(tick_index(
                    time,
                    history_length,
                    end_tick), i, num_objects)])));
            }
        }
    }
    float3 rgb;
    if (wavelength >= 380 && wavelength < 440) {
        rgb = (float3) (0.0166667 * (440 - wavelength), 0.0, 1.0);
    } else if (wavelength >= 440 && wavelength < 490) {
        rgb = (float3) (0.0, 0.02 * (wavelength - 440), 1.0);
    } else if (wavelength >= 490 && wavelength < 510) {
        rgb = (float3) (0.0, 1.0, 0.05 * (510 - wavelength));
    } else if (wavelength >= 510 && wavelength < 580) {
        rgb = (float3) (0.0142857 * (wavelength - 510), 1.0, 0.0);
    } else if (wavelength >= 580 && wavelength < 645) {
        rgb = (float3) (1.0, 0.0153846 * (645 - wavelength), 0.0);
    } else if (wavelength >= 645 && wavelength <= 780) {
        rgb = (float3) (1.0, 0.0, 0.0);
    } else {
        rgb = (float3) (0.0, 0.0, 0.0);
    }
    float factor;
    if (wavelength >= 380 && wavelength < 420) {
        factor = 0.3 + 0.0175 * (wavelength - 380);
    } else if (wavelength >= 420 && wavelength < 700) {
        factor = 1.0;
    } else if (wavelength >= 700 && wavelength <= 780) {
        factor = 0.3 + 0.00875 * (780 - wavelength);
    } else {
        factor = 0.0;
    }
    if (APPLY_INTENSITY) {
        float ratio = length(velocities[index_time_obj(
            tick_index(
                time,
                history_length,
                end_tick),
            0,
            num_objects)]) / SPEED_OF_LIGHT;
        factor *= sqrt(1 - ratio * ratio) /
            (1 - ratio * dot(normalize(initial_velocity),
            normalize(velocities[index_time_obj(
                tick_index(
                    time,
                    history_length,
                    end_tick),
                0,
                num_objects)])));

    }
    return (uchar3) ((rgb.x == 0.0 ? 0 : min(round(255 * powr(rgb.x * INTENSITY_FACTOR * factor, 0.8)), 255.0)),
        (rgb.y == 0.0 ? 0 : min(round(255 * powr(rgb.y * INTENSITY_FACTOR * factor, 0.8)), 255.0)),
        (rgb.z == 0.0 ? 0 : min(round(255 * powr(rgb.z * INTENSITY_FACTOR * factor, 0.8)), 255.0)));
}
