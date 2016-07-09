// TODO LIST:
// Make ray_time an int
//
__kernel void reals (
    __global float3* positions,
    __global float3* velocities,
    __global float3* orientations_r,
    __global float3* orientations_f,
    __global float3* orientations_u,
    __global float* local_times,
    __global float* masses,
    __global float* optical_radii,
    __global int* deprecated,
    const unsigned int width,
    const unsigned int height,
    const unsigned int num_objects,
    const unsigned int history_length,
    const unsigned int start_tick,
    const unsigned int end_tick,
    const float cur_time) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if (i < width && j < height) {
        1+1;
    }
}

// float wavelength_redshift(
//        float3 light_direction,
//        ,
//        float impact_time,
//        ) {
//
// }

float schwarzschild_radius(
        float mass,
        const float GRAVITATIONAL_CONSTANT,
        const float SPEED_OF_LIGHT);

int find_worldline_intersection (
        int min_time,
        int max_time,
        int time,
        int object_index,
        float3* positions,
        float3 ray_position,
        const float SPEED_OF_LIGHT,
        const int end_tick,
        const int history_length);

float3 acceleration (
        float3 relative_position,
        float3 velocity);

float factor(
        float ratio,
        const int USE_ALTERNATE_FACTOR,
        const float LIGHT_SLOWING_RATIO);

uchar3 ray_trace(
        __global float3* positions,
        __global float3* velocities,
        __global float3* orientations_r,
        __global float3* orientations_f,
        __global float3* orientations_u,
        __global float* masses,
        __global float* optical_radii,
        __global float* top_speeds,
        __global int* is_black_hole,
        __global int* is_sphere,
        __global int* deprecated,
        float3 start_point,
        float3 start_ray,
        const float cur_time,
        const float GRAVITATIONAL_CONSTANT,
        const float GRAVITATIONAL_RATIO,
        const float LIGHT_SLOWING_RATIO,
        const float MASSIVE_BOUND,
        const float SPEED_OF_LIGHT,
        const float VISIBLE_PARALLAX,
        const int num_objects,
        const int history_length,
        const int start_tick,
        const int end_tick,
        const int APPLY_LENSING,
        const int USE_ALTERNATE_FACTOR,
        const int USE_INSTANT_GRAVITY,
        const int USE_LONG_STEP
        ) {
    float gravitational_radii[num_objects];
    float relevant_radii[num_objects];
    int sizeable_objects[num_objects];
    float distance_traveled = 0.0;
    int num_sizeable = num_objects;
    float ray_time = 0.0;
    float3 ray_position = start_point;
    float combined_ratio = 0.0;

    for (int i = 0; i < num_objects; i++) {
        gravitational_radii[i] = GRAVITATIONAL_RATIO *
            schwarzschild_radius(masses[i], GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT);
        relevant_radii[i] = max(gravitational_radii[i], optical_radii[i]);
        sizeable_objects[i] = 1;
        combined_ratio += (
            schwarzschild_radius(masses[end_tick * num_objects + i], GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT) /
            distance(start_point, positions[end_tick * num_objects + i]));
    }
    float3 ray_velocity = SPEED_OF_LIGHT * combined_ratio * normalize(start_ray);

    while(num_sizeable > 0) {
        for (int i = 0; i < num_objects; i++) {
            if (sizeable_objects[i] == 1 && (deprecated[i] || distance_traveled
                > VISIBLE_PARALLAX * relevant_radii[i])) {
                sizeable_objects[i] = 0;
                num_sizeable--;
            }
        }

        int relevant_objects[num_objects];
        int num_relevant = 0;
        float3 relative_positions[num_objects];
        for (int i = 0; i < num_objects; i++) {
            relative_positions[i] = ray_position - positions[tick_index(ray_time,
                history_length, end_tick)];
        }
        float distance_squared[num_objects];
        for (int i = 0; i < num_objects; i++) {
            distance_squared[i] = dot(relative_positions[i],
                relative_positions[i]);
        }
        for (int i = 0; i < num_objects; i++) {
            if (sizeable_objects[i] == 1 && masses[tick_index(ray_time,
                history_length, end_tick)] >
                MASSIVE_BOUND && distance_squared[i] <
                gravitational_radii[tick_index(ray_time, history_length,
                end_tick)] * gravitational_radii[tick_index(ray_time,
                history_length, end_tick)]) {
                relevant_objects[num_objects] = 1;
                num_relevant++;
            }
        }

        if (!APPLY_LENSING || num_relevant == 0) {
            float min_step = MAXFLOAT;
            for (int i = 0; i < num_objects; i++) {
                if (sizeable_objects[i] == 1) {
                    float quadratic_term = (
                        SPEED_OF_LIGHT * SPEED_OF_LIGHT -
                        top_speeds[i] * top_speeds[i]);
                    float linear_term = (
                        SPEED_OF_LIGHT * dot(relative_positions[i],
                        normalize(ray_velocity)) + top_speeds[i] *
                        gravitational_radii[
                            tick_index(ray_time, history_length, end_tick)]);
                    min_step = min(USE_LONG_STEP ? (linear_term -
                        sqrt(linear_term * linear_term - quadratic_term *
                        (distance_squared[i] -
                        gravitational_radii[tick_index(ray_time, history_length,
                        end_tick)] * gravitational_radii[tick_index(ray_time,
                        history_length, end_tick)]))) / quadratic_term : 0.5 *
                        quadratic_term / linear_term, min_step);
                }
                ray_position += min_step * ray_velocity;
                ray_time += min_step;
                continue;
            }
        }

        float3 gravitational_positions[num_objects];
        for (int i = 0; i < num_objects; i++) {
            gravitational_positions[i] = (
                ray_position - positions[(USE_INSTANT_GRAVITY ?
                    tick_index(ray_time, history_length, end_tick) :
                    find_worldline_intersection(
                        floor(ray_time + (length(relative_positions[i]) /
                        (SPEED_OF_LIGHT + top_speeds[i]))),
                        ceil(ray_time + (length(relative_positions[i]) /
                        (SPEED_OF_LIGHT - top_speeds[i]))),
                        ray_time, i,
                        positions, ray_position,
                        SPEED_OF_LIGHT, end_tick, history_length))]);
        }

        float3 relative_velocities[num_objects];
        float3 acceleration1 = 0;
        float3 acceleration2 = 0;
        float3 acceleration3 = 0;
        float3 acceleration4 = 0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration1 += acceleration(gravitational_positions[i],
                    relative_velocities[i],
                    masses[tick_index(ray_time, history_length, end_tick)]);
            }
        }
        float3 velocity2 = ray_velocity + 0.5 * acceleration1;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration2 += acceleration(gravitational_positions[i],
                    length(velocity2) * normalize(velocity2 -
                    velocities[tick_index(ray_time, history_length, end_tick)],
                    masses[tick_index(ray_time, history_length, end_tick)]));
            }
        }
        float3 velocity3 = ray_velocity + 0.5 * acceleration2;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration3 += acceleration(gravitational_positions[i],
                    length(velocity3) * normalize(velocity3 -
                    velocities[tick_index(ray_time, history_length, end_tick)],
                    masses[tick_index(ray_time, history_length, end_tick)]));
            }
        }
        float3 velocity4 = ray_velocity + acceleration3;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration4 += acceleration(gravitational_positions[i],
                    length(velocity4) * normalize(velocity4 -
                    velocities[tick_index(ray_time, history_length, end_tick)],
                    masses[tick_index(ray_time, history_length, end_tick)]));
            }
        }

        float3 new_position = ray_position + 0.166667 * (ray_velocity + 2.0 *
            (velocity2 + velocity3) + velocity4);
        ray_time++;
        combined_ratio = 0.0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                combined_ratio += (
                    schwarzschild_radius(
                        masses[tick_index(ray_time, history_length, end_tick)],
                        GRAVITATIONAL_CONSTANT,
                        SPEED_OF_LIGHT)
                    / distance(
                        new_position,
                        positions[end_tick * num_objects + i]));
            }
        }

        ray_velocity = (
            SPEED_OF_LIGHT * factor(
                combined_ratio,
                USE_ALTERNATE_FACTOR,
                LIGHT_SLOWING_RATIO) *
            normalize(ray_velocity + 0.166667 * (
                    acceleration1 +
                    2 * (acceleration2 + acceleration3) +
                    acceleration4)));

        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                if (is_sphere[i]) {
                    if (distance(
                            new_position,
                            positions[tick_index(
                                ray_time,
                                history_length,
                                end_tick)]) < (
                            is_black_hole[i] ?
                            schwarzschild_radius(
                                masses[tick_index(
                                    ray_time,
                                    history_length,
                                    end_tick)],
                                GRAVITATIONAL_CONSTANT,
                                SPEED_OF_LIGHT) :
                            optical_radii[i])) {
                        //return colours[i];
                        return (uchar3) (0, 0, 0);
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
                                return colours[tick_index(ray_time,
                                history_length, end_tick)][j];
                            }
                        }
                    }
                }*/
            }
        }
        distance_traveled += distance(new_position, ray_position);
        ray_position = new_position;
        ray_velocity = new_position - ray_position;
    }

    return (uchar3) (0, 0, 0);
}

int tick_index(int ray_time, int history_length, int end_tick) {
    return (end_tick - ray_time + history_length) % history_length;
}

int index_time_obj (int tick, int object, int num_objects) {
    return tick*num_objects + object;
}

float schwarzschild_radius(
        float mass,
        const float GRAVITATIONAL_CONSTANT,
        const float SPEED_OF_LIGHT) {
    return 2 * GRAVITATIONAL_CONSTANT * mass / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
}

float factor(
        float ratio,
        const int USE_ALTERNATE_FACTOR,
        const float LIGHT_SLOWING_RATIO) {
    if (USE_ALTERNATE_FACTOR) {
        return (1 - LIGHT_SLOWING_RATIO * ratio) * (1 - LIGHT_SLOWING_RATIO *
        ratio);
    } else {
        return (1 - ratio) / (1 + ratio) / (1 + ratio) / (1 + ratio);
    }
}

float3 acceleration (float3 relative_position, float3 velocity, float mass) {
    return -1.5 * relative_position * schwarzschild_radius(mass) * dot(
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
int find_worldline_intersection (
        int min_time,
        int max_time,
        int time,
        int object_index,
        float3* positions,
        float3 ray_position,
        const float SPEED_OF_LIGHT,
        const int end_tick,
        const int history_length) {
    int lower_bound = max(min_time, 0);
    int upper_bound = min(max_time, history_length-1);
    const int iterations = 2*log2(history_length);
    for (int i = 0; i < iterations; i++) {
        if (upper_bound == lower_bound) {
             return lower_bound;
        }
        int half_time = (lower_bound + upper_bound) / 2;
        if (
            (
                SPEED_OF_LIGHT * (half_time - time)) >
            distance(
                positions[index_time_obj(
                    tick_index(
                        half_time,
                        history_length,
                        end_tick),
                    object_index)],
                ray_position)) {
            // If (ct)>d, we're looking too far in the past.
            upper_bound = half_time;
        } else {
            lower_bound = half_time+1;
        }
    }
    return lower_bound;
}
