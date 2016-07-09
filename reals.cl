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
    const float cur_time,
    ) {
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

float3 ray_trace(
    __global float* positions,
    __global float* velocities,
    __global float* orientations_r,
    __global float* orientations_f,
    __global float* orientations_u,
    __global float* masses,
    __global float* optical_radii,
    __global int* is_black_hole,
    __global int* deprecated,
    const float cur_time,
    const int num_objects,
    float3 start_point,
    float3 start_ray) {
    float gravitational_radii[num_objects];
    float relevant_radii[num_objects];
    int sizeable_objects[num_objects];
    float distance_traveled = 0.0;
    int num_sizable = num_objects;
    float ray_time = 0;
    float min_step;
    float3 ray_position = start_point;
    combined_ratio = 0;
    for (int i = 0; i < num_objects; i++) {
	gravitational_radii[i] = GRAVITATIONAL_RATIO *
            schwarzschild_radius(masses[i]);
	relevant_radii[i] = max(gravitational_radii[i], optical_radii[i]);
        sizeable_objects[i] = 1;
        combined_ratio += schwarzschild_radius(masses[tick_index(cur_time,
            history_length, end_tick)]) / distance(start_point,
            positions[tick_index(cur_time, history_length, end_tick)]);
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
            min_step = MAXFLOAT;
            for (int i = 0; i < num_objects; i++) {
                if (sizeable_objects[i] == 1) {
                    quadratic_term = SPEED_OF_LIGHT * SPEED_OF_LIGHT -
                    top_speeds[i] * top_speeds[i];
                    linear_term = SPEED_OF_LIGHT * dot(relative_positions[i],
                        normalize(ray_velocity)) + top_speed[i] *
                        gravitational_radii[tick_index(ray_time, history_length,
                        end_tick)];
                    min_step = min(USE_LONG_STEP ? (linear_term -
                        sqrt(linear_term * linear_term - quadratic_term *
                        (distance squared[i] -
                        gravitational_radii[tick_index(ray_time, history_length,
                        end_tick)] * gravitational_radii[tick_index(ray_time,
			history_length, end_tick)]))) / quadratic_term : 0.5 *
			quadratic_term / linear_term, min_step);
                }
                ray_position += min_step * velocity;
                ray_time += min_step;
                continue;
            }
        }
        int gravitational_positions[num_objects];
        for (int i = 0; i < num_objects; i++) {
            gravitational_positions[i] = ray_position -
                positions[(USE_INSTANT_GRAVITY ? tick_index(ray_time,
		history_length, end_tick) : find_worldline_intersection(time -
		relative_positions[i] / (SPEED_OF_LIGHT - top_speeds[i]), time -
		relative_positions[i] / (SPEED_OF_LIGHT + top_speeds[i]), time, i,
		ray_position))];
        }
        int relative_velocities[num_objects]
        float3 acceleration1 = 0;
        float3 acceleration2 = 0;
        float3 acceleration3 = 0;
        float3 acceleration4 = 0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration1 += acceleration(gravitational_positions[i],
                    relative_velocities[i]);
            }
        }
        float3 velocity2 = ray_velocity + 0.5 * acceleration1;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration2 += acceleration(gravitational_positions[i],
                    length(velocity2) * normalize(velocity2 -
                    velocities[tick_index(ray_time, history_length, end_tick)]));
            }
        }
        float3 velocity3 = ray_velocity + 0.5 * acceleration2;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration3 += acceleration(gravitational_positions[i],
                    length(velocity3) * normalize(velocity3 -
                    velocities[tick_index(ray_time, history_length, end_tick)]));
            }
        }
        float3 velocity4 = ray_velocity + acceleration3;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                acceleration4 += acceleration(gravitational_positions[i],
                    length(velocity4) * normalize(velocity4 -
                    velocities[tick_index(ray_time, history_length, end_tick)]));
            }
        }
        float3 new_position = ray_position + 0.166667 * (ray_velocity + 2 *
            (velocity2 + velocity3) + velocity4);
        ray_time++;
        combined_ratio = 0;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                combined_ratio +=
                    schwarzschild_radius(masses[tick_index(ray_time,
                    history_length, end_tick)]) / distance(new_position,
                    positions[cur_time * num_objects + i]);
            }
        }
        ray_velocity = SPEED_OF_LIGHT * factor(combined_ratio) *
        normalize(ray_velocity + 0.166667 * (acceleration1 + 2 * (acceleration2 +
            acceleration3) + acceleration4));
	distance_travelled += distance(new_position, ray_position);
        //float3 new_position = ray_position + ray_velocity;
	//ray_time++;
        for (int i = 0; i < num_objects; i++) {
            if (relevant_objects[i] == 1) {
                if (is_sphere[i]) {
                    if (distance(new_position, positions[tick_index(ray_time,
                        history_length, end_tick)]) < is_black_hole[i] == 1 ?
                        schwarzschild_radius(masses[tick_index(ray_time,
                        history_length, end_tick)]) : optical_radii[i]) {
                        return colours[i];
                    }
                } else {
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
                }
            }
        }
        ray_position = new_position;
        ray_velocity = new_position - ray_position;
    }

    return (float3)
}

int tick_index(int ray_time, int history_length, int end_tick) {
    return (end_tick - ray_time + history_length) % history_length;
}

float schwarzschild_radius(float mass) {
    return 2 * GRAVITATIONAL_CONSTANT * mass / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
}

float factor(float ratio) {
    if (USE_ALTERNATE_FACTOR) {
        return (1 - LIGHT_SLOWING_RATIO * ratio) * (1 - LIGHT_SLOWING_RATIO *
        ratio);
    } else {
        return (1 - ratio) / (1 + ratio) / (1 + ratio) / (1 + ratio);
    }
}
