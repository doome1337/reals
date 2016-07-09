__kernel void reals(
        __global float3* positions,
        __global float3* velocities,
        __global float3* orientations_r,
        __global float3* orientations_f,
        __global float3* orientations_u,
        __global float* local_times,
        __global float* masses,
        __global int* deprecated,
        const unsigned int width,
        const unsigned int height,
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
        __global int* deprecated,
	__global float* relevant_radii,
        const float cur_time,
	const int num_objects,
        float3 start_point,
        float3 start_ray,) {
	
	int sizeable_objects[num_objects];
	float distance_traveled = 0.0;
	for (int i = 0; i < num_objects; i++) {
		sizeable_objects[i] = 1;
	}

	int num_sizable = num_objects;
	float ray_time = cur_time;
	float3 ray_position = start_point;
	combined_ratio = 0;
	for (int i = 0; i < num_objects; i++) {
		combined_ratio += schwarzschild_radius(masses[cur_time *
			num_objects + i]) /
			distance(start_point,
			positions[cur_time * num_objects + i]);
	}
	float3 ray_velocity = SPEED_OF_LIGHT * combined_ratio *
		normalize(start_ray);
	
	
	while(num_sizeable > 0) {
		
		for (int i = 0; i < num_objects; i++) {
			if (distance_traveled > VISIBLE_PARALLAX *
				relevant_radii[i]) {
				sizeable_objects[i] = 0;
				num_sizeable--;
			}
		}
		
		int relevant_objects[num_objects];
		
		int num_relevant = num_objects;
		
		float3 relative_positions[num_objects];

		for (int i = 0; i < num_objects; i++) {
			relative_positions[i] = ray_position -
				positions[ray_time * num_objects + i];
		}
		float distance_squared[num_objects];
		for (int i = 0; i < num_objects; i++) {
			distance_squared[i] = dot(relative_positions[i],
				relative_positions[i]);
		}

		for (int i = 0; i < num_objects; i++) {
			if (sizeable_objects[i] == 1 &&
				masses[ray_time * num_objects + i] >
				MASSIVE_BOUND && distance_squared[i] <
				gravitational_radii[ray_time * num_objects + i] *
				gravitational_radii[ray_time * num_objects +
				i]) {
				relevant_objects[num_objects] = 1;
				num_relevant++;
			}
		}
		if (!APPLY_LENSING || num_relevant == 0) {
			min_step = MAXFLOAT;
			for (int i = 0; i < num_objects; i++) {
				if (sizeable_objects[i] == 1) {
					quadratic_term = SPEED_OF_LIGHT *
						SPEED_OF_LIGHT - top_speeds[i] *
						top_speeds[i];
					linear_term = SPEED_OF_LIGHT *
						dot(relative_positions[i],
						normalize(ray_velocity)) +
						top_speed[i] *
						gravitational_radii[ray_time *
						num_objects + i];
					min_step = min(USE_LONG_STEP ?
						(linear_term - sqrt(linear_term *
						linear_term - quadratic_term *
						(distance squared[i] -
						gravitational_radii[ray_time *
						num_objects + i] *
						gravitational_radii[ray_time *
						num_objects + i]))) /
						quadratic_term :
						0.5 * quadratic_term /
						linear_term, min_step);
				}
				ray_position += min_step * velocity;
				ray_time -= min_step;
				continue;
			}
		}
		int gravitational_positions[num_objects];
		for (int i = 0; i < num_objects; i++) {
			gravitational_positions[i] = ray_position -
				positions[(USE_INSTANT_GRAVITY ? ray_time :
				find_worldline_intersection(time -
				relative_position[i] / (SPEED_OF_LIGHT -
				top_speeds[i]), time - relative_position[i] /
				(SPEED_OF_LIGHT + top_speeds[i]), time, i,
				ray_position)) * num_objects + i];
		}
		
	}

        return (float3)
}
