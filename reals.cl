__kernel void reals(
        __global float* positions,
        __global float* velocities,
        __global float* orientations_r,
        __global float* orientations_f,
        __global float* orientations_u,
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
        const float cur_time,
        float3 start_point,
        float3 start_ray,) {

        return (float3)
}
