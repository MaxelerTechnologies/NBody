//*********************************************************************//
// N-Body Simulation
//
// Author:  Maxeler Technologies
//
// Imperial College Summer School, July 2012
//
//*********************************************************************//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
//#include <NBody.h>
#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#include "main.h"
#include "cmdline.h"


/**
 * \brief Count the number of lines in a file
 * \param [in] fp       File pointer
 * \return The number of lines
 */
static int count_lines(FILE *fp)
{
    int nl = 0;
    int el = 0;
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        if (strchr(buf, '\n')) {
            nl++;
            el = 0;
        } else {
            el = 1;
        }
    }
    return nl + el;
}

/**
 * \brief Run the N-body simulation on the CPU.
 * \param [in]  N               Number of particles
 * \param [in]  nt              Number of time-steps
 * \param [in]  EPS             Damping factor
 * \param [in]  m               Masses of the N particles
 * \param [in]  in_particles    Initial state of the N particles
 * \param [out] out_particles   Final state of the N particles after nt time-steps
 * \param [out] time            Execution time
 */
static void run_cpu(int N, int nt, float EPS, const float *m,
        const particle_t *in_particles, particle_t *out_particles,
        double *time)
{
    particle_t *p = (particle_t *) malloc(N * sizeof(particle_t));
    memcpy(p, in_particles, N * sizeof(particle_t));

    coord3d_t *a = (coord3d_t *) malloc(N * sizeof(coord3d_t));

    double wall_time;
    timer_init(&wall_time);

    timer_start(&wall_time);

    for (int t = 0; t < nt; t++) {

        memset(a, 0, N * sizeof(coord3d_t));

        for (int q = 0; q < N; q++) {
            for (int j = 0; j < N; j++) {
                //if (j != q) {
                    float rx = p[j].p.x - p[q].p.x;
                    float ry = p[j].p.y - p[q].p.y;
                    float rz = p[j].p.z - p[q].p.z;
                    float dd = rx*rx + ry*ry + rz*rz + EPS;
                    float d = 1/ (dd*sqrtf(dd));
                    float s = m[j] * d;
                    a[q].x += rx * s;
                    a[q].y += ry * s;
                    a[q].z += rz * s;
                //}
            }
        }

        for (int i = 0; i < N; i++) {
            p[i].p.x += p[i].v.x;
            p[i].p.y += p[i].v.y;
            p[i].p.z += p[i].v.z;
            p[i].v.x += a[i].x;
            p[i].v.y += a[i].y;
            p[i].v.z += a[i].z;
        }

    }

    timer_stop(&wall_time);

    *time = wall_time;

    memcpy(out_particles, p, N * sizeof(particle_t));

    free(p);
    free(a);
}

/**
 * \brief Run the N-body simulation on the CPU.
 * \param [in]  N               Number of particles
 * \param [in]  nt              Number of time-steps
 * \param [in]  EPS             Damping factor
 * \param [in]  m               Masses of the N particles
 * \param [in]  in_particles    Initial state of the N particles
 * \param [out] out_particles   Final state of the N particles after nt time-steps
 * \param [out] time            Execution time
 */
static void run_model(int N, int nt, float EPS, const float *m,
        const particle_t *in_particles, particle_t *out_particles,
        double *time)
{
    particle_t *p = (particle_t *) malloc(N * sizeof(particle_t));
    memcpy(p, in_particles, N * sizeof(particle_t));

    coord3d_t *a = (coord3d_t *) malloc(N * sizeof(coord3d_t));

    double wall_time;
    timer_init(&wall_time);

    timer_start(&wall_time);

    const int latency  = 12;

    for (int t = 0; t < nt; t++) {

        memset(a, 0, N * sizeof(coord3d_t));

        coord3d_t *acc = a;
        particle_t *pIn = p;

        for(int iTiled = 0; iTiled < N/latency; ++iTiled, acc += latency, pIn += latency) {
			for(int j = 0; j < N; ++j) {
				for(int latencyCount = 0; latencyCount < latency; ++latencyCount) {
                    float rx = p[j].p.x - pIn[latencyCount].p.x;
                    float ry = p[j].p.y - pIn[latencyCount].p.y;
                    float rz = p[j].p.z - pIn[latencyCount].p.z;
                    float dd = rx*rx + ry*ry + rz*rz + EPS;
                    float d = 1/ (dd*sqrtf(dd));
                    float s = m[j] * d;
                    acc[latencyCount].x += rx * s;
                    acc[latencyCount].y += ry * s;
                    acc[latencyCount].z += rz * s;
				}
			}
        }

        for (int i = 0; i < N; i++) {
            p[i].p.x += p[i].v.x;
            p[i].p.y += p[i].v.y;
            p[i].p.z += p[i].v.z;
            p[i].v.x += a[i].x;
            p[i].v.y += a[i].y;
            p[i].v.z += a[i].z;
        }

    }

    timer_stop(&wall_time);

    *time = wall_time;

    memcpy(out_particles, p, N * sizeof(particle_t));

    free(p);
    free(a);
}

/**
 * \brief Run the N-body simulation on the DFE.
 * \param [in]  N               Number of particles
 * \param [in]  nt              Number of time-steps
 * \param [in]  EPS             Damping factor
 * \param [in]  m               Masses of the N particles
 * \param [in]  in_particles    Initial state of the N particles
 * \param [out] out_particles   Final state of the N particles after nt time-steps
 * \param [out] time            Execution time
 */
static void run_dfe(int N, int nt, float EPS, const float *m,
        const particle_t *in_particles, particle_t *out_particles,
        double *time)
{
    /* Initialise the DFE */
    max_file_t *maxfile = NBody_init();

    const int latency = max_get_constant_uint64t(maxfile, "loopLatency");
    const int max_particles = max_get_constant_uint64t(maxfile, "maxParticles");
    const int num_particles_divisor = max_get_constant_uint64t(maxfile, "numParticlesDivisor");

    int newN;
    if (N < latency) {
        newN = next_multiple(latency, num_particles_divisor);
        fprintf(stderr, "info: The number of particles on the DFE has been padded to %d\n", newN);
    } else if (N > max_particles) {
        fprintf(stderr, "error: the number of particles must be less than or equal to %d\n", max_particles);
        exit(EXIT_FAILURE);
    } else if (N % num_particles_divisor != 0) {
        newN = next_multiple(N, num_particles_divisor);
        fprintf(stderr, "info: The number of particles on the DFE has been padded to %d\n", newN);
    } else {
        newN = N;
    }

    /* Copy initial state of particles */
    dfe_in_t *dfe_particles = (dfe_in_t *) malloc(newN * sizeof(dfe_in_t));
    for (int i = 0; i < N; i++) {
        dfe_particles[i].x = in_particles[i].p.x;
        dfe_particles[i].y = in_particles[i].p.y;
        dfe_particles[i].z = in_particles[i].p.z;
        dfe_particles[i].m = m[i];
    }
    for (int i = N; i < newN; i++) {
        /* Setting a null mass on the padded particles ensures they have no effect */
        dfe_particles[i].m = 0;
    }

    /* Set up output streams */
    dfe_out_t *acc = (dfe_out_t *) malloc(newN * sizeof(dfe_out_t));

    particle_t *cur_particles = (particle_t *) malloc(N * sizeof(particle_t));
    memcpy(cur_particles, in_particles, N * sizeof(particle_t));

    double wall_time, run_time, update_time;
    timer_init(&wall_time);
    timer_init(&run_time);
    timer_init(&update_time);

    timer_start(&wall_time);

    for (int t = 0; t < nt; t++) {

        timer_start(&run_time);
        NBody(newN, EPS, (float *) dfe_particles, (float *) acc);
        timer_stop(&run_time);

        timer_start(&update_time);

        /* Update the state of particles */
        for (int i = 0; i < N; i++) {
            cur_particles[i].p.x += cur_particles[i].v.x;
            cur_particles[i].p.y += cur_particles[i].v.y;
            cur_particles[i].p.z += cur_particles[i].v.z;
            cur_particles[i].v.x += acc[i].x;
            cur_particles[i].v.y += acc[i].y;
            cur_particles[i].v.z += acc[i].z;
        }

        for (int i = 0; i < N; i++) {
            dfe_particles[i].x = cur_particles[i].p.x;
            dfe_particles[i].y = cur_particles[i].p.y;
            dfe_particles[i].z = cur_particles[i].p.z;
        }
        timer_stop(&update_time);
    }

    timer_stop(&wall_time);

    printf("Wall clock time:   %-12gs\n", wall_time);
    printf("Run time:          %-12gs (%.1f%%)\n", run_time, run_time/wall_time*100);
    printf("Update time:       %-12gs (%.1f%%)\n", update_time, update_time/wall_time*100);

    *time = wall_time;

    max_file_free(maxfile);

    memcpy(out_particles, cur_particles, N * sizeof(particle_t));

    free(acc);
    free(cur_particles);
    free(dfe_particles);

}

/**
 * \brief Check if the relative error between two numbers is
 *        less than a given threshold.
 * \param [in]  val         Inferred value
 * \param [in]  ref         Reference value
 * \param [in]  threshold   Relative error threshold
 */
static int check_error(float val, float ref, float threshold)
{
    int error;
    if (ref == 0) {
        error = val != 0;
    } else if (isnan(val) || isnan(ref)) {
        error = 1;
    } else {
        error = fabs((val - ref) / ref) > threshold;
    }
    return error;
}

int main(int argc, char **argv)
{
    struct gengetopt_args_info args_info;

    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(EXIT_FAILURE);
    }

    int N = args_info.num_particles_arg;
    int nt = args_info.num_timesteps_arg;
    float EPS = args_info.EPS_arg;

    if (EPS == 0) {
        fprintf(stderr, "EPS cannot be set to zero\n");
        exit(EXIT_FAILURE);
    }

    particle_t *particles;
    float *m;

    if (!args_info.random_given && !args_info.file_given) {
        cmdline_parser_print_help();
        exit(EXIT_FAILURE);
    }

    if (args_info.random_given) {
        particles = (particle_t *) malloc(N * sizeof(particle_t));
        m = (float *) malloc(N * sizeof(float));

        srand(100);
        for (int i = 0; i < N; i++)
        {
            m[i] = (float)rand()/100000;
            particles[i].p.x = (float)rand()/100000;
            particles[i].p.y = (float)rand()/100000;
            particles[i].p.z = (float)rand()/100000;
            particles[i].v.x = (float)rand()/100000;
            particles[i].v.y = (float)rand()/100000;
            particles[i].v.z = (float)rand()/100000;
        }
    } else {
        const char *filename = args_info.file_arg;

        FILE *fp = fopen(args_info.file_arg, "r");
        if (fp == NULL) {
            fprintf(stderr, "Failed to open input file: `%s'\n", filename);
            exit(EXIT_FAILURE);
        }

        N = count_lines(fp) - 1;

        if (args_info.num_particles_given &&
            args_info.num_particles_arg < N) {
                N = args_info.num_particles_arg;
        }

        particles = (particle_t *) malloc(N * sizeof(particle_t));
        m = (float *) malloc(N * sizeof(float));

        rewind(fp);

        fscanf(fp, "m,x,y,z,vx,vy,vz\n");
        for (int i = 0; i < N; i++) {
            fscanf(fp, "%g,%g,%g,%g,%g,%g,%g", &m[i],
                    &particles[i].p.x, &particles[i].p.y, &particles[i].p.z,
                    &particles[i].v.x, &particles[i].v.y, &particles[i].v.z);
        }

        fclose(fp);
    }

    double dfeTime = 0;
    particle_t *dfe_particles = NULL;

    double cpuTime = 0;
    particle_t *cpu_particles = NULL;

    if (!args_info.cpu_given) {
        /* Not CPU-only: Run the DFE */
        dfe_particles = (particle_t *) malloc(N * sizeof(particle_t));

        if(!args_info.model_given) {
			puts("Running on DFE...\n");
			run_dfe(N, nt, EPS, m, particles, dfe_particles, &dfeTime);
        } else {
			puts("Running on DFE Model...\n");
        	run_model(N, nt, EPS, m, particles, dfe_particles, &dfeTime);
        }
        printf("DFE execution time: %.3gs\n", dfeTime);
    }

    if (!args_info.dfe_given) {
        /* Not DFE-only: Run the CPU */
        cpu_particles = (particle_t *) malloc(N * sizeof(particle_t));

        puts("Running on CPU...\n");
        run_cpu(N, nt, EPS, m, particles, cpu_particles, &cpuTime);
        printf("CPU execution time: %.3gs\n", cpuTime);
    }

    if (!args_info.dfe_given && !args_info.cpu_given) {

        printf("Speed-up (1 card vs. 1 thread): %.1fx\n", cpuTime / dfeTime);
        printf("Speed-up (node to node):        %.1fx\n", (cpuTime/12) / (dfeTime/4));

        printf("Checking results...\n");

        int error = 0;
        for (int i = 0; i < N; i++) {
            if (check_error(dfe_particles[i].p.x, cpu_particles[i].p.x, THRESHOLD) ||
                    check_error(dfe_particles[i].p.y, cpu_particles[i].p.y, THRESHOLD) ||
                    check_error(dfe_particles[i].p.z, cpu_particles[i].p.z, THRESHOLD) ||
                    check_error(dfe_particles[i].v.x, cpu_particles[i].v.x, THRESHOLD) ||
                    check_error(dfe_particles[i].v.y, cpu_particles[i].v.y, THRESHOLD) ||
                    check_error(dfe_particles[i].v.z, cpu_particles[i].v.z, THRESHOLD)) {
                error = 1;
                printf("Difference on particle %d\n", i);
                printf("              px                py              pz              vx              vy              vz\n");
                printf("CPU:  %12g\t%12g\t%12g\t%12g\t%12g\t%12g\n",
                        cpu_particles[i].p.x, cpu_particles[i].p.y, cpu_particles[i].p.z,
                        cpu_particles[i].v.x, cpu_particles[i].v.y, cpu_particles[i].v.z);
                printf("DFE:  %12g\t%12g\t%12g\t%12g\t%12g\t%12g\n",
                        dfe_particles[i].p.x, dfe_particles[i].p.y, dfe_particles[i].p.z,
                        dfe_particles[i].v.x, dfe_particles[i].v.y, dfe_particles[i].v.z);
                printf("Diff: %12g\t%12g\t%12g\t%12g\t%12g\t%12g\n",
                        fabs(dfe_particles[i].p.x-cpu_particles[i].p.x),
                        fabs(dfe_particles[i].p.y-cpu_particles[i].p.y),
                        fabs(dfe_particles[i].p.z-cpu_particles[i].p.z),
                        fabs(dfe_particles[i].v.x-cpu_particles[i].v.x),
                        fabs(dfe_particles[i].v.y-cpu_particles[i].v.y),
                        fabs(dfe_particles[i].v.z-cpu_particles[i].v.z));
                break;
            }
        }
        if (error) {
            printf("FAILED\n");
        } else {
            printf("PASSED\n");
        }
    }

    cmdline_parser_free (&args_info);

    free(particles);
    free(m);
    free(dfe_particles);
    free(cpu_particles);

    return 0;
}
