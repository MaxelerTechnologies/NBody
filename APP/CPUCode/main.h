//*********************************************************************//
// N-Body Simulation
//
// Author:  Maxeler Technologies
//
// Imperial College Summer School, July 2012
//
//*********************************************************************//

#ifndef MAIN_H_
#define MAIN_H_

/** Maximum allowance for the relative error between CPU and DFE */
#define THRESHOLD 0

/**
 * 3-D coordinates
 */
typedef struct {
    float x;
    float y;
    float z;
} coord3d_t;

/** Descriptor of a particle state */
typedef struct {
    coord3d_t p;
    coord3d_t v;
} particle_t;

/** Input to the DFE */
typedef struct {
    float x;
    float y;
    float z;
    float m;
} dfe_in_t;

/** Output from the DFE */
typedef struct {
    float x;
    float y;
    float z;
    float pad;
} dfe_out_t;


/**
 * \brief Return the smallest multiple of b greater than or equal to a
 * \param a
 * \param b
 * \return The smallest multiple of b greater than or equal to a
 */
static inline unsigned int next_multiple(unsigned int a, unsigned int b)
{
    return (a + b - 1) / b * b;
}

/**
 * \brief Initialise a timer
 * \param v The timer (in s)
 */
static inline void timer_init(double *v)
{
    *v = 0;
}

/**
 * \brief Start the timer
 * \param v The timer (in s)
 */
static inline void timer_start(double *v)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    *v -= tv.tv_sec + tv.tv_usec * 1e-6;
}
/**
 * \brief Stop the timer
 * \param v The timer (in s)
 */
static inline void timer_stop(double *v)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    *v += tv.tv_sec + tv.tv_usec * 1e-6;
}

#endif /* MAIN_H_ */
