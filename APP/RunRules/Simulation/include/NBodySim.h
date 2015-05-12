/**\file */
#ifndef SLIC_DECLARATIONS_NBodySim_H
#define SLIC_DECLARATIONS_NBodySim_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define NBodySim_loopLatency (12)
#define NBodySim_numPipes (1)
#define NBodySim_PCIE_ALIGNMENT (16)
#define NBodySim_maxParticles (65535)
#define NBodySim_numParticlesDivisor (12)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface advanced ----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'advanced'.
 * 
 * \param [in] param_N Interface Parameter "N".: Number of particles in system
 * \param [in] param_P Interface Parameter "P".: Number of particles to compute acceleration for
 * \param [in] inscalar_NBodyKernel_EPS Input scalar parameter "NBodyKernel.EPS".
 * \param [in] instream_fromhost The stream should be of size ((param_N * 4) * 4) bytes.
 * \param [out] outstream_acc The stream should be of size ((param_P * 4) * 4) bytes.
 */
void NBodySim_advanced(
	int32_t param_N,
	int32_t param_P,
	double inscalar_NBodyKernel_EPS,
	const float *instream_fromhost,
	float *outstream_acc);

/**
 * \brief Basic static non-blocking function for the interface 'advanced'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_N Interface Parameter "N".: Number of particles in system
 * \param [in] param_P Interface Parameter "P".: Number of particles to compute acceleration for
 * \param [in] inscalar_NBodyKernel_EPS Input scalar parameter "NBodyKernel.EPS".
 * \param [in] instream_fromhost The stream should be of size ((param_N * 4) * 4) bytes.
 * \param [out] outstream_acc The stream should be of size ((param_P * 4) * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *NBodySim_advanced_nonblock(
	int32_t param_N,
	int32_t param_P,
	double inscalar_NBodyKernel_EPS,
	const float *instream_fromhost,
	float *outstream_acc);

/**
 * \brief Advanced static interface, structure for the engine interface 'advanced'
 * 
 */
typedef struct { 
	int32_t param_N; /**<  [in] Interface Parameter "N".: Number of particles in system */
	int32_t param_P; /**<  [in] Interface Parameter "P".: Number of particles to compute acceleration for */
	double inscalar_NBodyKernel_EPS; /**<  [in] Input scalar parameter "NBodyKernel.EPS". */
	const float *instream_fromhost; /**<  [in] The stream should be of size ((param_N * 4) * 4) bytes. */
	float *outstream_acc; /**<  [out] The stream should be of size ((param_P * 4) * 4) bytes. */
} NBodySim_advanced_actions_t;

/**
 * \brief Advanced static function for the interface 'advanced'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void NBodySim_advanced_run(
	max_engine_t *engine,
	NBodySim_advanced_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'advanced'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_advanced_run_nonblock(
	max_engine_t *engine,
	NBodySim_advanced_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'advanced'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void NBodySim_advanced_run_group(max_group_t *group, NBodySim_advanced_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'advanced'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_advanced_run_group_nonblock(max_group_t *group, NBodySim_advanced_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'advanced'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void NBodySim_advanced_run_array(max_engarray_t *engarray, NBodySim_advanced_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'advanced'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_advanced_run_array_nonblock(max_engarray_t *engarray, NBodySim_advanced_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* NBodySim_advanced_convert(max_file_t *maxfile, NBodySim_advanced_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_N Interface Parameter "N".: Number of particles
 * \param [in] inscalar_NBodyKernel_EPS Input scalar parameter "NBodyKernel.EPS".
 * \param [in] instream_fromhost The stream should be of size ((param_N * 4) * 4) bytes.
 * \param [out] outstream_acc The stream should be of size ((param_N * 4) * 4) bytes.
 */
void NBodySim(
	int32_t param_N,
	double inscalar_NBodyKernel_EPS,
	const float *instream_fromhost,
	float *outstream_acc);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_N Interface Parameter "N".: Number of particles
 * \param [in] inscalar_NBodyKernel_EPS Input scalar parameter "NBodyKernel.EPS".
 * \param [in] instream_fromhost The stream should be of size ((param_N * 4) * 4) bytes.
 * \param [out] outstream_acc The stream should be of size ((param_N * 4) * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *NBodySim_nonblock(
	int32_t param_N,
	double inscalar_NBodyKernel_EPS,
	const float *instream_fromhost,
	float *outstream_acc);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	int32_t param_N; /**<  [in] Interface Parameter "N".: Number of particles */
	double inscalar_NBodyKernel_EPS; /**<  [in] Input scalar parameter "NBodyKernel.EPS". */
	const float *instream_fromhost; /**<  [in] The stream should be of size ((param_N * 4) * 4) bytes. */
	float *outstream_acc; /**<  [out] The stream should be of size ((param_N * 4) * 4) bytes. */
} NBodySim_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void NBodySim_run(
	max_engine_t *engine,
	NBodySim_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_run_nonblock(
	max_engine_t *engine,
	NBodySim_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void NBodySim_run_group(max_group_t *group, NBodySim_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_run_group_nonblock(max_group_t *group, NBodySim_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void NBodySim_run_array(max_engarray_t *engarray, NBodySim_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *NBodySim_run_array_nonblock(max_engarray_t *engarray, NBodySim_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* NBodySim_convert(max_file_t *maxfile, NBodySim_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* NBodySim_init(void);

/* Error handling functions */
int NBodySim_has_errors(void);
const char* NBodySim_get_errors(void);
void NBodySim_clear_errors(void);
/* Free statically allocated maxfile data */
void NBodySim_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int NBodySim_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int NBodySim_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_NBodySim_H */

