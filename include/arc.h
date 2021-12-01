#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <signal.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <execinfo.h>
#include <getopt.h>
#include <inttypes.h>
#include <omp.h>


// Information Variables
extern int PRINT;
extern int AVAIL_THREADS;
extern int NUM_CONFIGURATIONS;
// User Control Variables
extern double ARC_ANY_THREADS;
extern double ARC_ANY_SIZE;
extern double ARC_ANY_BW;
extern int ARC_ANY_ECC;
// ECC Method Identifier Values
extern int ARC_PARITY;
extern int ARC_HAMMING;
extern int ARC_SECDED;
extern int ARC_RS;
extern int ARC_DET_SPARSE;
extern int ARC_COR_SPARSE;
extern int ARC_COR_BURST;

// Functions 
void arc_help();
int arc_init(uint32_t max_threads);
int arc_save();
int arc_close();
int arc_encode(uint8_t* data, uint32_t data_size, double memory_constraint, double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint8_t** encoded_data, uint32_t* encoded_data_size);
int arc_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** decoded_data, uint32_t* decoded_data_size);
int arc_check(uint8_t* encoded_data, uint32_t encoded_data_size);
int arc_memory_optimizer(double memory_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* memory_choice_parameter_a, uint32_t* memory_choice_parameter_b);
int arc_throughput_optimizer(double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* throughput_choice_parameter_a, uint32_t* throughput_choice_parameter_b, uint32_t* num_threads);
int arc_joint_optimizer(double memory_constraint, double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* optimizer_parameter_a, uint32_t* optimizer_parameter_b, uint32_t* num_threads);
uint8_t arc_calculate_parity_uint64 (uint64_t byte);
uint8_t arc_calculate_parity_uint8 (uint8_t byte);
int arc_parity_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size);
int arc_parity_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size);
uint8_t arc_calculate_hamming_uint64 (uint64_t byte);
uint8_t arc_calculate_hamming_uint8 (uint8_t byte);
int arc_hamming_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size);
int arc_hamming_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size);
uint8_t arc_calculate_secded_uint64 (uint64_t byte);
uint8_t arc_calculate_secded_uint8 (uint8_t byte);
int arc_secded_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size);
int arc_secded_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size);
int arc_reed_solomon_encode(uint8_t* data, uint32_t data_size, uint32_t data_devices, uint32_t code_devices, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size);
int arc_reed_solomon_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t *data_size);

#ifdef __cplusplus
}
#endif
