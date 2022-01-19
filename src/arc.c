// Libraries
#include <sys/stat.h>
#include <errno.h>
#include "arc.h"
#include "jerasure.h"
#include "reed_sol.h"
// Type Malloc Macro
#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

// GLOBAL Variables Section
// ############################
// Determines whether initialization has occurred
int INIT = 0;
// Determines whether to print extra information or not
int PRINT = 1;
// Maximum number of threads ARC should use
int AVAIL_THREADS;
// Total number of configuratiosn to choose from
int NUM_CONFIGURATIONS;
// Identifier for free thread constraint
double ARC_ANY_THREADS = -1;
// Identifier for free memory constraint
double ARC_ANY_SIZE = -1;
// Identifier for free throughput constraint
double ARC_ANY_BW = -1;
// Identifier for free resiliency constraint
int ARC_ANY_ECC = 0;
// Identifier to use Parity ECC
int ARC_PARITY = -1;
// Identifier to use Hamming ECC
int ARC_HAMMING = -2;
// Identifier to use SECDED ECC
int ARC_SECDED = -3;
// Identifier to use Reed-Solomon ECC
int ARC_RS = -4;
// Identifier to use ECC capable of detecting sparse errors
int ARC_DET_SPARSE = -5;
// Identifier to use ECC capable of correcting sparse errors
int ARC_COR_SPARSE = -6;
// Identifier to use ECC capable of correcting burst errors
int ARC_COR_BURST = -7;

// ECC Method Identifier Variables Section
// ########################################
int NUM_METHODS = 4;
int PARITY_ID = 1;
int HAMMING_ID = 2;
int SECDED_ID = 3;
int RS_ID = 4;

// Resource Variables Section
// ###########################
// Resource Folder Location

// path to cache directories
char *cache_resource_location = NULL;
// Set configuration information cache string
char *thread_resource_file = "_information_cache.csv";
// Hamming & SECDED Resource Variables 
const uint8_t H_S_1_Parity_Matrix[4] = {
  0x5B,
  0x6D,
  0x8E,
  0xF0
};
const uint64_t H_S_8_Parity_Matrix[7] = {
  0xAB55555556AAAD5B,
  0xCD9999999B33366D,
  0xF1E1E1E1E3C3C78E,
  0x01FE01FE03FC07F0,
  0x01FFFE0003FFF800,
  0x01FFFFFFFC000000,
  0xFE00000000000000
};
const uint8_t H_1_Syndrome_Table[16] = {
  0x03, 0x05, 0x06, 0x07, 0x09, 0x0A, 0x0B, 0x0C,
  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80
};
const uint8_t H_8_Syndrome_Table[72] = {
  0x03, 0x05, 0x06, 0x07, 0x09, 0x0A, 0x0B, 0x0C,
  0x0D, 0x0E, 0x0F, 0x11, 0x12, 0x13, 0x14, 0x15,
  0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D,
  0x1E, 0x1F, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26,
  0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E,
  0x2F, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36,
  0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E,
  0x3F, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80
};
const uint8_t S_1_Syndrome_Table[16] = {
  0x13, 0x15, 0x16, 0x07, 0x19, 0x1A, 0x0B, 0x1C,
  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80
};
const uint8_t S_8_Syndrome_Table[72] = {
  0x83, 0x85, 0x86, 0x07, 0x89, 0x8A, 0x0B, 0x8C,
  0x0D, 0x0E, 0x8F, 0x91, 0x92, 0x13, 0x94, 0x15,
  0x16, 0x97, 0x98, 0x19, 0x1A, 0x9B, 0x1C, 0x9D,
  0x9E, 0x1F, 0xA1, 0xA2, 0x23, 0xA4, 0x25, 0x26,
  0xA7, 0xA8, 0x29, 0x2A, 0xAB, 0x2C, 0xAD, 0xAE,
  0x2F, 0xB0, 0x31, 0x32, 0xB3, 0x34, 0xB5, 0xB6,
  0x37, 0x38, 0xB9, 0xBA, 0x3B, 0xBC, 0x3D, 0x3E,
  0xBF, 0xC1, 0xC2, 0x43, 0xC4, 0x45, 0x46, 0xC7,
  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80
};

// ARC Decision Variables Section
// ###############################
// ARC Configuration Information Struct
struct configuration_information {
    int ecc_algorithm;
    uint32_t ecc_parameter_a;
    uint32_t ecc_parameter_b;
    uint32_t num_threads;
    double memory_overhead;
    double throughput_overhead;
};
// ARC Configuration Information Decision Array
struct configuration_information *arc_configurations;

// Utility Functions Section
// ############################
// print_bits:
// Print binary representation of data value
// params:
// size     -   size of the data to print
// ptr      -   pointer to data to print
void print_bits(size_t const size, void const * const ptr){
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i=size-1;i>=0;i--){
        for (j=7;j>=0;j--){
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
}

// concat:
// Concatenate two strings
// params:
// s1       -   first char array
// s2       -   second char array
// return:
// s3       -   char array with s1 followed by s2
char* concat(const char *s1, const char *s2){
    char *result = malloc(strlen(s1) + strlen(s2) + 1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

// starts_with:
// Checks to see if the first string starts with the second string
// params:
// a        -   first string
// b        -   second string
// return:
// x        -   True=1, False=0
int starts_with(const char *a, const char *b) {
   if(strncmp(a, b, strlen(b)) == 0) return 1;
   return 0;
}



// ARC Utility Section
// ###################
// arc_help:
// Prints usage information out for the user
void arc_help(){
    printf("ARC Help not supported at this time...\n");
    return;
}

// arc_resource_init:
// Loads required resource files for encoding
// return:
// x        -   Success=1, Failure=0
int arc_resource_init(){
    printf("ARC Resource Files Initialized\n");
    INIT = 1;
    return 1;
}

// arc_init:
// Loads resource files and training data. If training data is not available, run training tests
// params:
// max_threads  -   Maximum number of threads ARC should use
// return:
// x            -   Success=1, Failure=0
int arc_init(uint32_t max_threads){
    // Load all the required resources files needed for encoding and decoding
    printf("Initializing ARC Resource Files\n");
    int err = arc_resource_init();
    if(err == 0){
        return 0;
    }
    const char* base;
    const char* xdg_runtime_dir = getenv("XDG_RUNTIME_DIR");
    if( xdg_runtime_dir == NULL) {
      const char* tmp_dir = getenv("TMPDIR");
      if(tmp_dir == NULL) {
        base = tmp_dir;
      } else {
        base = "/tmp";
      }
    } else {
      base = xdg_runtime_dir;
    }
    cache_resource_location = concat(base, "/arc");
    if(mkdir(cache_resource_location, S_IRWXU | S_IRWXG) && errno != EEXIST) {
      return 0;
    }

    // Turn off print outs
    PRINT = 0;

    // Initialize Resources for Later Use
    FILE *fp;
    uint32_t i, j, k;
    // Set up training variables
    struct timeval total_start, total_stop;
    struct timeval method_start, method_stop;
    struct timeval start, stop;
    uint8_t *encoded;
    uint32_t encoded_size;
    double encode_time_taken;
    double memory_overhead;
    double throughput_overhead;
    // Set up simulated data set to train with
    uint32_t data_size = 100000;
    uint8_t *data = malloc(sizeof(uint8_t)*data_size);
    srand(0);
    for (i = 0; i < data_size; i++){
        data[i] = (uint8_t)(rand() % 255);
    }

    // Check for valid max_threads
    if (max_threads == ARC_ANY_THREADS || max_threads >= omp_get_max_threads()){
        AVAIL_THREADS = omp_get_max_threads();
    } else if (max_threads <= 0 && max_threads != ARC_ANY_THREADS){
        AVAIL_THREADS = 1;
    } else {
        AVAIL_THREADS = max_threads;
    }

    // Calculate number of configurations per num_threads
    int elements_per_num_threads = 2 + 2 + 128 + 128; // 2 SECDED, 2 Hamming, 128 RS, & 128 Parity Configs
    NUM_CONFIGURATIONS = AVAIL_THREADS * elements_per_num_threads;
    // Initialize arc_configurations array with enough space for each configuration
    arc_configurations = malloc(NUM_CONFIGURATIONS * sizeof(*arc_configurations));
    int current_config_array_index = 0;

    // For threads 1 to max_threads, check to see if we have the training data for these already
    printf("ARC Training and Loading Process Started\n");
    gettimeofday(&total_start, NULL);
    for(i = 1; i <= AVAIL_THREADS; i++){
        // Check resource folder for corresponding num_threads information cache
        char num_thread_string[8];
        sprintf(num_thread_string, "%d", i);
        char * thread_file_name = concat(num_thread_string, thread_resource_file);
        char * thread_file_location = concat(cache_resource_location, thread_file_name);
        fp = fopen(thread_file_location, "r");

        printf("\n%d Thread(s) Training and Loading Started...\n", i);
        // If the file does not exist,
        if(fp == NULL){

            // Open a write pointer
            fp = fopen(thread_file_location, "w");
            fprintf(fp, "ecc_algorithm,ecc_parameter_a,ecc_parameter_b,num_threads,memory_overhead,throughput_overhead\n");

            // Run training tests to populate both the arc_configurations array and file
            for(j = 1; j <= NUM_METHODS; j++){
                // Parity Training
                if(j == PARITY_ID){
                    // Start Training
                    printf("Parity Training Started...\n");
                    gettimeofday(&method_start, NULL);

                    // Train 1 parity bits over 1 byte to 128 bytes of data
                    for(k = 1; k <= 128; k++){
                        // Time the encoding process
                        gettimeofday(&start, NULL);
                        err = arc_parity_encode(data, data_size, k, i, &encoded, &encoded_size);
                        gettimeofday(&stop, NULL);
                        encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                        // Calculate Overhead Factors
                        memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                        throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                        // Add this information to configuration_information struct array
                        arc_configurations[current_config_array_index].ecc_algorithm = j;
                        arc_configurations[current_config_array_index].ecc_parameter_a = (uint32_t)k;
                        arc_configurations[current_config_array_index].ecc_parameter_b = (uint32_t)0;
                        arc_configurations[current_config_array_index].num_threads = (uint32_t)i;
                        arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                        arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                        // Write these results out to corresponding thread training file
                        fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, k, 0, i, memory_overhead, throughput_overhead);
                        fflush(fp);
                        // Increment current config array index
                        current_config_array_index++;
                        // Free encoded data
                        free(encoded);
                    }

                    // End Training
                    gettimeofday(&method_stop, NULL);
                    double method_time_taken = (double)(method_stop.tv_usec - method_start.tv_usec) / 1000000 + (double)(method_stop.tv_sec - method_start.tv_sec);
                    printf("Parity Training Completed in %lf Seconds!\n", method_time_taken);
                // Hamming Training 
                } else if(j == HAMMING_ID){  
                    // Start Training
                    printf("Hamming Training Started...\n");
                    gettimeofday(&method_start, NULL);

                    // Train hamming over 1 byte of data
                    uint32_t hamming_block_size = 1;
                    gettimeofday(&start, NULL);
                    err = arc_hamming_encode(data, data_size, hamming_block_size, i, &encoded, &encoded_size);
                    gettimeofday(&stop, NULL);
                    encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                    // Calculate Overhead Factors
                    memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                    throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                    // Add this information to configuration_information struct array
                    arc_configurations[current_config_array_index].ecc_algorithm = j;
                    arc_configurations[current_config_array_index].ecc_parameter_a = hamming_block_size;
                    arc_configurations[current_config_array_index].ecc_parameter_b = 0;
                    arc_configurations[current_config_array_index].num_threads = i;
                    arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                    arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                    // Write these results out to corresponding thread training file
                    fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, hamming_block_size, 0, i, memory_overhead, throughput_overhead);
                    fflush(fp);
                    // Increment current config array index
                    current_config_array_index++;
                    // Free encoded data
                    free(encoded);

                    // Train hamming over 8 bytes of data
                    hamming_block_size = 8;
                    gettimeofday(&start, NULL);
                    err = arc_hamming_encode(data, data_size, hamming_block_size, i, &encoded, &encoded_size);
                    gettimeofday(&stop, NULL);
                    encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                    // Calculate Overhead Factors
                    memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                    throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                    // Add this information to configuration_information struct array
                    arc_configurations[current_config_array_index].ecc_algorithm = j;
                    arc_configurations[current_config_array_index].ecc_parameter_a = hamming_block_size;
                    arc_configurations[current_config_array_index].ecc_parameter_b = 0;
                    arc_configurations[current_config_array_index].num_threads = i;
                    arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                    arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                    // Write these results out to corresponding thread training file
                    fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, hamming_block_size, 0, i, memory_overhead, throughput_overhead);
                    fflush(fp);
                    // Increment current config array index
                    current_config_array_index++;
                    // Free encoded data
                    free(encoded);

                    // End Training 
                    gettimeofday(&method_stop, NULL);
                    double method_time_taken = (double)(method_stop.tv_usec - method_start.tv_usec) / 1000000 + (double)(method_stop.tv_sec - method_start.tv_sec);
                    printf("Hamming Training Completed in %lf Seconds!\n", method_time_taken);

                // SECDED Training
                } else if(j == SECDED_ID){
                    // Start Training
                    printf("SECDED Training Started...\n");
                    gettimeofday(&method_start, NULL);

                    // Train secded over 1 byte of data
                    uint32_t secded_block_size = 1;
                    gettimeofday(&start, NULL);
                    err = arc_secded_encode(data, data_size, secded_block_size, i, &encoded, &encoded_size);
                    gettimeofday(&stop, NULL);
                    encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                    // Calculate Overhead Factors
                    memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                    throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                    // Add this information to configuration_information struct array
                    arc_configurations[current_config_array_index].ecc_algorithm = j;
                    arc_configurations[current_config_array_index].ecc_parameter_a = secded_block_size;
                    arc_configurations[current_config_array_index].ecc_parameter_b = 0;
                    arc_configurations[current_config_array_index].num_threads = i;
                    arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                    arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                    // Write these results out to corresponding thread training file
                    fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, secded_block_size, 0, i, memory_overhead, throughput_overhead);
                    fflush(fp);
                    // Increment current config array index
                    current_config_array_index++;
                    // Free encoded data
                    free(encoded);

                    // Train secded over 8 bytes of data
                    secded_block_size = 8;
                    gettimeofday(&start, NULL);
                    err = arc_secded_encode(data, data_size, secded_block_size, i, &encoded, &encoded_size);
                    gettimeofday(&stop, NULL);
                    encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                    // Calculate Overhead Factors
                    memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                    throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                    // Add this information to configuration_information struct array
                    arc_configurations[current_config_array_index].ecc_algorithm = j;
                    arc_configurations[current_config_array_index].ecc_parameter_a = secded_block_size;
                    arc_configurations[current_config_array_index].ecc_parameter_b = 0;
                    arc_configurations[current_config_array_index].num_threads = i;
                    arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                    arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                    // Write these results out to corresponding thread training file
                    fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, secded_block_size, 0, i, memory_overhead, throughput_overhead);
                    fflush(fp);
                    // Increment current config array index
                    current_config_array_index++;
                    // Free encoded data
                    free(encoded);

                    // End Training 
                    gettimeofday(&method_stop, NULL);
                    double method_time_taken = (double)(method_stop.tv_usec - method_start.tv_usec) / 1000000 + (double)(method_stop.tv_sec - method_start.tv_sec);
                    printf("SECDED Training Completed in %lf Seconds!\n", method_time_taken);

                // RS Training
                } else if(j == RS_ID){
                    // Start Training
                    printf("Reed-Solomon Training Started...\n");
                    printf("(This may take some time)\n");
                    gettimeofday(&method_start, NULL);

                    // Set initial data and code devices
                    uint32_t data_devices = 128;
                    uint32_t code_devices = 128;

                    // Train every 4th configuration and estimate the others
                    double previous_memory_overhead = 0;
                    double previous_throughput_overhead = 0;
                    uint32_t tmp_data_devices;
                    uint32_t tmp_code_devices;
                    for (k = 0; k <= 32; k++){
                        // Time the encoding process
                        gettimeofday(&start, NULL);
                        err = arc_reed_solomon_encode(data, data_size, data_devices, code_devices, i, &encoded, &encoded_size);
                        gettimeofday(&stop, NULL);
                        encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
                        // Calculate Overhead Factors
                        memory_overhead = ((double)encoded_size / (double)data_size) - 1;
                        throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s

                        // Estimate missing values and save current values for next estimation
                        if (k == 0){
                            // Store current calculations for next round of estimations
                            previous_memory_overhead = memory_overhead;
                            previous_throughput_overhead = throughput_overhead;
                            // Add this information to configuration_information struct array
                            arc_configurations[current_config_array_index].ecc_algorithm = j;
                            arc_configurations[current_config_array_index].ecc_parameter_a = data_devices;
                            arc_configurations[current_config_array_index].ecc_parameter_b = code_devices;
                            arc_configurations[current_config_array_index].num_threads = i;
                            arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                            arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                            // Write these results out to corresponding thread training file
                            fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, data_devices, code_devices, i, memory_overhead, throughput_overhead);
                            fflush(fp);
                            // Increment current config array index
                            current_config_array_index++;
                        } else {
                            // Estimate the mid point values
                            double mid_memory_overhead = (previous_memory_overhead + memory_overhead) / (double) 2;
                            double mid_throughput_overhead = (previous_throughput_overhead + throughput_overhead) / (double) 2;
                            // Estimate the first quartile point values
                            double first_quartile_memory_overhead = (previous_memory_overhead + mid_memory_overhead) / (double) 2;
                            double first_quartile_throughput_overhead = (previous_throughput_overhead + mid_throughput_overhead) / (double) 2;
                            // Estimate the third quartile point values
                            double third_quartile_memory_overhead = (memory_overhead + mid_memory_overhead) / (double) 2;
                            double third_quartile_throughput_overhead = (throughput_overhead + mid_throughput_overhead) / (double) 2;

                            // Store the newly estimated values in the arc configurations table
                            // First estimated value
                            tmp_data_devices = data_devices - 3;
                            tmp_code_devices = code_devices + 3;
                            arc_configurations[current_config_array_index].ecc_algorithm = j;
                            arc_configurations[current_config_array_index].ecc_parameter_a = tmp_data_devices;
                            arc_configurations[current_config_array_index].ecc_parameter_b = tmp_code_devices;
                            arc_configurations[current_config_array_index].num_threads = i;
                            arc_configurations[current_config_array_index].memory_overhead  = first_quartile_memory_overhead;
                            arc_configurations[current_config_array_index].throughput_overhead  = first_quartile_throughput_overhead;
                            // Write these results out to corresponding thread training file
                            fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, tmp_data_devices, tmp_code_devices, i, first_quartile_memory_overhead, first_quartile_throughput_overhead);
                            fflush(fp);
                            // Increment current config array index
                            current_config_array_index++;

                            // Dont include the mid point of the final iteration as this is not actually a valid configuration
                            // i.e: 253.5 data devices is not valid
                            if (k != 32){
                                // Second estimated value
                                tmp_data_devices = data_devices - 2;
                                tmp_code_devices = code_devices + 2;
                                arc_configurations[current_config_array_index].ecc_algorithm = j;
                                arc_configurations[current_config_array_index].ecc_parameter_a = tmp_data_devices;
                                arc_configurations[current_config_array_index].ecc_parameter_b = tmp_code_devices;
                                arc_configurations[current_config_array_index].num_threads = i;
                                arc_configurations[current_config_array_index].memory_overhead  = mid_memory_overhead;
                                arc_configurations[current_config_array_index].throughput_overhead  = mid_throughput_overhead;
                                // Write these results out to corresponding thread training file
                                fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, tmp_data_devices, tmp_code_devices, i, mid_memory_overhead, mid_throughput_overhead);
                                fflush(fp);
                                // Increment current config array index
                                current_config_array_index++;
                            }
                            
                            // Third estimated value
                            tmp_data_devices = data_devices - 1;
                            tmp_code_devices = code_devices + 1;
                            arc_configurations[current_config_array_index].ecc_algorithm = j;
                            arc_configurations[current_config_array_index].ecc_parameter_a = tmp_data_devices;
                            arc_configurations[current_config_array_index].ecc_parameter_b = tmp_code_devices;
                            arc_configurations[current_config_array_index].num_threads = i;
                            arc_configurations[current_config_array_index].memory_overhead  = third_quartile_memory_overhead;
                            arc_configurations[current_config_array_index].throughput_overhead  = third_quartile_throughput_overhead;
                            // Write these results out to corresponding thread training file
                            fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, tmp_data_devices, tmp_code_devices, i, third_quartile_memory_overhead, third_quartile_throughput_overhead);
                            fflush(fp);
                            // Increment current config array index
                            current_config_array_index++;

                            // Store current calculations for next round of estimations
                            previous_memory_overhead = memory_overhead;
                            previous_throughput_overhead = throughput_overhead;
                            // Add this information to configuration_information struct array
                            arc_configurations[current_config_array_index].ecc_algorithm = j;
                            arc_configurations[current_config_array_index].ecc_parameter_a = data_devices;
                            arc_configurations[current_config_array_index].ecc_parameter_b = code_devices;
                            arc_configurations[current_config_array_index].num_threads = i;
                            arc_configurations[current_config_array_index].memory_overhead  = memory_overhead;
                            arc_configurations[current_config_array_index].throughput_overhead  = throughput_overhead;
                            // Write these results out to corresponding thread training file
                            fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", j, data_devices, code_devices, i, memory_overhead, throughput_overhead);
                            fflush(fp);
                            // Increment current config array index
                            current_config_array_index++;
                        }

                        // Change data and code devices for next iteration
                        if (data_devices == 252){
                            data_devices = 255;
                            code_devices = 1;
                        } else {
                            data_devices = data_devices + 4;
                            code_devices = code_devices - 4;
                        }
                        // Free encoded data
                        free(encoded);
                    }

                    // End Training 
                    gettimeofday(&method_stop, NULL);
                    double method_time_taken = (double)(method_stop.tv_usec - method_start.tv_usec) / 1000000 + (double)(method_stop.tv_sec - method_start.tv_sec);
                    printf("Reed-Solomon Training Completed in %lf Seconds!\n", method_time_taken);
                }
            }
            // Free strings
            free(thread_file_location);
            free(thread_file_name);
        // If the file does exist, pull in information from file
        } else {
            int tmp_alg;
            uint32_t tmp_a;
            uint32_t tmp_b;
            uint32_t tmp_th;
            double tmp_mo;
            double tmp_to;
            char *buffer = NULL;
            size_t len = 0;
            ssize_t read;
            while ((read = getline(&buffer, &len, fp)) != -1) {
                // Skip the first line
                if(!starts_with(buffer, "ecc_algorithm")){
                    // Scan each line and pull the specific values
                    sscanf(buffer, "%d,%" SCNu32 ",%" SCNu32 ",%" SCNu32 ",%lf,%lf", &tmp_alg, &tmp_a, &tmp_b, &tmp_th, &tmp_mo, &tmp_to);
                    // Store these values in the configurations array
                    arc_configurations[current_config_array_index].ecc_algorithm = tmp_alg;
                    arc_configurations[current_config_array_index].ecc_parameter_a = tmp_a;
                    arc_configurations[current_config_array_index].ecc_parameter_b = tmp_b;
                    arc_configurations[current_config_array_index].num_threads = tmp_th;
                    arc_configurations[current_config_array_index].memory_overhead  = tmp_mo;
                    arc_configurations[current_config_array_index].throughput_overhead  = tmp_to;
                    // Increment current config array index
                    current_config_array_index++;
                }
            }
            fclose(fp);
            if (buffer){
                free(buffer);
            }
        }
        printf("%d Thread(s) Training and Loading Completed!\n", i);
    }
    // Return once everything has been loaded
    gettimeofday(&total_stop, NULL);
    double total_time_taken = (double)(total_stop.tv_usec - total_start.tv_usec) / 1000000 + (double)(total_stop.tv_sec - total_start.tv_sec);
    printf("ARC Training and Loading Process Finished!\n");
    printf("Total Training Time Taken: %lf Seconds\n", total_time_taken);

    // Turn on print outs
    PRINT = 1;
    // Free allocated simulated data
    free(data);

    return 1;
}

// arc_save:
// Saves training data collected during normal use to resource files. This process helps 
// ARC become smarter over time and should be done before arc_close is called.
// return:
// x            -   Success=1, Failure=0
int arc_save(){
    // Initialize Resources for Later Use
    int i, j;
    FILE *fp;

    // For each number of available threads, write the updated values to the file
    for(i = 1; i <= AVAIL_THREADS; i++){
        // Create the pointer to the file
        char num_thread_string[8];
        sprintf(num_thread_string, "%d", i);
        char * thread_file_name = concat(num_thread_string, thread_resource_file);
        char * thread_file_location = concat(cache_resource_location, thread_file_name);
        fp = fopen(thread_file_location, "w");
        fprintf(fp, "ecc_algorithm,ecc_parameter_a,ecc_parameter_b,num_threads,memory_overhead,throughput_overhead\n");

        // Iterate over all configurations and write the corresponding thread configurations to file
        for (j = 0; j < NUM_CONFIGURATIONS; j++){
            if (arc_configurations[j].num_threads == i){
                fprintf(fp, "%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%lf,%lf\n", arc_configurations[j].ecc_algorithm, arc_configurations[j].ecc_parameter_a, arc_configurations[j].ecc_parameter_b, arc_configurations[j].num_threads, arc_configurations[j].memory_overhead, arc_configurations[j].throughput_overhead);
            }
        }

        // Close file once all configurations have been written
        fclose(fp);
    }

    return 1;
}

// arc_close:
// Frees allocated resources needed for ARC to be used
// return:
// x            -   Success=1, Failure=0
int arc_close(){
    // Free arc configurations list
    free(arc_configurations);
    // Set init back to false 
    INIT = 0;
    printf("ARC successfully closed\n");
    return 1;
}

// ARC Main Section
// ################
// arc_encode:
// Encode given data using best ECC approach based on given constraints
// params:
// data                     -   uint8_t data stream
// data_size                -   size of data stream
// memory_constraint        -   maximum amount of memory overhead to introduce (1 - (encoded_size/original_size))
// throughput_constraint    -   minimum bandwidth ARC should have in MB/s (data_size / encode_time_taken)/1E6)
// resiliency_constraint    -   list of ECC methods to choose from / number of errors predicted to occur per MB of data
// resiliency_count         -   Number of values set in resiliency constraint list
// encoded_data             -   address of pointer to uint8_t arc encoded data stream
// encoded_data_size        -   address of pointer to size of arc encoded data stream
// return:
// err                      -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_encode(uint8_t* data, uint32_t data_size, double memory_constraint, double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint8_t** encoded_data, uint32_t* encoded_data_size){
    // Ensure initialization was called first
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    // Check that all constraints are valid
    if (memory_constraint <= 0 && memory_constraint != ARC_ANY_SIZE){
        printf("INVALID MEMORY CONSTRAINT: Please ensure memory constraint is ARC_ANY_SIZE or greater than 0. . .\n");
        return 0;
    }
    if (throughput_constraint <= 0 && throughput_constraint != ARC_ANY_BW){
        printf("INVALID TIME CONSTRAINT: Please ensure timing constraint is ARC_ANY_BW greater than 0. . .\n");
        return 0;
    }
    if (resiliency_count <= 0){
        printf("INVALID RESILIENCY CONSTRAINT: Please ensure resiliency constraint is at least [ARC_ANY_ECC] and resiliency_count is greater than 0. . .\n");
        return 0;
    }
    printf("ARC Encoding Started\n");

    // Given user constraints determine which encoding mode to use
    int err;
    int optimizer_choice;
    uint32_t optimizer_parameter_a;
    uint32_t optimizer_parameter_b;
    uint32_t num_threads;

    // Doesn't care about either
    if (memory_constraint == ARC_ANY_SIZE && throughput_constraint == ARC_ANY_BW){
        // When neither are required, use the method with highest memory overhead (aka most code bits)
        num_threads = 1;
        optimizer_choice = arc_memory_optimizer(1000000, resiliency_constraint, resiliency_count, &optimizer_parameter_a, &optimizer_parameter_b);
    // Only cares about storage
    } else if (memory_constraint != ARC_ANY_SIZE && throughput_constraint == ARC_ANY_BW) {
        num_threads = 1;
        optimizer_choice = arc_memory_optimizer(memory_constraint, resiliency_constraint, resiliency_count, &optimizer_parameter_a, &optimizer_parameter_b);
    // Only cares about time
    } else if (memory_constraint == ARC_ANY_SIZE && throughput_constraint != ARC_ANY_BW) {
        optimizer_choice = arc_throughput_optimizer(throughput_constraint, resiliency_constraint, resiliency_count, &optimizer_parameter_a, &optimizer_parameter_b, &num_threads);
    // Cares about both
    } else if (memory_constraint != ARC_ANY_SIZE && throughput_constraint != ARC_ANY_BW){
        optimizer_choice = arc_joint_optimizer(memory_constraint, throughput_constraint, resiliency_constraint, resiliency_count, &optimizer_parameter_a, &optimizer_parameter_b, &num_threads);
    } else {
        printf("Optimizer Error: This should not occur...\n");
        return 0;
    }

    // Given the optimizer choice, encode with the correct scheme
    struct timeval start, stop;
    if (optimizer_choice == PARITY_ID){
        printf("Utilizing Parity-%" PRIu32 " on %" PRIu32 " threads\n", optimizer_parameter_a, num_threads);
        gettimeofday(&start, NULL);
        err = arc_parity_encode(data, data_size, optimizer_parameter_a, num_threads, encoded_data, encoded_data_size);
        gettimeofday(&stop, NULL);
    } else if (optimizer_choice == HAMMING_ID){
        printf("Utilizing Hamming-%" PRIu32 " on %" PRIu32 " threads\n", optimizer_parameter_a, num_threads);
        gettimeofday(&start, NULL);
        err = arc_hamming_encode(data, data_size, optimizer_parameter_a, num_threads, encoded_data, encoded_data_size);
        gettimeofday(&stop, NULL);
    } else if (optimizer_choice == SECDED_ID){
        printf("Utilizing SECDED-%"PRIu32" on %" PRIu32 " threads\n", optimizer_parameter_a, num_threads);
        gettimeofday(&start, NULL);
        err = arc_secded_encode(data, data_size, optimizer_parameter_a, num_threads, encoded_data, encoded_data_size);
        gettimeofday(&stop, NULL);
    } else if (optimizer_choice == RS_ID){
        printf("Utilizing Reed Solomon-<k=%" PRIu32 ",m=%" PRIu32 ",w=%" PRIu32 "> on %" PRIu32 " threads\n", optimizer_parameter_a, optimizer_parameter_b, 8, num_threads);
        gettimeofday(&start, NULL);
        err = arc_reed_solomon_encode(data, data_size, optimizer_parameter_a, optimizer_parameter_b, num_threads, encoded_data, encoded_data_size);
        gettimeofday(&stop, NULL);
    } else {
        printf("Encoding Error: No encoding function meets required criteria\n");
        return 0;
    }
    // Check for encoding errors
    if (err == 0){
        printf("Encoding Error: Error has occurred during encoding process\n");
        return 0;
    }

    // Get throughput of the encoding process
    double encode_time_taken = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
    double throughput_overhead = (data_size / encode_time_taken)/1E6; // MB/s
                        
    // Update the corresponding configuration in configuration array
    int i;
    for (i = 0; i < NUM_CONFIGURATIONS; i++){
        // Update configuration that was used
        if (arc_configurations[i].ecc_algorithm == optimizer_choice && arc_configurations[i].num_threads == num_threads 
        && arc_configurations[i].ecc_parameter_a == optimizer_parameter_a && arc_configurations[i].ecc_parameter_b == optimizer_parameter_b){
            arc_configurations[i].throughput_overhead = (arc_configurations[i].throughput_overhead + throughput_overhead) / 2;
        }
    }

    // Return data
    printf("ARC Encoding Completed\n");
    return 1;
}

// arc_decode:
// Decode the given encoded data using the correct ECC approach
// params:
// encoded_data         -   uint8_t arc encoded data stream
// encoded_data_size    -   size of arc encoded data stream
// decoded_data         -   address of pointer to arc decoded uint8_t data stream
// decoded_data_size    -   address of pointer to arc decoded data stream size
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** decoded_data, uint32_t* decoded_data_size){
    // Ensure initialization was called first
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    printf("ARC Decoding Started\n");

    // Determine which method was used  
    uint8_t encoding_method = encoded_data[0];

    // Decompress with the correct decode method
    int err;
    // If it was encoded with parity
    if (encoding_method == 0x01){
        printf("Parity Encoding Found\n");
        err = arc_parity_decode(encoded_data, encoded_data_size, decoded_data, decoded_data_size);
    // If it was encoded with hamming
    } else if (encoding_method == 0x02){
        printf("Hamming Encoding Found\n");
        err = arc_hamming_decode(encoded_data, encoded_data_size, decoded_data, decoded_data_size);
    // If it was encoded with secded
    } else if (encoding_method == 0x03){
        printf("SECDED Encoding Found\n");
        err = arc_secded_decode(encoded_data, encoded_data_size, decoded_data, decoded_data_size);
    } else if (encoding_method == 0x04){
        printf("Reed Solomon Encoding Found\n");
        err = arc_reed_solomon_decode(encoded_data, encoded_data_size, decoded_data, decoded_data_size);
    } else {
        printf("INVALID ENCODING METHOD: No valid encoding method found to decode this data. . .\n");
        return 0;
    }

    // Return data
    printf("ARC Decoding Completed\n");
    return err;
}

// arc_check:
// Check the parity bits of the encoded data stream to ensure integrity
// params:
// encoded_data         -   uint8_t protected data stream
// encoded_data_size    -   size of protected data stream
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_check(uint8_t* encoded_data, uint32_t encoded_data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 1;
    };
    printf("ARC Check not supported at this time...\n");
    return 1;
}



// ARC Engine Section
// ###################
// arc_memory_optimizer:
// Given a memory constraint, determine optimal ECC approach
// params:
// memory_constraint            -   maximum amount of memory overhead to introduce (1 - (encoded_size/original_size))
// resiliency_constraint        -   list of ECC methods to choose from / number of errors predicted to occur per MB of data
// resiliency_count             -   Number of values set in resiliency constraint list
// memory_choice_parameter_a    -   first resulting ecc method parameter
// memory_choice_parameter_b    -   second resulting ecc method parameter (if applicable, 0 otherwise)
// return:
// identifier                   -   identifier to determine which ecc method to use
int arc_memory_optimizer(double memory_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* memory_choice_parameter_a, uint32_t* memory_choice_parameter_b){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i, j;
    // Initialize Helper Variables
    int optimal_index = -1;
    double optimal_overhead = 0;
    int minimal_index = -1;
    double minimal_overhead = 0;
    int *ecc_options = (int *)malloc(NUM_METHODS * sizeof(int));
    for (i = 0; i < NUM_METHODS; i++){
        ecc_options[i] = 0;
    }
    
    // Determine viable ECC methods
    for (i = 0; i < resiliency_count; i++){
        // If ARC_ANY_ECC, set all to true and break out of loop
        if (resiliency_constraint[i] == ARC_ANY_ECC){
            for (j = 0; j < NUM_METHODS; j++){
                ecc_options[j] = 1;
            }
            break;
        // If one of ARC_* ECC constants, set the corresponding location to true
        } else if (resiliency_constraint[i] < 0){
            // Enable Parity
            if (resiliency_constraint[i] == -1){
                ecc_options[0] = 1;
            // Enable Hamming
            } else if (resiliency_constraint[i] == -2){
                ecc_options[1] = 1;
            // Enable SECDED
            } else if (resiliency_constraint[i] == -3){
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -4){
                ecc_options[3] = 1;
            // Enable Parity, Hamming, and SECDED
            } else if (resiliency_constraint[i] == -5){
                ecc_options[0] = 1;
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable Hamming and SECDED
            } else if (resiliency_constraint[i] == -6){
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -7){
                ecc_options[3] = 1;
            }
        // If in terms of Errors/MB turn on corresponding ECC methods
        } else {
            // If less than an eighth of MB is errors, set two highest ECC
            if (resiliency_constraint[i] < 125000){
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
                // Enable SECDED
                ecc_options[NUM_METHODS - 2] = 1;
            // If more than an eighth of MB is errors, set only highest ECC
            } else {
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
            }

        }
    }

    // Check each available configuration for highest memory overhead that's below the given memory constraint
    for (i = 0; i < NUM_CONFIGURATIONS; i++){
        // Ensure ECC algorithm is viable given resiliency constraint
        if (ecc_options[(arc_configurations[i].ecc_algorithm - 1)] == 1){
            // Look for optimal approach
            if ((optimal_index == -1 && arc_configurations[i].memory_overhead <= memory_constraint) || 
                (arc_configurations[i].memory_overhead > optimal_overhead && arc_configurations[i].memory_overhead <= memory_constraint)){
                optimal_index = i;
                optimal_overhead = arc_configurations[i].memory_overhead;
            }

            // Look for minimal approach as a backup
            if ((minimal_index == -1) || (arc_configurations[i].memory_overhead <= minimal_overhead)){
                minimal_index = i;
                minimal_overhead = arc_configurations[i].memory_overhead;
            }
        }
    }
    // Free ECC Options array
    free(ecc_options);

    // Check to see if a valid configuration was found
    if (optimal_index == -1){
        printf("ARC Memory Optimization Warning: No such configuration satisfies provided memory constraint\n");
        printf("Using the lowest possible memory overhead approach possible. . .\n");
        // Set the parameters to the minimal index parameters
        *memory_choice_parameter_a = arc_configurations[minimal_index].ecc_parameter_a;
        *memory_choice_parameter_b = arc_configurations[minimal_index].ecc_parameter_b;
        return arc_configurations[minimal_index].ecc_algorithm;
    } else {
        // Set the parameters to the optimal index parameters
        *memory_choice_parameter_a = arc_configurations[optimal_index].ecc_parameter_a;
        *memory_choice_parameter_b = arc_configurations[optimal_index].ecc_parameter_b;
        return arc_configurations[optimal_index].ecc_algorithm;
    }
}

// arc_throughput_optimizer:
// Given a throughput constraint, determine optimal ECC approach
// params:
// throughput_constraint            -   minimum bandwidth ARC should have in MB/s (data_size / encode_time_taken)/1E6)
// resiliency_constraint            -   list of ECC methods to choose from / number of errors predicted to occur per MB of data
// resiliency_count                 -   Number of values set in resiliency constraint list
// throughput_choice_parameter_a    -   first resulting ecc method parameter
// throughput_choice_parameter_b    -   second resulting ecc method parameter (if applicable, 0 otherwise)
// num_threads                      -   number of threads required to obtain optimal bandwidth
// return:
// identifier                       -   identifier to determine which ecc method to use
int arc_throughput_optimizer(double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* throughput_choice_parameter_a, uint32_t* throughput_choice_parameter_b, uint32_t* num_threads){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i, j;

    // Initialize Helper Variables
    int optimal_index = -1;
    double optimal_overhead = 0;
    int maximum_index = -1;
    double maximum_overhead = 0;
    int *ecc_options = (int *)malloc(NUM_METHODS * sizeof(int));
    for (i = 0; i < NUM_METHODS; i++){
        ecc_options[i] = 0;
    }
    
    // Determine viable ECC methods
    for (i = 0; i < resiliency_count; i++){
        // If ARC_ANY_ECC, set all to true and break out of loop
        if (resiliency_constraint[i] == ARC_ANY_ECC){
            for (j = 0; j < NUM_METHODS; j++){
                ecc_options[j] = 1;
            }
            break;
        // If one of ARC_* ECC constants, set the corresponding location to true
        } else if (resiliency_constraint[i] < 0){
            // Enable Parity
            if (resiliency_constraint[i] == -1){
                ecc_options[0] = 1;
            // Enable Hamming
            } else if (resiliency_constraint[i] == -2){
                ecc_options[1] = 1;
            // Enable SECDED
            } else if (resiliency_constraint[i] == -3){
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -4){
                ecc_options[3] = 1;
            // Enable Parity, Hamming, and SECDED
            } else if (resiliency_constraint[i] == -5){
                ecc_options[0] = 1;
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable Hamming and SECDED
            } else if (resiliency_constraint[i] == -6){
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -7){
                ecc_options[3] = 1;
            }
        // If in terms of Errors/MB turn on corresponding ECC methods
        } else {
            // If less than an eighth of MB is errors, set two highest ECC
            if (resiliency_constraint[i] < 125000){
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
                // Enable SECDED
                ecc_options[NUM_METHODS - 2] = 1;
            // If more than an eighth of MB is errors, set only highest ECC
            } else {
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
            }

        }
    }

    // Check each available configuration for lowest throughput that's above the given throughput constraint
    for (i = 0; i < NUM_CONFIGURATIONS; i++){
        // Ensure ECC algorithm is viable given resiliency constraint
        if (ecc_options[(arc_configurations[i].ecc_algorithm - 1)] == 1){
            // Look for optimal approach
            if ((optimal_index == -1 && arc_configurations[i].throughput_overhead >= throughput_constraint) || 
                (arc_configurations[i].throughput_overhead >= throughput_constraint && arc_configurations[i].throughput_overhead < optimal_overhead)){
                optimal_index = i;
                optimal_overhead = arc_configurations[i].throughput_overhead;
            }

            // Look for maximum bandwidth approach as a backup
            if ((maximum_index == -1) || (arc_configurations[i].throughput_overhead >= maximum_overhead)){
                maximum_index = i;
                maximum_overhead = arc_configurations[i].throughput_overhead;
            }
        }
    }
    // Free ECC Options array
    free(ecc_options);

    // Check to see if a valid configuration was found
    if (optimal_index == -1){
        printf("ARC Throughput Optimization Warning: No such configuration satisfies provided throughput constraint\n");
        printf("Using the highest possible bandwidth approach possible. . .\n");
        // Set the parameters to the maximum index parameters
        *throughput_choice_parameter_a = arc_configurations[maximum_index].ecc_parameter_a;
        *throughput_choice_parameter_b = arc_configurations[maximum_index].ecc_parameter_b;
        *num_threads = arc_configurations[maximum_index].num_threads;
        return arc_configurations[maximum_index].ecc_algorithm;
    } else {
        // Set the parameters to the optimal index parameters
        *throughput_choice_parameter_a = arc_configurations[optimal_index].ecc_parameter_a;
        *throughput_choice_parameter_b = arc_configurations[optimal_index].ecc_parameter_b;
        *num_threads = arc_configurations[optimal_index].num_threads;
        return arc_configurations[optimal_index].ecc_algorithm;
    }
}

// arc_joint_optimizer:
// Given memory and throughput constraints, determine optimal ECC approach
// params:
// memory_constraint                -   maximum amount of memory overhead to introduce (1 - (encoded_size/original_size))
// throughput_constraint            -   minimum bandwidth ARC should have in MB/s (data_size / encode_time_taken)/1E6)
// resiliency_constraint            -   list of ECC methods to choose from / number of errors predicted to occur per MB of data
// resiliency_count                 -   Number of values set in resiliency constraint list
// throughput_choice_parameter_a    -   first resulting ecc method parameter
// throughput_choice_parameter_b    -   second resulting ecc method parameter (if applicable, 0 otherwise)
// num_threads                      -   number of threads required to obtain optimal bandwidth
// return:
// identifier                       -   identifier to determine which ecc method to use
int arc_joint_optimizer(double memory_constraint, double throughput_constraint, int *resiliency_constraint, int resiliency_count, uint32_t* optimizer_parameter_a, uint32_t* optimizer_parameter_b, uint32_t* num_threads){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i, j;

    // Initialize ECC Configurations Array
    int *ecc_options = (int *)malloc(NUM_METHODS * sizeof(int));
    for (i = 0; i < NUM_METHODS; i++){
        ecc_options[i] = 0;
    }
    
    // Determine viable ECC methods
    for (i = 0; i < resiliency_count; i++){
        // If ARC_ANY_ECC, set all to true and break out of loop
        if (resiliency_constraint[i] == ARC_ANY_ECC){
            for (j = 0; j < NUM_METHODS; j++){
                ecc_options[j] = 1;
            }
            break;
        // If one of ARC_* ECC constants, set the corresponding location to true
        } else if (resiliency_constraint[i] < 0){
            // Enable Parity
            if (resiliency_constraint[i] == -1){
                ecc_options[0] = 1;
            // Enable Hamming
            } else if (resiliency_constraint[i] == -2){
                ecc_options[1] = 1;
            // Enable SECDED
            } else if (resiliency_constraint[i] == -3){
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -4){
                ecc_options[3] = 1;
            // Enable Parity, Hamming, and SECDED
            } else if (resiliency_constraint[i] == -5){
                ecc_options[0] = 1;
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable Hamming and SECDED
            } else if (resiliency_constraint[i] == -6){
                ecc_options[1] = 1;
                ecc_options[2] = 1;
            // Enable RS
            } else if (resiliency_constraint[i] == -7){
                ecc_options[3] = 1;
            }
        // If in terms of Errors/MB turn on corresponding ECC methods
        } else {
            // If less than an eighth of MB is errors, set two highest ECC
            if (resiliency_constraint[i] < 125000){
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
                // Enable SECDED
                ecc_options[NUM_METHODS - 2] = 1;
            // If more than an eighth of MB is errors, set only highest ECC
            } else {
                // Enable Reed-Solomon
                ecc_options[NUM_METHODS - 1] = 1;
            }

        }
    }

    // Create two arrays of -1's of length NUM_CONFIGURATIONS
    int possible_throughput_configurations[NUM_CONFIGURATIONS];
    int possible_memory_configurations[NUM_CONFIGURATIONS];
    int possible_throughput_config_index = 0;
    int possible_memory_config_index = 0;
    // Check each available configuration and gather all that satisfy each constraint individually
    for (i = 0; i < NUM_CONFIGURATIONS; i++){
        // Set the current index of each array to -1 by default
        possible_throughput_configurations[i] = -1;
        possible_memory_configurations[i] = -1;

        // Ensure ECC algorithm is viable given resiliency constraint
        if (ecc_options[(arc_configurations[i].ecc_algorithm - 1)] == 1){
            // Check to see if it satisfies throughput constraint
            if (arc_configurations[i].throughput_overhead >= throughput_constraint){
                possible_throughput_configurations[possible_throughput_config_index] = i;
                possible_throughput_config_index++;
            }
            // Check to see if it satisfies memory constraint
            if (arc_configurations[i].memory_overhead <= memory_constraint){
                possible_memory_configurations[possible_memory_config_index] = i;
                possible_memory_config_index++;
            }
        }
    }

    // If there are some configurations that satisfy both constraints,
    if (possible_throughput_config_index != 0 && possible_memory_config_index != 0){
        // Iterate through all the possible throughput satisfying configurations to find the best memory configuration
        int optimal_index = -1;
        double optimal_overhead = 0;
        int minimal_index = -1;
        double minimal_overhead = 0;

        for (i = 0; i < NUM_CONFIGURATIONS; i++){
            // Look for optimal approach
            if (possible_throughput_configurations[i] != -1){
                if ((optimal_index == -1 && arc_configurations[possible_throughput_configurations[i]].memory_overhead <= memory_constraint) || 
                    (arc_configurations[possible_throughput_configurations[i]].memory_overhead <= memory_constraint && arc_configurations[possible_throughput_configurations[i]].memory_overhead > optimal_overhead)){
                    optimal_index = possible_throughput_configurations[i];
                    optimal_overhead = arc_configurations[possible_throughput_configurations[i]].memory_overhead;
                }

                // Look for minimal approach as a backup
                if ((minimal_index == -1) || (arc_configurations[possible_throughput_configurations[i]].memory_overhead <= minimal_overhead)){
                    minimal_index = possible_throughput_configurations[i];
                    minimal_overhead = arc_configurations[possible_throughput_configurations[i]].memory_overhead;
                }
            }
        }
        // Free ECC Options array
        free(ecc_options);

        // Check to see if a valid configuration was found
        if (optimal_index == -1){
            printf("ARC Joint Optimization Warning: No such configuration satisfies provided memory constraint\n");
            printf("Using the lowest possible memory overhead approach that satisfies the throughput constraint. . .\n");
            // Set the parameters to the minimal index parameters
            *optimizer_parameter_a = arc_configurations[minimal_index].ecc_parameter_a;
            *optimizer_parameter_b = arc_configurations[minimal_index].ecc_parameter_b;
            *num_threads = arc_configurations[minimal_index].num_threads;
            return arc_configurations[minimal_index].ecc_algorithm;
        } else {
            // Set the parameters to the optimal index parameters
            *optimizer_parameter_a = arc_configurations[optimal_index].ecc_parameter_a;
            *optimizer_parameter_b = arc_configurations[optimal_index].ecc_parameter_b;
            *num_threads = arc_configurations[optimal_index].num_threads;
            return arc_configurations[optimal_index].ecc_algorithm;
        }
    // If there are some configurations that satisfy only the throughput constraint
    } else if (possible_throughput_config_index != 0 && possible_memory_config_index == 0){
        // Use the lowest memory overhead configuration that satisfies the throughput constraint
        int minimal_index = -1;
        double minimal_overhead = 0;

        for (i = 0; i < NUM_CONFIGURATIONS; i++){
            // Look for maximum bandwidth approach as a backup
            if (possible_throughput_configurations[i] != -1){
                // Look for maximum bandwidth approach as a backup
                if ((minimal_index == -1) || (arc_configurations[possible_throughput_configurations[i]].memory_overhead <= minimal_overhead)){
                    minimal_index = possible_throughput_configurations[i];
                    minimal_overhead = arc_configurations[possible_throughput_configurations[i]].memory_overhead;
                }
            }
        }
        // Free ECC Options array
        free(ecc_options);

        printf("ARC Joint Optimization Warning: No such configuration satisfies provided memory constraint\n");
        printf("Using the lowest possible memory overhead approach that satisfies the throughput constraint. . .\n");
        // Set the parameters to the minimal index parameters
        *optimizer_parameter_a = arc_configurations[minimal_index].ecc_parameter_a;
        *optimizer_parameter_b = arc_configurations[minimal_index].ecc_parameter_b;
        *num_threads = arc_configurations[minimal_index].num_threads;
        return arc_configurations[minimal_index].ecc_algorithm;

    // If there are some configurations that satisfy only the memory constraint
    } else if (possible_throughput_config_index == 0 && possible_memory_config_index != 0){
        // Use the highest throughput configuration that satisfies the memory constraint
        int maximum_index = -1;
        double maximum_overhead = 0;

        for (i = 0; i < NUM_CONFIGURATIONS; i++){
            // Look for maximum bandwidth approach as a backup
            if (possible_memory_configurations[i] != -1){
                // Look for maximum bandwidth approach as a backup
                if ((maximum_index == -1) || (arc_configurations[possible_memory_configurations[i]].throughput_overhead >= maximum_overhead)){
                    maximum_index = possible_memory_configurations[i];
                    maximum_overhead = arc_configurations[possible_memory_configurations[i]].throughput_overhead;
                }
            }
        }
        // Free ECC Options array
        free(ecc_options);

        printf("ARC Joint Optimization Warning: No such configuration satisfies provided throughput constraint\n");
        printf("Using the highest possible bandwidth approach that satisfies the memory constraint. . .\n");
        // Set the parameters to the maximum index parameters
        *optimizer_parameter_a = arc_configurations[maximum_index].ecc_parameter_a;
        *optimizer_parameter_b = arc_configurations[maximum_index].ecc_parameter_b;
        *num_threads = arc_configurations[maximum_index].num_threads;
        return arc_configurations[maximum_index].ecc_algorithm;

    // If there are no configurations that satisfy either constraint
    } else {
        // Use the lowest memory configuration with the maximum number of threads to maximize throughput
        int minimal_index = -1;
        double minimal_overhead = 0;

        for (i = 0; i < NUM_CONFIGURATIONS; i++){
            // Ensure ECC algorithm is viable given resiliency constraint
            if (ecc_options[(arc_configurations[i].ecc_algorithm - 1)] == 1){
                // Look for minimal memory overhead approach
                if ((minimal_index == -1) || (arc_configurations[i].memory_overhead <= minimal_overhead)){
                    minimal_index = i;
                    minimal_overhead = arc_configurations[i].memory_overhead;
                }
            }
        }
        // Free ECC Options array
        free(ecc_options);

        printf("ARC Joint Optimization Warning: No such configuration satisfies either provided constraint\n");
        printf("Using the lowest possible memory overhead approach with maximum number of threads possible. . .\n");
        // Set the parameters to the minimal index parameters with max number of threads
        *optimizer_parameter_a = arc_configurations[minimal_index].ecc_parameter_a;
        *optimizer_parameter_b = arc_configurations[minimal_index].ecc_parameter_b;
        *num_threads = AVAIL_THREADS;
        return arc_configurations[minimal_index].ecc_algorithm;
    }
}

// TODO: Make all decoding functions use OpenMP

// ARC Encoding Section
// #####################
// arc_calculate_parity_uint64:
// Calculates parity bit over 8 bytes of data
// params:
// byte         -   uint64_t 8 bytes to calculate secded on
// return:
// parity_bits  -   parity bit for 8 bytes of data
uint8_t arc_calculate_parity_uint64(uint64_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    uint8_t parity;
    int shift = 6; // log2(64 Bits) = 6
    int i;
    for (i = 0; i < shift; i++){
        byte = byte ^ (byte >> (1 << i));
    }
    parity = (uint8_t) byte & 1;
    return parity;
}

// arc_calculate_parity_uint8:
// Calculates parity bit over 1 byte of data
// params:
// byte         -   uint8_t byte to calculate secded on
// return:
// parity_bits  -   parity bit for byte of data
uint8_t arc_calculate_parity_uint8(uint8_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    uint8_t parity;
    int shift = 3; // log2(8 Bits) = 3
    int i;
    for (i = 0; i < shift; i++){
        byte = byte ^ (byte >> (1 << i));
    }
    parity = byte & 1;
    return parity;
}

// arc_parity_encode:
// Encodes each data block of given block size using single-bit parity
// params:
// data                 -   uint8_t data stream
// data_size            -   size of data stream
// block_size           -   number of bytes in each protected block
// threads              -   number of OpenMP threads to split work across
// encoded_data         -   address of pointer to encoded uint8_t protected data stream
// encoded_data_size    -   address of pointer to encoded size of protected data stream
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_parity_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Parity Encoding\n");
    }
    // Determine number of blocks given data_size and block_size
    uint32_t block_count = (uint32_t)ceil((float)data_size / (float)block_size);
    uint32_t remainder_size = data_size % block_size;

    // Calculate and create resulting array
    uint32_t metadata_length = 17;
    *encoded_data_size = data_size + block_count + metadata_length;
    *encoded_data = (uint8_t*)malloc(sizeof(uint8_t) * *encoded_data_size);

    // Write specific metadata back to array
    // Store Approach Identifier
    (*encoded_data)[0] = 0x01;
    // Store Block Count in 4 Bytes
    (*encoded_data)[1] = (uint8_t)((block_count & 0xFF000000) >> 24);  
    (*encoded_data)[2] = (uint8_t)((block_count & 0x00FF0000) >> 16);
    (*encoded_data)[3] = (uint8_t)((block_count & 0x0000FF00) >> 8);
    (*encoded_data)[4] = (uint8_t)((block_count & 0x000000FF));
    // Store Block Size in 4 Bytes
    (*encoded_data)[5] = (uint8_t)((block_size & 0xFF000000) >> 24);  
    (*encoded_data)[6] = (uint8_t)((block_size & 0x00FF0000) >> 16);
    (*encoded_data)[7] = (uint8_t)((block_size & 0x0000FF00) >> 8);
    (*encoded_data)[8] = (uint8_t)((block_size & 0x000000FF));
    // Store Remainder Block Size in 4 Bytes 
    (*encoded_data)[9] = (uint8_t)((remainder_size & 0xFF000000) >> 24); 
    (*encoded_data)[10] = (uint8_t)((remainder_size & 0x00FF0000) >> 16);
    (*encoded_data)[11] = (uint8_t)((remainder_size & 0x0000FF00) >> 8);
    (*encoded_data)[12] = (uint8_t)((remainder_size & 0x000000FF));
    // Store Number of Threads used
    (*encoded_data)[13] = (uint8_t)((threads & 0xFF000000) >> 24); 
    (*encoded_data)[14] = (uint8_t)((threads & 0x00FF0000) >> 16);
    (*encoded_data)[15] = (uint8_t)((threads & 0x0000FF00) >> 8);
    (*encoded_data)[16] = (uint8_t)((threads & 0x000000FF));
    // Set current index for encoded_data after metadata
    uint32_t encoded_data_index = metadata_length;
    
    // Setup Loop Variables 
    uint32_t blocks_processed;
    // Determine number of threads per block
    int n_per_thread;
    if (block_count < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = block_count / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Calculate parity for each block of data
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < block_count; blocks_processed++){
        // Declare private loop variables
        uint32_t current_block_size;

        // Determine loop offsets
        // Remainder Block
        if (blocks_processed == block_count-1 && remainder_size != 0){
            current_block_size = remainder_size;
        // Full Block
        } else {
            current_block_size = block_size;
        }
        // Loop Specific Index Locations
        uint32_t current_data_index = block_size*blocks_processed;
        uint32_t current_encoded_data_index = encoded_data_index + (block_size+1)*blocks_processed;

        // Setup Current Block Variables
        uint32_t block_data_start = current_data_index;
        uint8_t block_parity = 0;
        uint32_t bytes_processed = 0;
        int first_iter = 0;

        // Setup temporary variables
        uint64_t bytes_8;
        uint8_t bytes_1;
        uint8_t bytes_parity;

        // Iterate over all bytes and calculate overall block parity
        while(bytes_processed != current_block_size){
            // Work on 8 bytes at a time
            if (current_block_size - bytes_processed >= 8){
                // Load bytes
                bytes_8 = ((uint64_t)data[current_data_index] << 56) | ((uint64_t)data[current_data_index+1] << 48) | \
                ((uint64_t)data[current_data_index+2] << 40) | ((uint64_t)data[current_data_index+3] << 32) | \
                ((uint64_t)data[current_data_index+4] << 24) | ((uint64_t)data[current_data_index+5] << 16) | \
                ((uint64_t)data[current_data_index+6] << 8) | ((uint64_t)data[current_data_index+7]);
                current_data_index = current_data_index+8;

                // Calculate parity for 8 bytes
                bytes_parity = arc_calculate_parity_uint64(bytes_8);

                // Combine block parity and bytes parity
                if (first_iter){
                    block_parity = bytes_parity;
                    first_iter = 1;
                } else {
                    if (block_parity == 0 && bytes_parity == 0){
                        block_parity = 0;
                    } else if (block_parity == 1 && bytes_parity == 1) {
                        block_parity = 0;
                    } else {
                        block_parity = 1;
                    }
                }
                // Increment bytes processed by 8
                bytes_processed = bytes_processed + 8;

            // Work on a single byte at a time
            } else {
                // Load byte
                bytes_1 = data[current_data_index];
                current_data_index++;

                // Calculate parity for byte
                bytes_parity = arc_calculate_parity_uint8(bytes_1);

                // Combine block parity and byte parity
                if (first_iter){
                    block_parity = bytes_parity;
                    first_iter = 1;
                } else {
                    if (block_parity == 0 && bytes_parity == 0){
                        block_parity = 0;
                    } else if (block_parity == 1 && bytes_parity == 1) {
                        block_parity = 0;
                    } else {
                        block_parity = 1;
                    }
                }
                // Increment bytes processed by 1
                bytes_processed = bytes_processed + 1;
            }
        }

        // Upon block completion, write parity and data to resulting array
        (*encoded_data)[current_encoded_data_index] = block_parity;
        current_encoded_data_index++;
        int i;
        for (i = block_data_start; i < current_data_index; i++){
            (*encoded_data)[current_encoded_data_index] = data[i];
            current_encoded_data_index++;
        }
    }
    if (PRINT){
        printf("Parity Encoding Finished!\n");
    }
    // Return resulting array
    return 1;
}

// arc_parity_decode:
// Decodes data encoded using single-bit parity
// params:
// encoded_data         -   uint8_t protected data stream
// encoded_data_size    -   size of protected data stream
// data                 -   address of pointer to decoded uint8_t data stream
// data_size            -   address of pointer to decoded data stream size
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_parity_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Parity Decoding\n");
    }
    // Pull metadata from encoded_data array
    uint32_t metadata_length = 17;
    uint32_t block_count = (((uint32_t)encoded_data[1] & 0x000000FF) << 24) | (((uint32_t)encoded_data[2] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[3] & 0x000000FF) << 8) | (((uint32_t)encoded_data[4] & 0x000000FF));
    uint32_t block_size = (((uint32_t)encoded_data[5] & 0x000000FF) << 24) | (((uint32_t)encoded_data[6] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[7] & 0x000000FF) << 8) | (((uint32_t)encoded_data[8] & 0x000000FF));
    uint32_t remainder_size = (((uint32_t)encoded_data[9] & 0x000000FF) << 24) | (((uint32_t)encoded_data[10] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[11] & 0x000000FF) << 8) | (((uint32_t)encoded_data[12] & 0x000000FF));
    uint32_t threads = (((uint32_t)encoded_data[13] & 0x000000FF) << 24) | (((uint32_t)encoded_data[14] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[15] & 0x000000FF) << 8) | (((uint32_t)encoded_data[16] & 0x000000FF));
    uint32_t encoded_data_index = metadata_length;

    // Determine size of original array and allocate space
    *data_size = encoded_data_size - block_count - metadata_length;
    *data = (uint8_t*)malloc(sizeof(uint8_t) * *data_size);

    // Setup Loop Variables 
    uint32_t blocks_processed;
    int decode_success = 1;
    // Determine number of threads per block
    int n_per_thread;
    if (block_count < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = block_count / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Recalculate Parity for each block and compare to original parity
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < block_count; blocks_processed++){
        // Declare private loop variables
        uint32_t current_block_size; 

        // Determine loop offsets
        // Remainder Block
        if (blocks_processed == block_count-1 && remainder_size != 0){
            current_block_size = remainder_size;
        // Full Block
        } else {
            current_block_size = block_size;
        }
        // Loop Specific Index Locations
        uint32_t current_encoded_data_index = encoded_data_index + (block_size+1)*blocks_processed;
        uint32_t current_data_index = block_size*blocks_processed;

        // Pull original parity and save for later
        uint8_t original_block_parity = encoded_data[current_encoded_data_index];
        current_encoded_data_index++;

        // Setup Current Block Variables
        uint32_t block_data_start = current_encoded_data_index;
        uint8_t block_parity = 0;
        uint32_t bytes_processed = 0;
        int first_iter = 0;

        // Setup temporary variables
        uint64_t bytes_8;
        uint8_t bytes_1;
        uint8_t bytes_parity;

        // Iterate over all bytes and calculate overall block parity
        while(bytes_processed != current_block_size){
            // Work on 8 bytes at a time
            if (current_block_size - bytes_processed >= 8){
                // Load bytes
                bytes_8 = ((uint64_t)encoded_data[current_encoded_data_index] << 56) | ((uint64_t)encoded_data[current_encoded_data_index+1] << 48) | \
                ((uint64_t)encoded_data[current_encoded_data_index+2] << 40) | ((uint64_t)encoded_data[current_encoded_data_index+3] << 32) | \
                ((uint64_t)encoded_data[current_encoded_data_index+4] << 24) | ((uint64_t)encoded_data[current_encoded_data_index+5] << 16) | \
                ((uint64_t)encoded_data[current_encoded_data_index+6] << 8) | ((uint64_t)encoded_data[current_encoded_data_index+7]);
                current_encoded_data_index = current_encoded_data_index+8;

                // Calculate parity for 8 bytes
                bytes_parity = arc_calculate_parity_uint64(bytes_8);

                // Combine block parity and byte parity
                if (first_iter){
                    block_parity = bytes_parity;
                    first_iter = 1;
                } else {
                    if (block_parity == 0 && bytes_parity == 0){
                        block_parity = 0;
                    } else if (block_parity == 1 && bytes_parity == 1) {
                        block_parity = 0;
                    } else {
                        block_parity = 1;
                    }
                }
                // Increment bytes processed by 8
                bytes_processed = bytes_processed + 8;

            // Work on a single byte at a time
            } else {
                // Load byte
                bytes_1 = encoded_data[current_encoded_data_index];
                current_encoded_data_index++;

                // Calculate parity for byte
                bytes_parity = arc_calculate_parity_uint8(bytes_1);

                // Combine block parity and byte parity
                if (first_iter){
                    block_parity = bytes_parity;
                    first_iter = 1;
                } else {
                    if (block_parity == 0 && bytes_parity == 0){
                        block_parity = 0;
                    } else if (block_parity == 1 && bytes_parity == 1) {
                        block_parity = 0;
                    } else {
                        block_parity = 1;
                    }
                }
                // Increment bytes processed by 1
                bytes_processed = bytes_processed + 1;
            }
        }
        
        // Compare Original Block Parity and New Block Parity
        if (original_block_parity == block_parity){
            // If parity is correct, add data to array
            int i;
            for (i = block_data_start; i < current_encoded_data_index; i++){
                (*data)[current_data_index] = encoded_data[i];
                current_data_index++;
            }
        } else {
            // If parity is incorrect, print error and exit
            printf("DATA INTEGRITY ERROR: Single Bit Parity Error Found. . .\n");
            #pragma omp critical 
            {
                decode_success = 0;
            }
        }
    }
    
    // Free decoded data and set data_size to 0 if decoding process failed
    if (decode_success == 0){
        free(*data);
        *data_size = 0;
        if (PRINT){
            printf("Parity Decoding Failed!\n");
        }
    } else{
        if (PRINT){
            printf("Parity Decoding Finished!\n");
        }
    }
    // Return resulting array and decode success value
    return decode_success;
}

// arc_calculate_hamming_uint64:
// Calculate hamming over 8 bytes of data
// params:
// byte         -   uint64_t 8 bytes to calculate secded on
// return:
// parity_bits  -   parity bits for 8 bytes of data
uint8_t arc_calculate_hamming_uint64(uint64_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i;
    // Setup parity variable
    uint8_t parity_result = 0;    
    uint64_t current_parity;
    uint64_t parity_holder;

    // Calculate each parity bit
    for (i = 0; i < 7; i++){
        // Gather all 
        parity_holder = H_S_8_Parity_Matrix[i];
        // AND data with parity row
        parity_holder = parity_holder & byte;
        current_parity = 0;
        // XOR all resulting bits to get parity bit
        while (parity_holder) {
            current_parity ^= parity_holder & 1;
            parity_holder >>= 1;
        }
        // Store parity bit in parity_result
        parity_result = parity_result | ((uint8_t)current_parity << i);
    }
    
    // Return parity result
    return parity_result;
}

// arc_calculate_hamming_uint8:
// Calculates hamming over 1 byte of data
// params:
// byte         -   uint8_t byte to calculate secded on
// return:
// parity_bits  -   parity bits for byte of data
uint8_t arc_calculate_hamming_uint8(uint8_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i;
    // Setup parity variable
    uint8_t parity_result = 0;    
    uint8_t current_parity;
    uint8_t parity_holder;

    // Calculate each parity bit
    for (i = 0; i < 4; i++){
        // Gather all 
        parity_holder = H_S_1_Parity_Matrix[i];
        // AND data with parity row
        parity_holder = parity_holder & byte;
        current_parity = 0;
        // XOR all resulting bits to get parity bit
        while (parity_holder) {
            current_parity ^= parity_holder & 1;
            parity_holder >>= 1;
        }
        // Store parity bit in parity_result
        parity_result = parity_result | (current_parity << i);
    }
    
    // Return parity result
    return parity_result;
}


// arc_hamming_encode:
// Encodes each data block of given block size using hamming encoding
// params:
// data                 -   uint8_t data stream
// data_size            -   size of data stream
// block_size           -   number of bytes in each protected block (1 or 8)
// threads              -   number of OpenMP threads to split work across
// encoded_data         -   address of pointer to encoded uint8_t protected data stream
// encoded_data_size    -   address of pointer to encoded size of protected data stream
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_hamming_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Hamming Encoding\n");
    }
    // Check that block_size is valid and supported
    if (block_size != 1 && block_size != 8){
        printf("INVALID BLOCK SIZE: Please select a supported block size. . .\n");
        return 0;
    }
    // Set metadata length
    uint32_t metadata_length = 17;

    // Determine number of blocks given data_size and block_size
    uint32_t block_count = (uint32_t)floor((float)data_size / (float)block_size);
    uint32_t remainder_blocks = data_size % block_size;

    // Given block size, determine number of parity bytes for each block
    if (block_size == 1 || block_size == 8){
        *encoded_data_size = data_size + block_count + remainder_blocks + metadata_length;
    }
    // Create encoded_data array from resulting 
    *encoded_data = (uint8_t*)malloc(sizeof(uint8_t) * *encoded_data_size);

    // Write specific metadata back to array
    // Store Approach Identifier
    (*encoded_data)[0] = 0x02;
    // Store Block Count in 4 Bytes 
    (*encoded_data)[1] = (uint8_t)((block_count & 0xFF000000) >> 24);  
    (*encoded_data)[2] = (uint8_t)((block_count & 0x00FF0000) >> 16);
    (*encoded_data)[3] = (uint8_t)((block_count & 0x0000FF00) >> 8);
    (*encoded_data)[4] = (uint8_t)((block_count & 0x000000FF));
    // Store Block Size in 4 Bytes 
    (*encoded_data)[5] = (uint8_t)((block_size & 0xFF000000) >> 24);  
    (*encoded_data)[6] = (uint8_t)((block_size & 0x00FF0000) >> 16);
    (*encoded_data)[7] = (uint8_t)((block_size & 0x0000FF00) >> 8);
    (*encoded_data)[8] = (uint8_t)((block_size & 0x000000FF));    
    // Store Remainder Blocks in 4 Bytes 
    (*encoded_data)[9] = (uint8_t)((remainder_blocks & 0xFF000000) >> 24);  
    (*encoded_data)[10] = (uint8_t)((remainder_blocks & 0x00FF0000) >> 16);
    (*encoded_data)[11] = (uint8_t)((remainder_blocks & 0x0000FF00) >> 8);
    (*encoded_data)[12] = (uint8_t)((remainder_blocks & 0x000000FF));
    // Store Number of Threads used
    (*encoded_data)[13] = (uint8_t)((threads & 0xFF000000) >> 24); 
    (*encoded_data)[14] = (uint8_t)((threads & 0x00FF0000) >> 16);
    (*encoded_data)[15] = (uint8_t)((threads & 0x0000FF00) >> 8);
    (*encoded_data)[16] = (uint8_t)((threads & 0x000000FF));
    // Set current index for encoded_data after metadata
    uint32_t encoded_data_index = metadata_length;

    // Setup Loop Variables 
    uint32_t blocks_processed;
    // Determine number of threads per block
    int n_per_thread;
    if ((block_count+remainder_blocks) < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = (block_count+remainder_blocks) / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Calculate hamming for each block of data
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < (block_count + remainder_blocks); blocks_processed++){
        // Declare private loop variables
        uint32_t current_data_index;
        uint32_t current_encoded_data_index;
        uint32_t current_block_size;

        // Determine loop offsets
        // Full Blocks
        if (blocks_processed < block_count){
            current_data_index = blocks_processed * block_size;
            current_encoded_data_index = encoded_data_index + (blocks_processed * (block_size + 1));
            current_block_size = block_size;
        // Remainder Blocks
        } else {
            current_data_index = block_count * block_size + (blocks_processed - block_count);
            current_encoded_data_index = (encoded_data_index + (block_count * (block_size + 1))) + ((blocks_processed - block_count) * 2);
            current_block_size = 1;
        }

        // Setup temporary variables
        uint8_t block_1;
        uint64_t block_8;
        uint8_t block_parity;

        // Gather data based on current_block_size
        if (current_block_size == 8){
            // Load data
            block_8 = ((uint64_t)data[current_data_index] << 56) | ((uint64_t)data[current_data_index+1] << 48) | \
            ((uint64_t)data[current_data_index+2] << 40) | ((uint64_t)data[current_data_index+3] << 32) | \
            ((uint64_t)data[current_data_index+4] << 24) | ((uint64_t)data[current_data_index+5] << 16) | \
            ((uint64_t)data[current_data_index+6] << 8) | ((uint64_t)data[current_data_index+7]);

            // Calculate Hamming on 8 Byte block
            block_parity = arc_calculate_hamming_uint64(block_8);

            // Upon completion, wriet parity bits and 8 bytes of data to resulting array
            (*encoded_data)[current_encoded_data_index] = block_parity;
            current_encoded_data_index++;
            int i;
            for (i = current_data_index; i < current_data_index+block_size; i++){
                (*encoded_data)[current_encoded_data_index] = data[i];
                current_encoded_data_index++;
            }

        } else if (current_block_size == 1){
            // Load data
            block_1 = data[current_data_index];

            // Calculate Hamming on 1 Byte block
            block_parity = arc_calculate_hamming_uint8(block_1);

            // Upon completion, write parity bits and 1 byte of data to resulting array
            (*encoded_data)[current_encoded_data_index] = block_parity;
            current_encoded_data_index++;
            (*encoded_data)[current_encoded_data_index] = data[current_data_index];
            current_encoded_data_index++;
        }
    }
    
    // Return resulting array
    if (PRINT){
        printf("Hamming Encoding Finished\n");
    }
    return 1;
}

// arc_hamming_decode:
// Decodes data encoded using hamming encoding
// params:
// encoded_data         -   uint8_t protected data stream
// encoded_data_size    -   size of protected data stream
// data                 -   address of pointer to decoded uint8_t data stream
// data_size            -   address of pointer to decoded data stream size
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_hamming_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Hamming Decoding\n");
    }
    // Set metadata length
    uint32_t metadata_length = 17;

    // Pull metadata from encoded_data array
    uint32_t block_count = (((uint32_t)encoded_data[1] & 0x000000FF) << 24) | (((uint32_t)encoded_data[2] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[3] & 0x000000FF) << 8) | (((uint32_t)encoded_data[4] & 0x000000FF));
    uint32_t block_size = (((uint32_t)encoded_data[5] & 0x000000FF) << 24) | (((uint32_t)encoded_data[6] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[7] & 0x000000FF) << 8) | (((uint32_t)encoded_data[8] & 0x000000FF));
    uint32_t remainder_blocks = (((uint32_t)encoded_data[9] & 0x000000FF) << 24) | (((uint32_t)encoded_data[10] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[11] & 0x000000FF) << 8) | (((uint32_t)encoded_data[12] & 0x000000FF));
    uint32_t threads = (((uint32_t)encoded_data[13] & 0x000000FF) << 24) | (((uint32_t)encoded_data[14] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[15] & 0x000000FF) << 8) | (((uint32_t)encoded_data[16] & 0x000000FF));
    uint32_t encoded_data_index = metadata_length;

    // Check that block_size is valid and supported
    if (block_size != 1 && block_size != 8){
        printf("INVALID BLOCK SIZE: Please select a supported block size. . .\n");
        return 0;
    } 

    // Determine size of original array and allocate it
    if (block_size == 1 || block_size == 8){
        *data_size = encoded_data_size - block_count - remainder_blocks - metadata_length;
    }
    *data = (uint8_t*)malloc(sizeof(uint8_t) * *data_size);

    // Setup Loop Variables 
    uint32_t blocks_processed;
    int decode_success = 1;
    // Determine number of threads per block
    int n_per_thread;
    if ((block_count+remainder_blocks) < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = (block_count+remainder_blocks) / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Recalculate hamming for each block of data
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < (block_count + remainder_blocks); blocks_processed++){
        // Declare private loop variables
        uint32_t current_encoded_data_index;
        uint32_t current_data_index;
        uint32_t current_block_size;

        // Determine loop offsets
        // Full Blocks
        if (blocks_processed < block_count){
            current_encoded_data_index = encoded_data_index + (blocks_processed * (block_size + 1));
            current_data_index = blocks_processed * block_size;
            current_block_size = block_size;
        // Remainder Blocks
        } else {
            current_encoded_data_index = (encoded_data_index + (block_count * (block_size + 1))) + ((blocks_processed - block_count) * 2);
            current_data_index = block_count * block_size + (blocks_processed - block_count);
            current_block_size = 1;
        }

        // Pull original parity bits and save for later
        uint8_t original_block_parity;
        if (current_block_size == 1 || current_block_size == 8){
            original_block_parity = encoded_data[current_encoded_data_index];
            current_encoded_data_index++;
        }

        // Set current block variables
        uint32_t block_data_start = current_encoded_data_index;

        // Setup temporary variables
        uint8_t block_1;
        uint64_t block_8;
        uint8_t block_parity;
        uint8_t xor_parity;

        // Gather data for current block 
        if (current_block_size == 8){
            // Load data 
            block_8 =  ((uint64_t)encoded_data[current_encoded_data_index] << 56) | ((uint64_t)encoded_data[current_encoded_data_index+1] << 48) | \
            ((uint64_t)encoded_data[current_encoded_data_index+2] << 40) | ((uint64_t)encoded_data[current_encoded_data_index+3] << 32) | \
            ((uint64_t)encoded_data[current_encoded_data_index+4] << 24) | ((uint64_t)encoded_data[current_encoded_data_index+5] << 16) | \
            ((uint64_t)encoded_data[current_encoded_data_index+6] << 8) | ((uint64_t)encoded_data[current_encoded_data_index+7]);
            current_encoded_data_index = current_encoded_data_index+8;

            // Calculate Hamming on 8 Byte block
            block_parity = arc_calculate_hamming_uint64(block_8);

            // Upon completion, compare both parity bits to ensure data integrity
            xor_parity = original_block_parity ^ block_parity;

            if (xor_parity == 0){
               // If parity is correct, add data to array
               int i;
               for (i = block_data_start; i < current_encoded_data_index; i++){
                   (*data)[current_data_index] = encoded_data[i];
                   current_data_index++;
               }
            } else {
                printf("Incorrect Parity Found...\nAttempting to Fix Now...\n");
                // Using xor_parity value attempt to fix data and parity
                int j;
                for (j = 0; j < 72; j++){
                    if (xor_parity == H_8_Syndrome_Table[j]){
                        if (j < 64){
                            block_8 = block_8 ^ ((uint64_t)1 << j);
                            block_parity = arc_calculate_hamming_uint64(block_8);
                        } else {
                            block_parity = block_parity ^ (1 << (j-64));
                        }
                        break;
                    }
                }

                // Recalculate Hamming on "fixed" data and parity
                xor_parity = original_block_parity ^ block_parity;
                if (xor_parity == 0){
                    printf("Data Corrected!\n");
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0xFF00000000000000) >> 56);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x00FF000000000000) >> 48);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x0000FF0000000000) >> 40);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x000000FF00000000) >> 32);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x00000000FF000000) >> 24);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x0000000000FF0000) >> 16);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x000000000000FF00) >> 8);
                    current_data_index++;
                    (*data)[current_data_index] = (uint8_t)((block_8 & 0x00000000000000FF));
                    current_data_index++;
                } else {
                    printf("DATA INTEGRITY ERROR: Parity Error Found, Correction Failed...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                }
            }
        } else if (current_block_size == 1){
            // Load data 
            block_1 = encoded_data[current_encoded_data_index];
            current_encoded_data_index++;

            // Calculate Hamming on 1 Byte block
            block_parity = arc_calculate_hamming_uint8(block_1);

            // Upon completion, compare both parity bits to ensure data integrity
            xor_parity = original_block_parity ^ block_parity;

            if (xor_parity == 0){
               // If parity is correct, add data to array
               (*data)[current_data_index] = encoded_data[block_data_start];
               current_data_index++;
            } else {
                printf("Incorrect Parity Found...\nAttempting to Fix Now...\n");
                // Using xor_parity value attempt to fix data and parity
                int j;
                for (j = 0; j < 16; j++){
                    if (xor_parity == H_1_Syndrome_Table[j]){
                        if (j < 8){
                            block_1 = block_1 ^ ( 1 << j);
                            block_parity = arc_calculate_hamming_uint8(block_1);
                        } else {
                            block_parity = block_parity ^ ( 1 << (j-8));
                        }
                        break;
                    }
                }

                // Recalculate Hamming on "fixed" data and parity
                xor_parity = original_block_parity ^ block_parity;
                if (xor_parity == 0){
                    printf("Data Corrected!\n");
                    (*data)[current_data_index] = block_1;
                    current_data_index++;
                } else {
                    printf("DATA INTEGRITY ERROR: Parity Error Found, Correction Failed...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                }
            }
        }
    }

    // Free decoded data and set data_size to 0 if decoding process failed
    if (decode_success == 0){
        free(*data);
        *data_size = 0;
        if (PRINT){
            printf("Hamming Decoding Failed!\n");
        }
    } else{
        if (PRINT){
            printf("Hamming Decoding Finished!\n");
        }
    }
    // Return resulting array and decode success value
    return decode_success;
}

// arc_calculate_secded_uint64:
// Calculates secded over 8 bytes of data
// params:
// byte         -   uint64_t 8 bytes to calculate secded on
// return:
// parity_bits  -   parity bits for 8 bytes of data
uint8_t arc_calculate_secded_uint64(uint64_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i;
    // Setup parity variable
    uint8_t parity_result = 0;    
    uint64_t current_parity;
    uint64_t parity_holder;

    // Calculate each parity bit
    for (i = 0; i < 7; i++){
        // Gather all 
        parity_holder = H_S_8_Parity_Matrix[i];
        // AND data with parity row
        parity_holder = parity_holder & byte;
        current_parity = 0;
        // XOR all resulting bits to get parity bit
        while (parity_holder) {
            current_parity ^= parity_holder & 1;
            parity_holder >>= 1;
        }
        // Store parity bit in parity_result
        parity_result = parity_result | ((uint8_t)current_parity << i);
    }

    // Calculate overall parity of parity and data
    uint8_t parity_parity = arc_calculate_parity_uint8(parity_result);
    uint8_t byte_parity = arc_calculate_parity_uint64(byte);

    // Combine all parity to find final parity bit
    if (parity_parity != byte_parity){
        parity_result = parity_result | (1 << 7);
    }
    
    // Return parity result
    return parity_result;
}

// arc_calculate_secded_uint8:
// Calculates secded over 1 byte of data
// params:
// byte         -   uint8_t byte to calculate secded on
// return:
// parity_bits  -   parity bits for byte of data
uint8_t arc_calculate_secded_uint8(uint8_t byte){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    int i;
    // Setup parity variable
    uint8_t parity_result = 0;    
    uint8_t current_parity;
    uint8_t parity_holder;

    // Calculate each parity bit
    for (i = 0; i < 4; i++){
        // Gather all 
        parity_holder = H_S_1_Parity_Matrix[i];
        // AND data with parity row
        parity_holder = parity_holder & byte;
        current_parity = 0;
        // XOR all resulting bits to get parity bit
        while (parity_holder) {
            current_parity ^= parity_holder & 1;
            parity_holder >>= 1;
        }
        // Store parity bit in parity_result
        parity_result = parity_result | (current_parity << i);
    }

    // Calculate overall parity of parity and data
    uint8_t parity_parity = arc_calculate_parity_uint8(parity_result);
    uint8_t byte_parity = arc_calculate_parity_uint8(byte);

    // Combine all parity to find final parity bit
    if (parity_parity != byte_parity){
        parity_result = parity_result | (1 << 4);
    }
    
    // Return parity result
    return parity_result;
}

// arc_secded_encode:
// Encodes each data block of given block size using secded encoding
// params:
// data                 -   uint8_t data stream
// data_size            -   size of data stream
// block_size           -   number of bytes in each protected block
// threads              -   number of OpenMP threads to split work across
// encoded_data         -   address of pointer to encoded uint8_t protected data stream
// encoded_data_size    -   address of pointer to encoded size of protected data stream
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_secded_encode(uint8_t* data, uint32_t data_size, uint32_t block_size, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting SECDED Encoding\n");
    }

    // Check that block_size is valid and supported
    if (block_size != 1 && block_size != 8){
        printf("INVALID BLOCK SIZE: Please select a supported block size. . .\n");
        return 0;
    }
    // Set metadata length
    uint32_t metadata_length = 17;

    // Determine number of blocks given data_size and block_size
    uint32_t block_count = (uint32_t)floor((float)data_size / (float)block_size);
    uint32_t remainder_blocks = data_size % block_size;

    // Given block size, determine number of parity bytes for each block
    if (block_size == 1 || block_size == 8){
        *encoded_data_size = data_size + block_count + remainder_blocks + metadata_length;
    }
    // Create encoded_data array from resulting 
    *encoded_data = (uint8_t*)malloc(sizeof(uint8_t) * *encoded_data_size);

    // Write specific metadata back to array
    // Store Approach Identifier
    (*encoded_data)[0] = 0x03;
    // Store Block Count in 4 Bytes 
    (*encoded_data)[1] = (uint8_t)((block_count & 0xFF000000) >> 24);  
    (*encoded_data)[2] = (uint8_t)((block_count & 0x00FF0000) >> 16);
    (*encoded_data)[3] = (uint8_t)((block_count & 0x0000FF00) >> 8);
    (*encoded_data)[4] = (uint8_t)((block_count & 0x000000FF));
    // Store Block Size in 4 Bytes 
    (*encoded_data)[5] = (uint8_t)((block_size & 0xFF000000) >> 24);  
    (*encoded_data)[6] = (uint8_t)((block_size & 0x00FF0000) >> 16);
    (*encoded_data)[7] = (uint8_t)((block_size & 0x0000FF00) >> 8);
    (*encoded_data)[8] = (uint8_t)((block_size & 0x000000FF));    
    // Store Remainder Blocks in 4 Bytes 
    (*encoded_data)[9] = (uint8_t)((remainder_blocks & 0xFF000000) >> 24);  
    (*encoded_data)[10] = (uint8_t)((remainder_blocks & 0x00FF0000) >> 16);
    (*encoded_data)[11] = (uint8_t)((remainder_blocks & 0x0000FF00) >> 8);
    (*encoded_data)[12] = (uint8_t)((remainder_blocks & 0x000000FF));
    // Store Number of Threads used
    (*encoded_data)[13] = (uint8_t)((threads & 0xFF000000) >> 24); 
    (*encoded_data)[14] = (uint8_t)((threads & 0x00FF0000) >> 16);
    (*encoded_data)[15] = (uint8_t)((threads & 0x0000FF00) >> 8);
    (*encoded_data)[16] = (uint8_t)((threads & 0x000000FF));
    // Set current index for encoded_data after metadata
    uint32_t encoded_data_index = metadata_length;

     // Setup Loop Variables 
    uint32_t blocks_processed;
    // Determine number of threads per block
    int n_per_thread;
    if ((block_count+remainder_blocks) < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = (block_count+remainder_blocks) / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Calculate SECDED for each block of data
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < (block_count + remainder_blocks); blocks_processed++){
        // Declare private loop variables
        uint32_t current_data_index;
        uint32_t current_encoded_data_index;
        uint32_t current_block_size;

        // Determine loop offsets
        // Full Blocks
        if (blocks_processed < block_count){
            current_data_index = blocks_processed * block_size;
            current_encoded_data_index = encoded_data_index + (blocks_processed * (block_size + 1));
            current_block_size = block_size;
        // Remainder Blocks
        } else {
            current_data_index = block_count * block_size + (blocks_processed - block_count);
            current_encoded_data_index = (encoded_data_index + (block_count * (block_size + 1))) + ((blocks_processed - block_count) * 2);
            current_block_size = 1;
        }

        // Setup temporary variables
        uint8_t block_1;
        uint64_t block_8;
        uint8_t block_parity;

        // Gather data based on current_block_size
        if (current_block_size == 8){
            // Load data
            block_8 = ((uint64_t)data[current_data_index] << 56) | ((uint64_t)data[current_data_index+1] << 48) | \
            ((uint64_t)data[current_data_index+2] << 40) | ((uint64_t)data[current_data_index+3] << 32) | \
            ((uint64_t)data[current_data_index+4] << 24) | ((uint64_t)data[current_data_index+5] << 16) | \
            ((uint64_t)data[current_data_index+6] << 8) | ((uint64_t)data[current_data_index+7]);

            // Calculate SECDED on 8 Byte block
            block_parity = arc_calculate_secded_uint64(block_8);

            // Upon completion, wriet parity bits and 8 bytes of data to resulting array
            (*encoded_data)[current_encoded_data_index] = block_parity;
            current_encoded_data_index++;
            int i;
            for (i = current_data_index; i < current_data_index+block_size; i++){
                (*encoded_data)[current_encoded_data_index] = data[i];
                current_encoded_data_index++;
            }

        } else if (current_block_size == 1){
            // Load data
            block_1 = data[current_data_index];

            // Calculate SECDED on 1 Byte block
            block_parity = arc_calculate_secded_uint8(block_1);

            // Upon completion, write parity bits and 1 byte of data to resulting array
            (*encoded_data)[current_encoded_data_index] = block_parity;
            current_encoded_data_index++;
            (*encoded_data)[current_encoded_data_index] = data[current_data_index];
            current_encoded_data_index++;
        }
    }
    
    // Return resulting array
    if (PRINT){
        printf("SECDED Encoding Finished\n");
    }
    return 1;
}

// arc_secded_decode:
// Decodes data encoded using secded encoding
// params:
// encoded_data         -   uint8_t protected data stream
// encoded_data_size    -   size of protected data stream
// data                 -   address of pointer to decoded uint8_t data stream
// data_size            -   address of pointer to decoded data stream size
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_secded_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t* data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting SECDED Decoding\n");
    }
    // Set metadata length
    uint32_t metadata_length = 17;

    // Pull metadata from encoded_data array
    uint32_t block_count = (((uint32_t)encoded_data[1] & 0x000000FF) << 24) | (((uint32_t)encoded_data[2] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[3] & 0x000000FF) << 8) | (((uint32_t)encoded_data[4] & 0x000000FF));
    uint32_t block_size = (((uint32_t)encoded_data[5] & 0x000000FF) << 24) | (((uint32_t)encoded_data[6] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[7] & 0x000000FF) << 8) | (((uint32_t)encoded_data[8] & 0x000000FF));
    uint32_t remainder_blocks = (((uint32_t)encoded_data[9] & 0x000000FF) << 24) | (((uint32_t)encoded_data[10] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[11] & 0x000000FF) << 8) | (((uint32_t)encoded_data[12] & 0x000000FF));
    uint32_t threads = (((uint32_t)encoded_data[13] & 0x000000FF) << 24) | (((uint32_t)encoded_data[14] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[15] & 0x000000FF) << 8) | (((uint32_t)encoded_data[16] & 0x000000FF));
    uint32_t encoded_data_index = 17;

    // Check that block_size is valid and supported
    if (block_size != 1 && block_size != 8){
        printf("INVALID BLOCK SIZE: Please select a supported block size. . .\n");
        return 0;
    } 

    // Determine size of original array and allocate it
    if (block_size == 1 || block_size == 8){
        *data_size = encoded_data_size - block_count - remainder_blocks - 17;
    }
    *data = (uint8_t*)malloc(sizeof(uint8_t) * *data_size);

    // Setup Loop Variables 
    uint32_t blocks_processed;
    int decode_success = 1;
    // Determine number of threads per block
    int n_per_thread;
    if ((block_count+remainder_blocks) < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = (block_count+remainder_blocks) / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);

    // Recalculate secded for each block of data
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < (block_count + remainder_blocks); blocks_processed++){
        // Declare private loop variables
        uint32_t current_encoded_data_index;
        uint32_t current_data_index;
        uint32_t current_block_size;

        // Determine loop offsets
        // Full Blocks
        if (blocks_processed < block_count){
            current_encoded_data_index = encoded_data_index + (blocks_processed * (block_size + 1));
            current_data_index = blocks_processed * block_size;
            current_block_size = block_size;
        // Remainder Blocks
        } else {
            current_encoded_data_index = (encoded_data_index + (block_count * (block_size + 1))) + ((blocks_processed - block_count) * 2);
            current_data_index = block_count * block_size + (blocks_processed - block_count);
            current_block_size = 1;
        }

        // Pull original parity bits and save for later
        uint8_t original_block_parity;
        if (current_block_size == 1 || current_block_size == 8){
            original_block_parity = encoded_data[current_encoded_data_index];
            current_encoded_data_index++;
        }

        // Set current block variables
        uint32_t block_data_start = current_encoded_data_index;

        // Setup temporary variables
        uint8_t block_1;
        uint64_t block_8;
        uint8_t block_parity;
        uint8_t xor_parity;
        uint8_t secded_parity_bit;
        uint8_t temp_original_block_parity;
        uint8_t block_parity_parity;
        uint8_t block_data_parity;
        uint8_t other_secded_parity_bit;

        // Gather data for current block 
        if (current_block_size == 8){
            // Load data 
            block_8 =  ((uint64_t)encoded_data[current_encoded_data_index] << 56) | ((uint64_t)encoded_data[current_encoded_data_index+1] << 48) | \
            ((uint64_t)encoded_data[current_encoded_data_index+2] << 40) | ((uint64_t)encoded_data[current_encoded_data_index+3] << 32) | \
            ((uint64_t)encoded_data[current_encoded_data_index+4] << 24) | ((uint64_t)encoded_data[current_encoded_data_index+5] << 16) | \
            ((uint64_t)encoded_data[current_encoded_data_index+6] << 8) | ((uint64_t)encoded_data[current_encoded_data_index+7]);
            current_encoded_data_index = current_encoded_data_index+8;

            // Calculate secded on 8 Byte block
            block_parity = arc_calculate_secded_uint64(block_8);

            // Upon completion, compare both parity bits to ensure data integrity
            xor_parity = original_block_parity ^ block_parity;

            if (xor_parity == 0){
               // If parity is correct, add data to array
               int i;
               for (i = block_data_start; i < current_encoded_data_index; i++){
                   (*data)[current_data_index] = encoded_data[i];
                   current_data_index++;
               }
            } else {
                printf("Incorrect Parity Found...\nChecking for Double Bit Errors...\n");
                // With syndrome being non-zero, check parity to see if double bit error occured
                // Grab current secded parity bit
                secded_parity_bit = ((original_block_parity & 0b10000000) >> 7);
                // Remove secded parity bit from block parity
                temp_original_block_parity = (original_block_parity & 0b01111111);

                // Calculate actual parity on faulty codeword data
                block_parity_parity = arc_calculate_parity_uint8(temp_original_block_parity);
                block_data_parity = arc_calculate_parity_uint64(block_8);
                // Combine found parity
                other_secded_parity_bit = 0;
                if (block_parity_parity != block_data_parity){
                    other_secded_parity_bit = 1;
                }
                // If secded parity is the same as previous block secded parity, then Double Error Detected
                if (secded_parity_bit == other_secded_parity_bit){
                    printf("DATA INTEGRITY ERROR: Double Bit Error Found, Correction Impossible...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                // If parity is not the same then a single bit error can be attempted to be fixed
                } else {
                    printf("Double Bit Error Not Found...\nAttempting to Fix...\n");
                    // Using xor_parity value attempt to fix data and parity
                    int j;
                    for (j = 0; j < 72; j++){
                        if (xor_parity == S_8_Syndrome_Table[j]){
                            if (j < 64){
                                block_8 = block_8 ^ ((uint64_t)1 << j);
                                block_parity = arc_calculate_secded_uint64(block_8);
                            } else {
                                block_parity = block_parity ^ (1 << (j-64));
                            }
                            break;
                        }
                    }

                    // Recalculate secded on "fixed" data and parity
                    xor_parity = original_block_parity ^ block_parity;
                    if (xor_parity == 0){
                        printf("Data Corrected!\n");
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0xFF00000000000000) >> 56);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x00FF000000000000) >> 48);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x0000FF0000000000) >> 40);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x000000FF00000000) >> 32);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x00000000FF000000) >> 24);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x0000000000FF0000) >> 16);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x000000000000FF00) >> 8);
                        current_data_index++;
                        (*data)[current_data_index] = (uint8_t)((block_8 & 0x00000000000000FF));
                        current_data_index++;
                    } else {
                        printf("DATA INTEGRITY ERROR: Parity Error Found, Correction Failed...\n");
                        #pragma omp critical 
                        {
                            decode_success = 0;
                        }
                    }
                }
            }
        } else if (current_block_size == 1){
            // Load data 
            block_1 = encoded_data[current_encoded_data_index];
            current_encoded_data_index++;

            // Calculate secded on 1 Byte block
            block_parity = arc_calculate_secded_uint8(block_1);

            // Upon completion, compare both parity bits to ensure data integrity
            xor_parity = original_block_parity ^ block_parity;

            if (xor_parity == 0){
               // If parity is correct, add data to array
               (*data)[current_data_index] = encoded_data[block_data_start];
               current_data_index++;
            } else {
                printf("Incorrect Parity Found...\nChecking for Double Bit Errors...\n");
                // With syndrome being non-zero, check parity to see if double bit error occured
                // Grab current secded parity bit
                secded_parity_bit = ((original_block_parity & 0b00010000) >> 4);
                // Remove secded parity bit from block parity
                temp_original_block_parity = (original_block_parity & 0b11101111);

                // Calculate actual parity on faulty codeword data
                block_parity_parity = arc_calculate_parity_uint8(temp_original_block_parity);
                block_data_parity = arc_calculate_parity_uint8(block_1);
                // Combine found parity
                other_secded_parity_bit = 0;
                if (block_parity_parity != block_data_parity){
                    other_secded_parity_bit = 1;
                }
                // If secded parity is the same as previous block secded parity, then Double Error Detected
                if (secded_parity_bit == other_secded_parity_bit){
                    printf("DATA INTEGRITY ERROR: Double Bit Error Found, Correction Impossible...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                // If parity is not the same then a single bit error can be attempted to be fixed
                } else {
                    printf("Double Bit Error Not Found...\nAttempting to Fix...\n");
                    // Using xor_parity value attempt to fix data and parity
                    int j;
                    for (j = 0; j < 16; j++){
                        if (xor_parity == S_1_Syndrome_Table[j]){
                            if (j < 8){
                                block_1 = block_1 ^ ( 1 << j);
                                block_parity = arc_calculate_secded_uint8(block_1);
                            } else {
                                block_parity = block_parity ^ ( 1 << (j-8));
                            }
                            break;
                        }
                    }

                    // Recalculate secded on "fixed" data and parity
                    xor_parity = original_block_parity ^ block_parity;
                    if (xor_parity == 0){
                        printf("Data Corrected!\n");
                        (*data)[current_data_index] = block_1;
                        current_data_index++;
                    } else {
                        printf("DATA INTEGRITY ERROR: Unknown Bit Error Found, Correction Failed...\n");
                        #pragma omp critical 
                        {
                            decode_success = 0;
                        }
                    }
                }
            }
        }
    }

    // Free decoded data and set data_size to 0 if decoding process failed
    if (decode_success == 0){
        free(*data);
        *data_size = 0;
        if (PRINT){
            printf("SECDED Decoding Failed!\n");
        }
    } else{
        if (PRINT){
            printf("SECDED Decoding Finished!\n");
        }
    }
    // Return resulting array and decode success value
    return decode_success;
}

// arc_reed_solomon_encode:
// Encodes each data block of (8*data_devices) bytes using (8*code_devices) with reed-solomon encoding
// params:
// data                 -   uint8_t data stream
// data_size            -   size of data stream
// data_devices         -   number of data devices per block
// code_devices         -   number of code devices per block
// threads              -   number of OpenMP threads to split work across
// encoded_data         -   address of pointer to encoded uint8_t protected data stream
// encoded_data_size    -   address of pointer to encoded size of protected data stream
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
// note:
// data_devices + code_devices <= 2^8
int arc_reed_solomon_encode(uint8_t* data, uint32_t data_size, uint32_t data_devices, uint32_t code_devices, uint32_t threads, uint8_t** encoded_data, uint32_t* encoded_data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Reed Solomon Encoding\n");
    }
    // Set metadata length 
    uint32_t metadata_length = 25;

    // Ensure k+m <= 2^w
    if (data_devices + code_devices > ((long)1 << 8)){
        printf("INVALID PARAMETER SIZES: Please select inputs where data devices + code devices <= 2^8. . .\n");
        return 0;
    }

    // Determine number of blocks given data_size and data_devices
    // Calculate number of longs within the overall data
    uint32_t longs_per_data = (uint32_t)ceil((float)data_size / (float)8);
    // Calculate remainder bytes in final long of data
    uint32_t remainder_long_in_data = data_size % 8;
    // Calculate number of ecc blocks that are needed given longs per data
    uint32_t block_count = (uint32_t)ceil((float)longs_per_data / (float)data_devices);
    // Calculate number of data devices in final ecc block
    uint32_t remainder_data_devices = longs_per_data % data_devices;

    // Create encoded_data array given sizes
    uint32_t data_bytes;
    uint32_t code_bytes;
    uint32_t parity_bytes;

    // Calculate bytes needed for data
    if (remainder_data_devices == 0){
        data_bytes = data_devices * 8 * block_count;
    } else {
        data_bytes = (data_devices * 8 * (block_count-1)) + (remainder_data_devices * 8);
    }
    // Calculate bytes needed for coding
    code_bytes = code_devices * 8 * block_count;
    // Calculate bytes needed for parity
    if (remainder_data_devices == 0){
        parity_bytes = (data_devices + code_devices) * block_count;
    } else {
        parity_bytes = ((data_devices + code_devices) * (block_count-1)) + (code_devices + remainder_data_devices);
    }

    // Create encoded_data array from resulting 
    *encoded_data_size = metadata_length + data_bytes + code_bytes + parity_bytes;
    *encoded_data = (uint8_t*)malloc(sizeof(uint8_t) * *encoded_data_size);

    // Write specific metadata back to array
    // Store Approach Identifier
    (*encoded_data)[0] = 0x04;
    // Store Data Devices
    (*encoded_data)[1] = (uint8_t)((data_devices & 0xFF000000) >> 24);  
    (*encoded_data)[2] = (uint8_t)((data_devices & 0x00FF0000) >> 16);
    (*encoded_data)[3] = (uint8_t)((data_devices & 0x0000FF00) >> 8);
    (*encoded_data)[4] = (uint8_t)((data_devices & 0x000000FF));
    // Store Code Devices
    (*encoded_data)[5] = (uint8_t)((code_devices & 0xFF000000) >> 24);  
    (*encoded_data)[6] = (uint8_t)((code_devices & 0x00FF0000) >> 16);
    (*encoded_data)[7] = (uint8_t)((code_devices & 0x0000FF00) >> 8);
    (*encoded_data)[8] = (uint8_t)((code_devices & 0x000000FF));
    // Store Block Count
    (*encoded_data)[9] = (uint8_t)((block_count & 0xFF000000) >> 24);  
    (*encoded_data)[10] = (uint8_t)((block_count & 0x00FF0000) >> 16);
    (*encoded_data)[11] = (uint8_t)((block_count & 0x0000FF00) >> 8);
    (*encoded_data)[12] = (uint8_t)((block_count & 0x000000FF));
    // Store Remainder Data Devices
    (*encoded_data)[13] = (uint8_t)((remainder_data_devices & 0xFF000000) >> 24);  
    (*encoded_data)[14] = (uint8_t)((remainder_data_devices & 0x00FF0000) >> 16);
    (*encoded_data)[15] = (uint8_t)((remainder_data_devices & 0x0000FF00) >> 8);
    (*encoded_data)[16] = (uint8_t)((remainder_data_devices & 0x000000FF));
    // Store Remainder Long in Data
    (*encoded_data)[17] = (uint8_t)((remainder_long_in_data & 0xFF000000) >> 24);  
    (*encoded_data)[18] = (uint8_t)((remainder_long_in_data & 0x00FF0000) >> 16);
    (*encoded_data)[19] = (uint8_t)((remainder_long_in_data & 0x0000FF00) >> 8);
    (*encoded_data)[20] = (uint8_t)((remainder_long_in_data & 0x000000FF));
    // Store Number of Threads used
    (*encoded_data)[21] = (uint8_t)((threads & 0xFF000000) >> 24); 
    (*encoded_data)[22] = (uint8_t)((threads & 0x00FF0000) >> 16);
    (*encoded_data)[23] = (uint8_t)((threads & 0x0000FF00) >> 8);
    (*encoded_data)[24] = (uint8_t)((threads & 0x000000FF));
    // Set current index for 
    uint32_t encoded_data_index = metadata_length;

    // Setup Loop Variables
    uint32_t blocks_processed;

    // Determine number of threads per block
    int n_per_thread;
    if (block_count < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = block_count / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);
    
    // Calculate Reed-Solomon Encoding for each block
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < block_count; blocks_processed++){
        // Declare private loop variables
        int i, j;
        int remainder_block;
        int longs_processed = 0;
        int current_blocks_longs;
        uint64_t temp_l;
        uint32_t current_block_data_devices;

        // Set current block size (data devices for this block)
        if (blocks_processed == block_count-1 && remainder_data_devices != 0){
            current_block_data_devices = remainder_data_devices;
            current_blocks_longs = current_block_data_devices * 8;
            remainder_block = 1;
        } else {
            current_block_data_devices = data_devices;
            current_blocks_longs = current_block_data_devices * 8;
            remainder_block = 0;
        }

        // Calculate index offsets
        uint32_t current_block_data_index = data_devices * 8 * blocks_processed;
        uint32_t current_block_encoded_data_index = encoded_data_index + ((data_devices+code_devices)*8 + (data_devices+code_devices))*blocks_processed;
        
        // Inflate data to the correct form for Reed-Solomon Encoding
        // Create appropriate size data matrix
        char** rs_data;
        rs_data = talloc(char *, current_block_data_devices);
        for (i = 0; i < current_block_data_devices; i++){
            rs_data[i] = talloc(char, 8);

            if (longs_processed != current_blocks_longs-1 || 
               (longs_processed == current_blocks_longs-1 && blocks_processed != block_count-1) || 
               (longs_processed == current_blocks_longs-1 && blocks_processed == block_count-1 && remainder_long_in_data == 0)){
                temp_l = 0;
                for (j = 0; j < 8; j++){
                    temp_l = temp_l | ((uint64_t)data[current_block_data_index] << (64-(8*(j+1)))); 
                    current_block_data_index++;
                }
                longs_processed++;
            } else {
                temp_l = 0;
                for (j = 0; j < remainder_long_in_data; j++){
                    temp_l = temp_l | ((uint64_t)data[current_block_data_index] << (64-(8*(j+1))));
                    current_block_data_index++;
                }
                longs_processed++;
            }
            // Save inflated data
            memcpy(rs_data[i], &temp_l, 8);
        }

        // Create Coding Data Matrix
        char** rs_code;
        rs_code = talloc(char*, code_devices);
        for (i = 0; i < code_devices; i++){
            rs_code[i] = talloc(char, 8);
        }

        // Setup Vandermonde Matrix
        int* matrix = reed_sol_vandermonde_coding_matrix(current_block_data_devices, code_devices, 8);

        // Encode with classic Reed-Solomon Encoding
        jerasure_matrix_encode(current_block_data_devices, code_devices, 8, matrix, rs_data, rs_code, 8);

        // Flatten rs_code and rs_data back to a uint8_t array
        uint32_t block_stream_size = (current_block_data_devices + code_devices) * 8;
        uint8_t* block_stream = (uint8_t*)malloc(sizeof(uint8_t) * block_stream_size);

        uint32_t block_stream_index = 0;
        for (i = 0; i < code_devices; i++){
            for (j = 7; j >= 0; j--){
                block_stream[block_stream_index] = (uint8_t)rs_code[i][j]; 
                rs_code[i][j] = 0;
                block_stream_index++;
            }
        }
        for (i = 0; i < current_block_data_devices; i++){
            for (j = 7; j >= 0; j--){
                block_stream[block_stream_index] = (uint8_t)rs_data[i][j]; 
                rs_data[i][j] = 0;
                block_stream_index++;
            }
        }

        // Apply Parity protection to each device and save to protected data stream
        uint32_t bytes_processed = 0;
        uint32_t device_data_bytes = 0;
        uint8_t device_parity;
        uint64_t device_data;

        while(bytes_processed != block_stream_size){
            // Gather eight bytes of device data
            device_data = 0;
            for (j = 0; j < 8; j++){
                device_data = device_data | ((uint64_t)block_stream[device_data_bytes] << (64-(8*(j+1)))); 
                device_data_bytes++;
            }
            // Calculate parity
            device_parity = arc_calculate_parity_uint64(device_data);

            // Write parity and data out to protected data stream
            (*encoded_data)[current_block_encoded_data_index] = device_parity;
            current_block_encoded_data_index++;
            for (i = bytes_processed; i < bytes_processed + 8; i++){
                (*encoded_data)[current_block_encoded_data_index] = block_stream[i];
                current_block_encoded_data_index++;
            }
            bytes_processed = bytes_processed + 8; 
        }

        // Free allocated variables
        for (i = 0; i < current_block_data_devices; i++){
            free(rs_data[i]);
        }
        free(rs_data);
        for (i = 0; i < code_devices; i++){
            free(rs_code[i]);
        }
        free(rs_code);
        free(block_stream);
        free(matrix);
    }

    // Return resulting array
    if (PRINT){
        printf("Reed Solomon Encoding Finished\n");
    }
    return 1;
}


// arc_reed_solomon_decode:
// Decode data encoded using reed-solomon encoding
// params:
// encoded_data         -   uint8_t protected data stream
// encoded_data_size    -   size of protected data stream
// data                 -   address of pointer to decoded uint8_t data stream
// data_size            -   address of pointer to decoded data stream size
// return:
// err                  -   error code to determine if sucessful (1=sucessful,0=unsucessful)
int arc_reed_solomon_decode(uint8_t* encoded_data, uint32_t encoded_data_size, uint8_t** data, uint32_t *data_size){
    if (!INIT){
        printf("Initialization Error: Please initialize ARC before further use\n");
        return 0;
    }
    if (PRINT){
        printf("Starting Reed Solomon Decoding\n");
    }
    // Set metadata length
    uint32_t metadata_length = 25;

    // Pull metadata from encoded_data array
    uint32_t data_devices = (((uint32_t)encoded_data[1] & 0x000000FF) << 24) | (((uint32_t)encoded_data[2] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[3] & 0x000000FF) << 8) | (((uint32_t)encoded_data[4] & 0x000000FF));
    uint32_t code_devices = (((uint32_t)encoded_data[5] & 0x000000FF) << 24) | (((uint32_t)encoded_data[6] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[7] & 0x000000FF) << 8) | (((uint32_t)encoded_data[8] & 0x000000FF));
    uint32_t block_count = (((uint32_t)encoded_data[9] & 0x000000FF) << 24) | (((uint32_t)encoded_data[10] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[11] & 0x000000FF) << 8) | (((uint32_t)encoded_data[12] & 0x000000FF));
    uint32_t remainder_data_devices = (((uint32_t)encoded_data[13] & 0x000000FF) << 24) | (((uint32_t)encoded_data[14] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[15] & 0x000000FF) << 8) | (((uint32_t)encoded_data[16] & 0x000000FF));
    uint32_t remainder_long_in_data = (((uint32_t)encoded_data[17] & 0x000000FF) << 24) | (((uint32_t)encoded_data[18] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[19] & 0x000000FF) << 8) | (((uint32_t)encoded_data[20] & 0x000000FF));
    uint32_t threads = (((uint32_t)encoded_data[21] & 0x000000FF) << 24) | (((uint32_t)encoded_data[22] & 0x000000FF) << 16) | \
    (((uint32_t)encoded_data[23] & 0x000000FF) << 8) | (((uint32_t)encoded_data[24] & 0x000000FF));
    uint32_t encoded_data_index = metadata_length;

    // Check that all inputs are still valid
    // Ensure k+m <= 2^w
    if (data_devices + code_devices > ((long)1 << 8)){
        printf("INVALID PARAMETER SIZES: Please select inputs where data devices + code devices <= 2^device length. . .\n");
        return 0;
    }


    // Determine size of the original array 
    if (remainder_data_devices == 0){
        if (remainder_long_in_data == 0){
            *data_size = (data_devices * 8 * block_count);
        } else {
            *data_size = (data_devices * 8 * block_count) - (8 - remainder_long_in_data);
        }
    } else {
        if (remainder_long_in_data == 0){
            *data_size = (data_devices * 8 * (block_count-1)) + (remainder_data_devices * 8);
        } else {
            *data_size = (data_devices * 8 * (block_count-1)) + (remainder_data_devices * 8) - (8 - remainder_long_in_data);
        }
    }
    *data = (uint8_t*)malloc(sizeof(uint8_t) * *data_size);

    // Setup Loop Variables
    uint32_t blocks_processed;
    int decode_success = 1;

    // Determine number of threads per block
    int n_per_thread;
    if (block_count < threads){
        n_per_thread = 1;
    } else {
        n_per_thread = block_count / threads;
    }
    // Set number of threads
	omp_set_num_threads(threads);
    
    // Calculate Reed-Solomon Encoding for each block
    #pragma omp parallel for schedule(static, n_per_thread)
    for (blocks_processed = 0; blocks_processed < block_count; blocks_processed++){
        // Declare private loop variables
        uint32_t current_block_data_devices;
        uint32_t current_block_encoded_data_start; 
        uint32_t current_block_total_devices;
        uint32_t current_devices_processed;
        int i, j;

        // Set up erasures array
        int* erasures = talloc(int, (code_devices+1));
        int erasure_count;

        // Set number of data devices in block
        if (blocks_processed == block_count-1 && remainder_data_devices != 0){
            current_block_data_devices = remainder_data_devices;
        } else {
            current_block_data_devices = data_devices;
        }

        // Calculate index offsets
        uint32_t current_block_encoded_data_index = encoded_data_index + ((data_devices+code_devices)*8 + (data_devices+code_devices))*blocks_processed;
        uint32_t current_block_data_index = data_devices * 8 * blocks_processed;

        // Store original encoded data index for resetting later
        current_block_encoded_data_start = current_block_encoded_data_index;
        // Set total number of devices in block
        current_block_total_devices = current_block_data_devices + code_devices;

        // Check parity for each device in block
        uint8_t previous_device_parity;
        uint8_t new_device_parity;
        uint8_t xor_parity;
        uint64_t device_data;
        int parity_correct = 1;
        current_devices_processed = 0;
        erasure_count = 0;
        while(current_devices_processed != current_block_total_devices){
            // Save previous datas parity
            previous_device_parity = encoded_data[current_block_encoded_data_index];
            current_block_encoded_data_index++;
            // Gather data for current device
            device_data = 0;
            for (j = 0; j < 8; j++){
                device_data = device_data | ((uint64_t)encoded_data[current_block_encoded_data_index] << (64-(8*(j+1)))); 
                current_block_encoded_data_index++;
            }
            // Calculate parity
            new_device_parity = arc_calculate_parity_uint64(device_data);

            // Compare parities using XOR
            xor_parity = previous_device_parity ^ new_device_parity;
            // If the parities are different, set incorrect parity and count needed erasure
            if (xor_parity != 0){
                if (parity_correct == 1){
                    parity_correct = 0;
                }
                // If erasures needed pass number of coding devices, return since data cannot be fixed
                if (erasure_count < code_devices){
                    // NOTE: (The below arithmetic is due to the fact that jerasure expects data then code, 
                    // but in this code, we put code then data)

                    // If the erasure is in a data device
                    if (current_devices_processed >= code_devices){
                        erasures[erasure_count] = current_devices_processed - code_devices;
                    // If the erasure is in a code device
                    } else {
                        erasures[erasure_count] = current_devices_processed + current_block_data_devices;
                    }
                    erasure_count++;
                } else {
                    printf("DATA INTEGRITY ERROR: Too Many Erasures Needed To Fix Data, Correction Impossible...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                }
            }
            current_devices_processed++;
        }

        if (decode_success == 1){
            // If parity was correct for all items in the block, then write data out to data stream
            if (parity_correct == 1){
                // Reset data index to beginning of block and devices processed to 0
                current_block_encoded_data_index = current_block_encoded_data_start;
                current_devices_processed = 0;
                while(current_devices_processed != current_block_total_devices){
                    if (current_devices_processed >= code_devices){
                        // Skip past parity byte
                        current_block_encoded_data_index++;

                        // Write out bytes to data stream
                        for (i = 0; i < 8; i++){
                            // Make sure not to write off the end of the data array
                            if (current_block_data_index < *data_size){
                                (*data)[current_block_data_index] = encoded_data[current_block_encoded_data_index];
                                current_block_data_index++;
                                current_block_encoded_data_index++;
                            }
                        }
                    } else {
                        // Skip coding devices
                        current_block_encoded_data_index = current_block_encoded_data_index + 9;
                    }
                    current_devices_processed++;
                }
            // If parity was incorrect for any items in the block, then erase faulty devices, decode with jerasure, and write fixed data out
            } else {
                printf("Incorrect Parity Found. . .\nFixing Now\n");
                // Reset data index to beginning of block and devices processed to 0
                current_block_encoded_data_index = current_block_encoded_data_start;
                current_devices_processed = 0;
                // Set negative one to last space in erasures array
                erasures[erasure_count] = -1;

                // Create Coding Data Matrix
                char** rs_code;
                rs_code = talloc(char*, code_devices);
                for (i = 0; i < code_devices; i++){
                    rs_code[i] = talloc(char, 8);
                }

                // Inflate data and code 
                char** rs_data;
                rs_data = talloc(char *, current_block_data_devices);
                int rs_data_i = 0;
                int rs_code_i = 0;
                long temp_l = 0;
                while(current_devices_processed != current_block_total_devices){
                    // place data into rs_data
                    if (current_devices_processed >= code_devices){
                        // Initialize rs_data row 
                        rs_data[rs_data_i] = talloc(char, 8);
                        // Determine if device needs to be erased
                        int erase = 0;
                        for (i = 0; i < code_devices; i++){
                            if (erasures[i] == -1){
                                break;
                            } else if (erasures[i] == current_devices_processed - code_devices){
                                erase = 1;
                            }
                        }
                        if (erase == 1){
                            // Erase device
                            temp_l = 0;
                            memcpy(rs_data[rs_data_i], &temp_l, 8);
                            rs_data_i++;
                            // Move past current device data
                            current_block_encoded_data_index = current_block_encoded_data_index + 9;
                        } else {
                            // Skip parity byte
                            current_block_encoded_data_index++;
                            // Pull all eight bytes together 
                            temp_l = 0;
                            for (j = 0; j < 8; j++){
                                temp_l = temp_l | ((long)encoded_data[current_block_encoded_data_index] << (64-(8*(j+1)))); 
                                current_block_encoded_data_index++;
                            }
                            // Save long to rs_data
                            memcpy(rs_data[rs_data_i], &temp_l, 8);
                            rs_data_i++;
                        }
                    } else {
                        // Initialize rs_code row
                        temp_l = 0;
                        memcpy(rs_code[rs_code_i], &temp_l, 8);

                        // Determine if device needs to be erased
                        int erase = 0;
                        for (i = 0; i < code_devices; i++){
                            if (erasures[i] == -1){
                                break;
                            } else if (erasures[i] == current_devices_processed + current_block_data_devices){
                                erase = 1;
                            }
                        }
                        if (erase == 1){
                            // Erase device
                            temp_l = 0;
                            memcpy(rs_code[rs_code_i], &temp_l, 8);
                            rs_code_i++;
                            // Move past current device data
                            current_block_encoded_data_index = current_block_encoded_data_index + 9;
                        } else {
                            // Skip parity byte
                            current_block_encoded_data_index++;
                            // Pull all eight bytes together 
                            temp_l = 0;
                            for (j = 0; j < 8; j++){
                                temp_l = temp_l | ((long)encoded_data[current_block_encoded_data_index] << (64-(8*(j+1)))); 
                                current_block_encoded_data_index++;
                            }
                            // Save long to rs_code
                            memcpy(rs_code[rs_code_i], &temp_l, 8);
                            rs_code_i++;
                        }
                    }
                    // Move to next device
                    current_devices_processed++;
                }

                // Setup Vandermonde Matrix
                int* matrix = reed_sol_vandermonde_coding_matrix(current_block_data_devices, code_devices, 8);

                // Send inflated data code and erasures to jerasure decode
                int err = jerasure_matrix_decode(current_block_data_devices, code_devices, 8, matrix, 1, erasures, rs_data, rs_code, 8);

                // Check for a successful decode
                if (err == -1){
                    printf("DATA INTEGRITY ERROR: Reed Solomon Decode Has Failed...\n");
                    #pragma omp critical 
                    {
                        decode_success = 0;
                    }
                }

                if (decode_success == 1){
                    // Flatten data and write out to result data stream
                    for (i = 0; i < current_block_data_devices; i++){
                        // Write bytes to data stream
                        for (j = 7; j >= 0; j--){
                            // Make sure not to write off the end of the data array
                            if (current_block_data_index < *data_size){
                                (*data)[current_block_data_index] = (uint8_t)rs_data[i][j];
                                rs_data[i][j] = 0;
                                current_block_data_index++;
                            }
                        }
                    }
                } 

                // Free allocated variables
                for (i = 0; i < current_block_data_devices; i++){
                    free(rs_data[i]);
                }
                free(rs_data);
                for (i = 0; i < code_devices; i++){
                    free(rs_code[i]);
                }
                free(matrix);
            }
        }  

        // Free allocated variables
        free(erasures);
    }

    // Free decoded data and set data_size to 0 if decoding process failed
    if (decode_success == 0){
        free(*data);
        *data_size = 0;
        if (PRINT){
            printf("Reed Solomon Decoding Failed!\n");
        }
    } else{
        if (PRINT){
            printf("Reed Solomon Decoding Finished!\n");
        }
    }
    // Return resulting array and decode success value
    return decode_success;
}
