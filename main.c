#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#ifdef __linux
#include <unistd.h>
#endif
#include <errno.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#include <CL/cl.h>

#define CLEAR_LINE "\x1b[K"

#include "clutil.h"

#define TYPE CL_DEVICE_TYPE_ALL

size_t flp2 (size_t x) {
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    return x - (x >> 1);
}

static inline uint64_t get_timer() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((uint64_t) ts.tv_sec) * 1000000000ULL + ts.tv_nsec;
}

static inline uint64_t start_timer() {
    return get_timer();
}

static inline uint64_t stop_timer(uint64_t start) {
    return get_timer() - start;
}

#define LOCAL_SIZE 65536
#define N_STATES 5
#define N_RULES (N_STATES * (N_STATES+1) / 2 * (N_STATES-1))
#define max(a, b) ((a) > (b) ? (a) : (b))

const char DAN_AUTOMATON[N_STATES-1][N_STATES][N_STATES] = {{{0,1,0,1,0},{1,3,2,0,3},{0,2,0,3,0},{1,0,3,0,0},{0,3,0,0,0}},{{1,1,0,3,2},{1,1,0,3,0},{0,0,3,0,0},{3,3,0,0,0},{2,0,0,0,0}},{{2,3,2,0,0},{3,2,3,3,2},{2,3,1,0,0},{0,3,0,3,3},{0,2,0,3,0}},{{2,2,2,2,2},{2,0,2,2,0},{2,2,3,0,0},{2,2,0,4,4},{2,0,0,4,0}}};

void encode_automaton(const char *multi_dimensional_array, char *dest) {
    for (int middle = 0; middle < N_STATES - 1; middle++) {
        for (int left = 0; left < N_STATES; left++) {
            for (int right = 0; right < N_STATES; right++) {
                int l = min(left, right);
                int r = max(left, right);
                int index = middle * (N_RULES / (N_STATES - 1));
                index += ((N_STATES * (N_STATES+1)) - ((N_STATES-l) * (N_STATES-l+1))) / 2;
                index += r - l;
                dest[index] = *(multi_dimensional_array + right + N_STATES * (left + N_STATES * middle));
            }
        }
    }
}

// Returns factorial of n
size_t factorial(size_t n) {
    size_t res = 1;
    for (size_t i = 2; i <= n; i++)
        res = res * i;
    return res;
}

size_t nCr(size_t n, size_t r) {
    return factorial(n) / (factorial(r) * factorial(n - r));
}

size_t pow(size_t x, size_t n) {
    size_t result = 1;
    for (size_t i = 0; i < n; i++)
        result *= x;
    return result;
}


int main(int argc, char **argv) {

    char automaton[N_RULES];
    encode_automaton(&DAN_AUTOMATON[0][0][0], automaton);

    // Settings
    // Whenever you add a CL file, remember to also edit the bottom of CMakeLists.txt
    const char *kernel_file = "kernel.cl";
    const char *kernel_name = "start";
    const char *header_names[] = {};


#ifdef __linux
    printf("PID: %u\n", getpid());
#endif

    // Find devices
    cl_uint num_platforms;
    check(clGetPlatformIDs(0, NULL, &num_platforms), "getPlatformIDs (num)");
    printf("%d platforms:\n", num_platforms);

    cl_platform_id platforms[num_platforms];
    check(clGetPlatformIDs(num_platforms, platforms, NULL), "getPlatformIDs");

    cl_device_id device = NULL;
    cl_uint cus = 0;

    printf("Available platforms:\n");
    for (int i = 0; i < num_platforms; i++) {
        char *info = getPlatformInfo(platforms[i]);
        printf("%d: %s\n", i, info);
        free(info);
        cl_uint num_devices;
        check(clGetDeviceIDs(platforms[i], TYPE, 0, NULL, &num_devices), "getDeviceIDs (num)");

        cl_device_id devices[num_devices];
        check(clGetDeviceIDs(platforms[i], TYPE, num_devices, devices, NULL), "getDeviceIDs");

        printf("  %d available devices:\n", num_devices);
        for (int j = 0; j < num_devices; j++) {
            device_info *infos = getDeviceInfo(devices[j]);
            printf("    %s\n", infos->info_str);
            if (infos->compute_units >= cus) {
                cus = infos->compute_units;
                device = devices[j];
            }

        }
        putchar('\n');
    }

    if (!device) {
        fprintf(stderr, "No devices found.\n");
        return -1;
    }

    // Create CL context
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, NULL);
    device_info *info = getDeviceInfo(device);
    printf("Using %s\n", info->info_str);
    free(info);
    cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, 0, NULL);

    // Compile headers
    const unsigned header_count = sizeof(header_names) / sizeof(*header_names);
    cl_int error;
    cl_program headers[header_count];
    for (int i = 0; i < header_count; i++) {
        const char *header_src = readFile(header_names[i]);
        headers[i] = clCreateProgramWithSource(context, header_count, &header_src, NULL, &error);
        check(error, header_names[i]);
        error = clCompileProgram(headers[i], 0, NULL, NULL, 0, NULL, NULL, NULL, NULL);
        if (error == CL_COMPILE_PROGRAM_FAILURE || error == CL_BUILD_PROGRAM_FAILURE) {
            size_t log_size;
            clGetProgramBuildInfo(headers[i], device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            char *log = (char *) malloc(log_size);
            clGetProgramBuildInfo(headers[i], device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            printf("%s\n", log);
            exit(error);
        } else if (error) {
            printf("Error compiling: %s\n", getErrorString(error));
            exit(error);
        }
    }

    // Compile main source
    const char *main_src = readFile(kernel_file);
    cl_program main_cl = clCreateProgramWithSource(context, 1, &main_src, NULL, &error);
    check(error, "Creating program");
    printf("Compiling...\n");
    error = clCompileProgram(main_cl, 0, NULL, NULL, header_count, header_count ? headers : NULL, header_count ? header_names : NULL, NULL, NULL);
    if (error == CL_COMPILE_PROGRAM_FAILURE || error == CL_BUILD_PROGRAM_FAILURE) {
        size_t log_size;
        clGetProgramBuildInfo(main_cl, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char *log = (char *) malloc(log_size);
        clGetProgramBuildInfo(main_cl, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        printf("%s j\n", log);
        exit(error);
    } else if (error) {
        fprintf(stderr, "Error compiling: %s\n", getErrorString(error));
        exit(error);
    }

    // Link
    const cl_program programs[] = {main_cl};
    printf("Linking...\n");
    cl_program program = clLinkProgram(context, 0, NULL, NULL, 1, programs, NULL, NULL, &error);
    check(error, "Linking");
    cl_kernel kernel = clCreateKernel(program, kernel_name, &error);
    check(error, "Create kernel");



    // Main program, host side
    cl_uint address_bits;
    check(clGetDeviceInfo(device, CL_DEVICE_ADDRESS_BITS, sizeof(address_bits), &address_bits, NULL), "Get device address bits");
    size_t max_global_size = 1LLU << min(address_bits, 32);

    printf("Local size: %ld; Global size: %ld\n", LOCAL_SIZE, max_global_size);

    cl_mem mem_base_automaton = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(automaton), NULL, NULL);

    cl_mem mem_global_results = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char) * (max_global_size / LOCAL_SIZE), NULL, NULL);
    cl_mem mem_global_successes = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(size_t) * (max_global_size / LOCAL_SIZE), NULL, NULL);

    check(clEnqueueWriteBuffer(queue, mem_base_automaton, CL_FALSE, 0, sizeof(automaton), automaton, 0, NULL, NULL), "Write Base Automaton");

    check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &mem_base_automaton), "Base Automaton");
    check(clSetKernelArg(kernel, 4, sizeof(cl_mem), &mem_global_results), "Global Results");
    check(clSetKernelArg(kernel, 5, sizeof(cl_mem), &mem_global_successes), "Global Successes");

    char *kf = malloc(strlen(kernel_file));
    strcpy(kf, kernel_file);
    kf[strlen(kf) - 3] = 0;
    printf("Executing %s.%s\n", kf, kernel_name);

    for (int rule_changes = 1; rule_changes < 7; rule_changes++) {

        size_t combination_count = pow(N_STATES - 1, (size_t)rule_changes) * nCr(N_RULES - 1, (size_t)rule_changes);
        check(clSetKernelArg(kernel, 1, sizeof(int), &rule_changes), "Rule changes");
        check(clSetKernelArg(kernel, 2, sizeof(size_t), &combination_count), "Combination count");
        size_t global_size = min(max_global_size, (combination_count & (combination_count-1)) == 0 ? combination_count : flp2(combination_count * 2));

        //printf("Work Load possible: %d %d %d \n",CL_DEVICE_MAX_WORK_GROUP_SIZE,CL_DEVICE_MAX_WORK_ITEM_SIZES,CL_DEVICE_LOCAL_MEM_SIZE);

        char best_result = 0;
        size_t best_combination_id = 0;

        uint64_t t = start_timer();
        for (size_t offset = 0; offset < combination_count; offset += global_size) {
            float perc = offset * 100.f / combination_count;
            uint64_t t2 = start_timer();
            check(clSetKernelArg(kernel, 0, sizeof(offset), &offset), "Argument offset");
            printf("\rx  %3.3f%%", perc);
            fflush(stdout);
            size_t local_size = LOCAL_SIZE;
            check(clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL), "\nExecute");
            check(clFinish(queue), "\nFinish execute");
            printf("\r<- %3.3f%%", perc);
            fflush(stdout);
            char result;
            size_t combination_id;
            check(clEnqueueReadBuffer(queue, mem_global_results, CL_TRUE, 0, sizeof(result), &result, 0, NULL, NULL), "\nRead");
            check(clEnqueueReadBuffer(queue, mem_global_successes, CL_TRUE, 0, sizeof(combination_id), &combination_id, 0, NULL, NULL),
                  "\nRead");
            uint64_t d2 = stop_timer(t2);
            printf("\rw  %3.3f%% %.4f" CLEAR_LINE, perc, d2 / 1000000.f);
            fflush(stdout);

            if (result > best_result) {
                best_result = result;
                best_combination_id = combination_id;
            }

            uint64_t d = get_timer() - t;
            double per_item = (double) d / (offset + global_size);
            double eta = ((double) d / ((offset + global_size) * 1000000000.)) * (combination_count - offset - global_size);
            uint64_t d2ms = d2 / 1000000;
            printf("  %fs / %llu items = %lfns/item, %llums/batch, ETA: %lfs (%dh%dm%ds)", d / 1000000000.f,
                   offset + global_size,
                   per_item, d2ms, eta,
                   (int) (eta / 3600), ((int) (eta / 60)) % 60, ((int) eta) % 60);
        }

        printf("Best combination for %d rule changes:\n", rule_changes);
        printf("Combination ID: %ld\n", best_combination_id);
        size_t combination_id = best_combination_id;
        int min_index = 0;
        for (int i = 0; i < rule_changes; i++) {
            int possible_indices = (N_RULES - 1) - min_index - (rule_changes - 1 - i);
            int index = 1 + min_index + combination_id % possible_indices;
            combination_id /= possible_indices;
            int dest_state = combination_id % (N_STATES - 1);
            dest_state += dest_state >= automaton[index];
            combination_id /= (N_STATES - 1);
            min_index = index + 1;

            int working_index = index;
            int middle = working_index / (N_RULES / (N_STATES - 1));
            working_index %= N_RULES / (N_STATES - 1);
            int left, right;
            for (int row_size = 5; row_size > 0; row_size--) {
                if (row_size > working_index) {
                    left = 5 - row_size;
                    right = left + working_index;
                }
                working_index -= row_size;
            }
            printf("(%d, %d, %d) (index %d) -> %d\n", left, middle, right, index, dest_state);
        }

        puts("\rDone.                                                       ");
        uint64_t d = stop_timer(t);
        printf("%fs / %llu items = %lfns/item\n", d / 1000000000.f, combination_count, (double) d / combination_count);
    }

    return 0;
}
