
#define LOCAL_SIZE 65536
#define N_STATES 5
#define N_RULES (N_STATES * (N_STATES+1) / 2 * (N_STATES-1))
#define MAX_SIZE 12
#define MAX_ITERS 50

static inline char apply_automaton(char *automaton, char left, char middle, char right) {
    int l = min(left, right);
    int r = max(left, right);
    int index = (middle - (middle == N_STATES - 1)) * (N_RULES / (N_STATES - 1));
    index += ((N_STATES * (N_STATES+1)) - ((N_STATES-l) * (N_STATES-l+1))) / 2;
    index += r - l;
    return automaton[index];
}

static inline int simulate(int size, char *automaton) {
    char a[MAX_SIZE];
    char b[MAX_SIZE];
    char *system = a;
    char *last_system = b;
    char *temp;

    last_system[0] = 1;
    for (int i = 1; i < size; i++)
        last_system[i] = 0;

    int first_end_state = 0;
    int first_all_end_state = 0;

    for (int y = 0; y < MAX_ITERS; y++) {
        int end_state_count = 0;
        system[0] = apply_automaton(automaton, N_STATES-1, last_system[0], last_system[1]);
        end_state_count += system[0] == N_STATES-1;
        for (int x = 1, e = size-1; x < e; x++) {
            system[x] = apply_automaton(automaton, last_system[x-1], last_system[x], last_system[x+1]);
            end_state_count += system[x] == N_STATES-1;
        }
        system[size-1] = apply_automaton(automaton, last_system[size-2], last_system[size-1], N_STATES-1);
        end_state_count += system[size-1] == N_STATES-1;
        first_end_state = (first_end_state != 0) * first_end_state + (first_end_state == 0) * y * (end_state_count != 0);
        first_all_end_state = (first_all_end_state != 0) * first_all_end_state + (first_all_end_state == 0) * y * (end_state_count == size);
        temp = system;
        system = last_system;
        last_system = temp;
    }

    return (first_all_end_state != 0) * (first_all_end_state == first_end_state);
}

// combination_count = (N_STATES - 1) ^ rule_changes * nCr(N_RULES - 1, rule_changes)
__kernel void start(ulong offset, int rule_changes, ulong combination_count, __constant char *base_automaton,
       __global char *global_results, __global size_t *global_successes) {

    size_t global_id = get_global_id(0);
    size_t local_id = get_local_id(0);
    size_t local_size = get_local_size(0);
    size_t global_index = global_id / local_size;

    // Combination to rule changes
    size_t combination_id = offset + global_id;
    if (combination_id > combination_count)
        return;

    char automaton[N_RULES];
    for (int i = 0; i < N_RULES; i++)
        automaton[i] = base_automaton[i];

    int min_index = 0;
    for (int i = 0; i < rule_changes; i++) {
        int possible_indices = (N_RULES - 1) - min_index - (rule_changes - 1 - i);
        int index = 1 + min_index + combination_id % possible_indices;
        combination_id /= possible_indices;
        int dest_state = combination_id % (N_STATES - 1);
        dest_state += dest_state >= automaton[index];
        combination_id /= (N_STATES - 1);
        min_index = index + 1;
        automaton[index] = dest_state;
    }

    // Simulate
    int successes = simulate(7, automaton);
    successes += simulate(8, automaton);
    successes += simulate(10, automaton);
    successes = (successes + 1) * simulate(9, automaton);

    __local char local_results[LOCAL_SIZE];
    __local size_t local_successes[LOCAL_SIZE];

    // Search for the result with the most successes
    local_results[local_id] = successes;
    local_successes[local_id] = offset + global_id;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i = local_size / 2; i > 0; i >>= 1) {
        if (local_id < i && local_results[local_id + i] > local_results[local_id]) {
            local_results[local_id] = local_results[local_id + i];
            local_successes[local_id] = local_successes[local_id + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (local_id == 0) {
        global_results[global_index] = local_results[0];
        global_successes[global_index] = local_successes[0];
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
    for (int i = get_global_size(0) / local_size / 2; i > 0; i >>= 2) {
        if (local_id == 0 && global_index < i && global_results[global_index + i] > global_results[global_index]) {
            global_results[global_index] = global_results[global_index + i];
            global_successes[global_index] = global_successes[global_index + i];
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}
