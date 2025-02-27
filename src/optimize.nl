/* file: optimize.nl
 * purpose: namelists for optimizing beamline
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

static char *optimize_mode[N_OPTIM_MODES] = {
    "minimize", "maximize"
    } ;

static char *optimize_method[N_OPTIM_METHODS] = {
    "simplex", "grid", "sample", "powell", "randomsample", "randomwalk", "genetic", "hybridsimplex", "swarm", "1dscans", "rcds",
    } ;

      
#if USE_MPI
/* This is needed by the parallel version. There will be a warning that it is unused when compiling the serial version */
static char *crossover_type[N_CROSSOVER_TYPES] = {
    "onepoint", "twopoint", "uniform"
    } ;
#endif
	
#namelist optimization_term,struct
    STRING term = NULL;
    double weight = 1.0;
    STRING field_string = NULL;
    long field_initial_value = 0;
    long field_final_value = 0;
    long field_interval = 1;
    STRING input_file = NULL;
    STRING input_column = NULL;
    long verbose = 0;
#end

#namelist optimization_setup static
    STRING equation = NULL;
    STRING mode = "minimize";
    STRING method = "simplex";
    STRING statistic = "sum-absolute-value";
    double tolerance = -0.01;
    double hybrid_simplex_tolerance = -0.01;
    double hybrid_simplex_tolerance_count = 2;
    long hybrid_simplex_comparison_interval = 0;
    double target = -DBL_MAX;
    long center_on_orbit = 0;
    long center_momentum_also = 1;
    long soft_failure = 1;
    long n_passes = 2;
    long n_evaluations = 500;
    long n_restarts = 0;
    double restart_reset_threshold = 0;
    double restart_worst_term_factor = 1;
    long restart_worst_terms = 1;
    long matrix_order = 1;
    STRING log_file = NULL;
    STRING term_log_file = NULL;
    long verbose = 1;
    long output_sparsing_factor = 1;
    STRING crossover = "twopoint";
    long balance_terms = 0;
    double simplex_divisor = 3;
    double simplex_pass_range_factor = 1;
    long include_simplex_1d_scans = 1;
    long start_from_simplex_vertex1 = 0;
    long restart_random_numbers = 0;
    double random_factor = 1.0;		
    double rcds_step_factor = 0.1;
    long n_iterations = 10000;
    long max_no_change = 1000;
    long population_size = 100;
    long print_all_individuals = 0;	
    STRING population_log = NULL;	
    STRING interrupt_file = "%s.interrupt";
    double interrupt_file_check_interval = 10;
    long inhibit_tracking = 0;
    STRING simplex_log = NULL;
    long simplex_log_interval = 1;
#end

#namelist optimization_variable static
    STRING name = NULL;
    STRING item = NULL;
    double lower_limit = 0;
    double upper_limit = 0;
    double step_size = 1;
    double fractional_step_size = -1;
    long disable = 0;
    long force_inside = 0;
    long differential_limits = 0;
    long no_element = 0;
    double initial_value = 0;
#end

#namelist optimization_constraint static
    STRING quantity = NULL;
    double lower = 0;
    double upper = 0;
#end

#namelist optimize static
     long summarize_setup = 0;
#end

#namelist set_reference_particle_output,struct
          STRING match_to = NULL;
          double weight[6] = {1, 1, 1, 1, 0, 1};
          STRING comparison_mode = NULL;
#end
