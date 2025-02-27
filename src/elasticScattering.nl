/* file: elasticScattering.nl
 * contents: namelist for elastic scattering simulation
 * 
 * Michael Borland, 2017
 */
#include "namelist.h"

#namelist elastic_scattering static
    STRING losses = NULL;
    STRING output = NULL;
    STRING log_file = NULL;
    double theta_min = 0.001;
    double theta_max = 0.010;
    long n_theta = 11;
    long quadratic_theta_spacing = 0;
    long n_phi = 37;
    long twiss_scaling = 0;
    double s_start = 0;
    double s_end = DBL_MAX;
    STRING include_name_pattern = NULL;
    STRING include_type_pattern = NULL;
    long verbosity = 1;        
    long soft_failure = 0;
    long allow_watch_file_output = 0;
#end


