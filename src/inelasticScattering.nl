/* file: inelasticScattering.nl
 * contents: namelist for inelastic scattering simulation
 * 
 * Michael Borland, 2017
 */
#include "namelist.h"

#namelist inelastic_scattering static
    STRING losses = NULL;
    STRING output = NULL;
    STRING log_file = NULL;
    double k_min = 0.001;
    STRING momentum_aperture = NULL;
    double momentum_aperture_scale = 0.90;
    double momentum_aperture_periodicity = 0;
    long n_k = 101;
    double s_start = 0;
    double s_end = DBL_MAX;
    STRING include_name_pattern = NULL;
    STRING include_type_pattern = NULL;
    long verbosity = 1;        
    long soft_failure = 0;
    long allow_watch_file_output = 0;
#end


