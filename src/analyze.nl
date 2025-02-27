/* file: analyze.nl
 * purpose: namelist definition for transport analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist analyze_map static
    STRING output = NULL;
    long output_order = 1;
    STRING printout = NULL;
    long printout_order = 2;
    STRING printout_format = "%22.15e ";
    double delta_x = 5e-5;
    double delta_xp = 5e-5;
    double delta_y = 5e-5;
    double delta_yp = 5e-5;
    double delta_s  = 5e-5;
    double delta_dp = 5e-5;
    double accuracy_factor = 1e-12;
    long center_on_orbit = 0;
    long verbosity = 0;
    long n_points = 9;
    long max_fit_order = 8;
    long canonical_variables = 0;
    long periodic = 1;
    double beta_x = 1;
    double alpha_x = 0;
    double eta_x = 0;
    double etap_x = 0;
    double beta_y = 1;
    double alpha_y = 0;
    double eta_y = 0;
    double etap_y = 0;
#end

