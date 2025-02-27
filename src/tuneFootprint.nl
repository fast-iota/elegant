/* file: tuneFootprint.nl
 * purpose: namelist definition for tune footprint analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist tune_footprint static
    STRING delta_output = NULL;
    STRING xy_output = NULL;
    double xmin =  -0.02;
    double xmax =  0.02;
    double ymin =  1e-6;
    double ymax =  0.02;
    double x_for_delta = 1e-6;
    double y_for_delta = 1e-6;
    double separate_xy_for_delta = 0;
    double delta_min = 0;
    double delta_max = 0;
    long nx = 20;
    long ny = 21;
    long ndelta = 21;
    long verbosity = 1;
    long quadratic_spacing = 1;
    long compute_diffusion = 1;
    long diffusion_rate_limit = -5;
    long immediate = 0;
    long filtered_output = 1;
    long ignore_half_integer = 0;
    long chromaticity_fit_order = 5;
#end

