/* file: chaosMap.nl
 * purpose: namelist definition for chaos map analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist chaos_map static
    STRING output = NULL;
    double xmin =  -0.1;
    double xmax =  0.1;
    double ymin =  1e-6;
    double ymax =  0.1;
    double delta_min = 0;
    double delta_max = 0;
    long nx = 20;
    long ny = 21;
    long ndelta = 1;
    long forward_backward = 0;
    double change_x = 1e-6;
    double change_y = 1e-6;
    long verbosity = 1;
#end

