#!julia

#
#  This is the main file for data assimilation runs.
#
#
# 1. loads up a configuration file,
# 2. obtains data from a WRF model,
# 3. reads in observations and metadata for a list of stations,
# 4. runs the moisture model and the assimilation mechanism.
#
#

using Storage
using Kriging
using WRF
using FM




function main(args)

    # the arguments passed to the julia program do not include program name
    if length(args) != 1
        println("Usage: julia run_data_assimilation.jl cfg_file")
        exit(1)
    end
    
    cfg = evalfile(args[1])

    # configure Storage mechanism
    Storage.sopen(cfg["output_dir"], "moisture_model_v2_diagnostics.txt", "frame")

    # Storage tags
    Storage.setup_tag("assim_K0", false, true, true)
    Storage.setup_tag("assim_K1", true, true, true)
    Storage.setup_tag("assim_data", false, false, true)

    Storage.setup_tag("obs_residual_var", true, true, true)

    Storage.setup_tag("fm10_model_residual_var", true, true, true)
    Storage.setup_tag("fm10_model_var", false, true, true)
    Storage.setup_tag("fm10_kriging_var", false, true, true)

    #read in data from the WRF output file pointed to by cfg
    

end


main(ARGS)
