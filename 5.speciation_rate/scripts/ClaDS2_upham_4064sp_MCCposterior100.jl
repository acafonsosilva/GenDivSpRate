##/usr/local/src/julia-1.1.0/bin/julia
using ArgParse
using RCall

function parse_commandline(args)
    s = ArgParseSettings()

    @add_arg_table s begin
        "arg1"
            help = "for julia script"
            required = true
    end
    pa = parse_args(args, s)
    h = pa["arg1"]
    return(h)
end

a = parse_commandline(ARGS)
@rput a
println(a)

folder_path = /data/biodiv/asilva/5.speciation_rate/
@rput folder_path

include(joinpath(folder_path,"scripts/ClaDS_Julia_scripts/load_ClaDS2_functions.jl"))
include(joinpath(folder_path,"scripts/ClaDS_Julia_scripts/mcmc_ClaDS2.jl"))
include(joinpath(folder_path,"scripts/ClaDS_Julia_scripts/load_ETR2.jl"))
include(joinpath(folder_path,"scripts/ClaDS_Julia_scripts/LTT.jl"))


reval("""library(RPANDA)
    load(paste0(folder_path,"outputs/upham_4064sp_FR_MCCposterior100.rdata"))
    tree_full=TreeSet[[as.numeric(a)]]
""")

@rget tree_full

tree = ape2Tree(tree_full);

reval("""
        
        f_tips = c()
        for (sp in tree_full[[4]]){
            if ( !sp %in% species_SF[,2]){
                f_tips = c(f_tips,1)
            }else{
                f_tips = c(f_tips,species_SF[which(species_SF[,2]==sp),3])
            }
        }
    """)
@rget f_tips
fs = sample_fractions_allDeep(tree, f_tips)[2:end]# to have 1. as a sampling fraction in the backbone tree.

pathSet = joinpath(folder_path,"outputs/upham_tree$(a).Rdata");

sampler = run_ClaDS2(tree, 400, f=fs, 
    plot_chain = false, save_as_R=true, Rfile=pathSet,
    plot_tree=0, goal_gelman = 1.05, end_it = 50_000)

