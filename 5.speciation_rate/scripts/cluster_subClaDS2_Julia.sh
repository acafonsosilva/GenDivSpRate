universe = vanilla
executable = /usr/local/src/julia-1.1.0/bin/julia
log = /data/biodiv/asilva/5.speciation_rate/outputs/log/log$(Cluster)_$(Process).txt
error = /data/biodiv/asilva/5.speciation_rate/outputs/log/e$(Cluster)_$(Process).txt
output = /data/biodiv/asilva/5.speciation_rate/outputs/log/o$(Cluster)_$(Process).txt
arguments = /data/biodiv/asilva/5.speciation_rate/scripts/ClaDS2_upham_4064sp_MCCposterior100.jl $(item)
notification = Never
request_cpus = 1
request_memory = 10G
nice_user = TRUE
concurrency_limits = anaTest:500
queue from seq 1 1 101 |