universe = vanilla
executable = /bin/bash
arguments = /data/biodiv/asilva/mutationRate/baseml.sh $(item)
log = /data/biodiv/asilva/mutationRate/log/log$(Cluster)_$(Process).txt
error = /data/biodiv/asilva/mutationRate/log/e$(Cluster)_$(Process).txt
output = /data/biodiv/asilva/mutationRate/log/o$(Cluster)_$(Process).txt
notification = Never
request_cpus = 1
request_memory = 1G
nice_user = TRUE
queue from seq 1 1 100 |
#queue from seq 1 1 25 |  ##To parallelise MCC tree for each of the 25 clades