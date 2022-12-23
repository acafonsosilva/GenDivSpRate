universe = vanilla
arguments = --vanilla --args $(item)
log = /data/biodiv/asilva/SuperCrunchClean/bmlm/global/log/$(Cluster)_$(Process)log.txt
error = /data/biodiv/asilva/SuperCrunchClean/bmlm/global/log/$(Cluster)_$(Process)e.txt
output = /data/biodiv/asilva/SuperCrunchClean/bmlm/global/log/$(Cluster)_$(Process)o.txt
executable = /users/biodiv/asilva/R-3.5.0/bin/R
input = /data/biodiv/asilva/SuperCrunchClean/bmlm/scripts/runBMLM_globalAll.R
notification = Never
request_cpus = 1
request_memory = 3G
nice_user = TRUE
getenv = True
Queue from seq 1 1 101 |