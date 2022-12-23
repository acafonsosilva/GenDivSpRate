universe = vanilla
executable = /users/biodiv/asilva/R-3.6.3/bin/R
input = /data/biodiv/asilva/SuperCrunchClean/pgls/scripts/runPGLS.R
arguments = --vanilla --args $(item)
log = /data/biodiv/asilva/SuperCrunchClean/pgls/log/$(Cluster)_$(Process)log.txt
error = /data/biodiv/asilva/SuperCrunchClean/pgls/log/$(Cluster)_$(Process)e.txt
output = /data/biodiv/asilva/SuperCrunchClean/pgls/log/$(Cluster)_$(Process)o.txt
notification = Never
request_cpus = 1
request_memory = 2G
nice_user = TRUE
queue from seq 1 1 101 |
