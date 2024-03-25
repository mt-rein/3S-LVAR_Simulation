# This script defines auxiliary functions to be used in the simulation
# specifically, these functions catch errors and warnings
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))
run_NFS <- quietly(safely(sem))
run_SAM <- quietly(safely(sam))
run_SEM <- quietly(safely(sem))