# Some remarks:
- The SOR original MPI version can be signficantly improved in terms of readbility;
- The SOR MPI distribution can be optimized in the way the iterations are assigned to the
processes. Moreover, the distribution can be improved in terms of robutness:
	- it runs into problems when executed with more than 32 processes;
	- it was fully tested with 1, 2, 4, 6, 8, 12, 16, 24 and 32 processes. It might not 
validate correctly with other numbers of processes.
