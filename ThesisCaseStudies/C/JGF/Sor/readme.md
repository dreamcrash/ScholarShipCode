# Some remarks:
- The readbility of the SOR adapted versions (e.g., original and improved versions) can be signficantly improved;
- The SOR MPI distribution can be improved: - in the way the iterations are assigned to the
processes; and regarding robutness:
	- it runs into problems when executed with more than 32 processes;
	- it was fully tested with 1, 2, 4, 6, 8, 12, 16, 24 and 32 processes. It might not 
validate correctly with other numbers of processes.
