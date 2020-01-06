# Data Redundacy approach (the fastest of the SM versions):

* This approach significantly increases the performance of the SM parallelization by reducing the synchronization overhead and improving load balacing;
* This approach instead of using synchronization to deal with the variables susceptible to race conditions (i.e., fx, fy, fz, epot, vir and interactions) replicates those variables per thread. At the end of the parallel region the values of those variables are reduced across all threads.
	- In the 'ManualReduction' version the reduction process is done manually, whereas in the 'OpenMPReductions' version is done using OpenMP 4.0 and 4.5 reduction capabilities.
	- The 'OpenMPReductions' versions: 
		- needs to be compiled with GCC 6.1+;
		- for the largest inputs it might be necessary to tune the threads stack size:
			-> For instance with "ulimit -s <stacksize>" and OMP_STACKSIZE (e.g., 10m);
	- The 'ManualReduction' version uses less memory than the 'OpenMPReductions' since the master thread uses the original arrays of forces;
	- The 'ManualReduction' version is more verbose and less readable than the 'OpenMPReductions' version;

* The load balancing was improved by dividing the force calculation between particles among threads in a dynamic fashion.
