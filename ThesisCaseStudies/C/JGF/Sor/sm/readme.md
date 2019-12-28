# First Touch Memory Allocation Policy:
- Applying the first touch approach into the C Shared Memory SOR version allowed the C code to surpass the equivalente version in JAVA (using XX:+UseNUMA flag);
- In our bechmarks we used the following machine:
	- Intel Nehalem CPU Hex-Core, 2.67 GHz, 12 MB Cache L3, corresponding to 12 cores and 24 threads.
- For the three largest inputs (i.e., 2k, 10k and 15k) the first touch approach reduced the execution time by approximately 2x;
- Besides the first touch we also improved the C version by allocating the SOR matrix continuously in memory;

# Important Note:
- The validation method does not work correctly for multi-threads with the first touch technique since the input matrix
will be generated in a non-deterministic way.
	- Workaround : When running in validation mode one should initialize the matrix G sequentially (i.e., using #pragma omp master instead of #pragma omp for nowait).

