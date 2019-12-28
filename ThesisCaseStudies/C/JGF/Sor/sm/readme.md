# First Touch Improvement:
- Applying the first touch into the C Shared Memory SOR version, allowed C code to surpass the equivalente version in JAVA (using XX:+UseNUMA flag);
- In our tests using the following machines:
	- Intel Nehalem CPU Hex-Core, 2.67 GHz, 12 MB Cache L3, corresponding to 12 cores and 24 threads.
- For the three largest inputs (i.e., 2k, 10k and 15k) the first touch approach reduced the execution time by approximately 2x;
- Besides the first touch we also improved the C version by allocating the SOR matrix continuously in memory;

- The validation method does not work correctly for multi-threads with the first touch technique since the input matrix
will be generated in a non-deterministic way.

