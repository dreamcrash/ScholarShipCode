* The C version uses usigned char to replace the byte type used in the JGF JAVA original implementation;
* Because of the above point, operations like >>> (of the original JAVA version)  where replace to >> (in the C version);
* The SM version with the first touch memory allocation policy provided gains ranging from 1.17 to 1.02x over the
the version SM without that policy:
	* tested in machine with two E5-2650 v2 processors5 (NUMA) fixed frequency of 2.6 GHz;
	* machine from the SeARCH cluster; run with CentOS 6.3 and GCC 4.9.3;
	* More informations:
		* https://ark.intel.com/products/75269/Intel-Xeon-Processor-E5-2650-v2-20M-Cache-2_60-GHz.
		* search.di.uminho.pt.
