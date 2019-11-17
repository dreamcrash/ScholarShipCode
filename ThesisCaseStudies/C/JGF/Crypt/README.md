* The C version uses usigned char to replace the byte type used in the JGF JAVA original implementation;
* Because of the above point, operations like >>> (of the original JAVA version)  where replace to >> (in the C version);
* The SM version with the first touch memory allocation policy provided gains ranging from 1.17 to 1.02x over the
the version SM without that policy. 
