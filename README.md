# Investigating priority scheduling algorithms
Project by Ido Kessler, Tal Leibovitch and Gilad Fudim.

Priority queues are data structures which store keys in an ordered fashion to allow efficient access to the minimal (maximal) key. They are a data structure for maintaining a set of keys (potentially stored along with some data) such that the minimum (or maximum) key can be efficiently accessed and removed.

In this project we implemented the CA-PQ as priority queue in Galois.
Benchmarked and investigated why the k-LSM present such unfavorite results in real graphs.
We suggested few ideas on how to improve the efficiency of the queue, and implemented the idea we thought will have the most impactful results, swapping the inner shared data-structure which was being used with a MultiQueue with “peek” operation we implemented.
Benchmarked all of the above against known queues on real graphs.
Compared the queues against each other and summarized our conclusions for our project and additonal future work.