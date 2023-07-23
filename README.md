# lock-free-concurrent-cuckoo-filter
We propose a novel lock free concurrent cuckoo fiter, which can improve throughput on multi core platforms.
  
## Build and run
```sh
mkdir build
cd build
cmake ..
make
./test/lf_multislot_kcas_bfs_cuckoo_test
```

## Evaluation
|Algorithm| Description|
|:----:|----|
|Segment Lock Cuckoo Filter (SLCF)|Rong Gu, Simian Li, Haipeng Dai, Hancheng Wang, Yili Luo, Bin Fan, Ran Ben Basat, Ke Wang, Zhenyu Song, Shouwei Chen, Beinan Wang, Yihua Huang, Guihai Chen. 2023. Adaptive Online Cache Capacity Optimization via Lightweight Working Set Size Estimation at Scale. In Proceeding of USENIX Annual Technical Conference 2023. USENIX, 467-484.|
|Global Lock Cuckoo Filter (GLCF)|Bin Fan, David G. Andersen, Michael Kaminsky, and Michael Mitzenmacher. 2014. Cuckoo Filter: Practically Better Than Bloom. In Proceedings of ACM International Conference on Emerging Networking Experiments and Technologies. ACM, 75–88. Implementation: https://github.com/efficient/cuckoofilter|
|Vector Quotient Filter (VQF)|P. Pandey, A. Conway, J. Durie, M. A. Bender, M. Farach-Colton, and R. Johnson, “Vector quotient filters: Overcoming the time/space trade-off in filter design,” in Proceedings of International Conference on Management of Data. ACM, 2021, pp. 1386–1399. Implementation: https://github.com/splatlab/vqf|
|One Hash Blocked Bloom Filters (OHBBF)|Elakkiya Prakasam and Arun Manoharan. 2022. A Cache Efficient One Hashing Blocked Bloom Filter (OHBB) for Random Strings and the K-mer Strings in DNA Sequence. Symmetry 14, 9 (2022), 1–24.|


## To generate YCSB workloads
```sh
cd YCSB
wget https://github.com/brianfrankcooper/YCSB/releases/download/0.17.0/ycsb-0.17.0.tar.gz
tar -xvf ycsb-0.17.0.tar.gz
mv ycsb-0.17.0 YCSB
#Then run workload generator
mkdir workloads
./generate_all_workloads.sh
```