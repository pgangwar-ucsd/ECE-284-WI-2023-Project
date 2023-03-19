## Note: Might need a supercomputer cluster to run this code as it requires more than 200 GB memory. This is the reason placement of reads could not complete on AWS (tried on t3.micro, t2.small).

Steps to run:
1. git clone https://github.com/pgangwar-ucsd/ECE-284-WI-2023-Project.git
2. cd ECE-284-WI-2023-Project
3. rm -rf build

Option-1: Follow the build steps mentioned in https://usher-wiki.readthedocs.io/en/latest/Installation.html

Option 2:
4. ./install/installUbuntu.sh
5. rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
6. chmod +x faToVcf
7. mv faToVcf scripts/

8. source run.sh

NOTE: run.sh can be edited to generate different results by the following switches
-r 150 : read_len of 150
-w 20 : using 20 samples from a clade
-s 2 : sequencing depth of 2
-T * : currenty not mentioned in run.sh. Including this value allows us to run the code on specified number of cores 



Key Points: There is a file named "experiment.cpp" in "src/matUtils/" that contains the main code. I have implemented two function for the tree search "place_reads_nodes_sequential" and "place_reads_nodes_parallel". Unfortunately both cannot be run at same time. Doing so consumes a lot of memory. Currently "place_reads_nodes_parallel" is commented on line number 950 in "experiment.cpp". That means the current code gives timing of "place_reads_nodes_sequential" that is the main proposed algorithm of this project. If you need to get placement timing of "place_reads_nodes_parallel", you need to comment line 949 where it calls "place_reads_nodes_sequential" and uncomment line 950 where it calls "place_reads_nodes_parallel".
