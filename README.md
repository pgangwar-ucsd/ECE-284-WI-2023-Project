## Note: Might need a supercomputer cluster to run this code as it requires close to 210 GB memory. This is the reason placement of reads could not complete on AWS (tried on t3.micro, t2.small).

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
