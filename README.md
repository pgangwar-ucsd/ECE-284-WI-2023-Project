# ECE-284-WI-2023-Project

Steps to run:
1. git clone https://github.com/pgangwar-ucsd/usher.git
2. cd usher
3. rm -rf build
4. ./install/installUbuntu.sh
5. rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
6. chmod +x faToVcf
7. mv faToVcf scripts/
8. source run.sh
