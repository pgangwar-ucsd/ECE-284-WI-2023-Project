export PATH=$PATH:$PWD/build
cd build
make -j
cd ../
matUtils place_read -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.1.198 -d 1 -v my_vcf -r 150 -w 20 -e 0 -s 2 -f test/NC_045512v2.fa
