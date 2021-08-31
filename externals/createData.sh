#!/bin/bash

#Get Paths to binaries
echo Please input prodigal binary path:
read prodPath

echo Please input the fetchMG Perl script path:
read fetchmgPath

#Output paths for check
echo Prodigal: $prodPath
echo FetchMG: $fetchmgPath

echo Please input contig fasta filepath:
read contigPath

echo Name the prodigal result file:
read prodigalFilename

echo Name the fetchMG result Dictionary:
read fetchmgDictname

# Testing
echo
echo $contigPath 
echo $prodigalFilename
echo $fetchmgDictname
echo

./$prodPath -a $prodigalFilename -i $contigPath -p meta


./$fetchmgPath -d $contigPath -o $fetchmgDictname -m extraction $prodigalFilename

