#!/bin/bash

module use /tungstenfs/groups/gbioinfo/Appz/easybuild/modules/all
module use /tungstenfs/groups/gbioinfo/Appz/modules
module load dimsum

DiMSum --experimentDesignPath 000-data/experimentDesign_Jun.txt --countPath 000-data/dimsum_input.txt --projectName 004-dimsum_output --startStage 4 --numCores 8 --sequenceType noncoding --maxSubstitutions 202 --wildtypeSequence CAGGAGCGGATCAAGGCGGAGAGGAAGCGCATGAGGAACCGCATCGCTGCCTCCAAGTGCCGAAAAAGGAAGCTGGAGAGAATCGCCCGGCTGGAGGAAAAAGTGAAAACGTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTTAAACAGAAAGTCATGAACTAACCTA 1>004-dimsum.out 2>004-dimsum.err
