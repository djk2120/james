#!/bin/bash

rm *~
rm all.bib
echo "" > all.tmp
for file in *.bib
do
    echo $file
    cat $file >> all.tmp
    echo "" >> all.tmp
done
mv all.tmp all.bib