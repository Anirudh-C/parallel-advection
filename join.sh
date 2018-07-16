#!/bin/sh
for i in {0..99}
do
    cat "src/data/sol0-$i.dat" "src/data/sol1-$i.dat" "src/data/sol2-$i.dat" "src/data/sol3-$i.dat" >> "result/$1/$2/sol$i.dat"
done
