# 
seqs=("WWWW" "WRWW")
straight=("0" "1")
antipara=("0" "1")


source activate py38
for seq in "${seqs[@]}"; do
    for st in "${straight[@]}"; do
        for ap in "${antipara[@]}"; do
            python 01-split_layer.py $seq $ap $st
            python 02-add_sidechain.py $seq $ap $st
            python 03-martinize.py $seq $ap $st
        done
    done
done
rm \#*
