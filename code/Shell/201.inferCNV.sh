# create a function to infer CNV
inferCNV() {
    python -u tools/SequencingCancerFinder/infer.py --ckp=data/203.inferCNV/model_epoch92.pkl \
        --matrix=$1 \
        --out=data/203.inferCNV/result/$(basename $1 .tsv).csv
}
export -f inferCNV
# use parallel --progress --keep-order --line-buffer to run the function on all files in data/203.inferCNV/raw_count_matrix
ls data/203.inferCNV/raw_count_matrix/*.tsv | parallel --progress --keep-order --line-buffer -j10 inferCNV
