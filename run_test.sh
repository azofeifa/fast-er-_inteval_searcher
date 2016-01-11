src=/Users/joazofeifa/Lab/fast-er_inteval_searcher/src/FIS
name=test
out=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_output_directory/
IN=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_input_directory/
query=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-19_divergent_classifications.bed
pad=0
upad=1500
pairwise=0
MIN=0
$src -N $name -i $IN -q $query -o $out  -pairwise $pairwise -pad $pad -upad $upad -min $MIN