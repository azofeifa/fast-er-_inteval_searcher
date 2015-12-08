src=/Users/joazofeifa/Lab/fast-er_inteval_searcher/src/FIS
name=test
out=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_output_directory/
IN=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_input_directory/
query=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_query_file.bed
pad=0
upad=1500
pairwise=0
$src -N $name -i $IN -q $query -o $out  -pairwise $pairwise -pad $pad -upad $upad