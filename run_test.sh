src=/Users/joazofeifa/Lab/fast-er_inteval_searcher/src/FIS
name=test
out=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_output_directory/
IN=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_input_directory/
query=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_query_file.bed
log_out=/Users/joazofeifa/Lab/fast-er_inteval_searcher/test_output_directory/
pad=0
upad=100
rebuild=1
pairwise=0
$src -N $name -i $IN -q $query -o $out  -rebuild $rebuild -pairwise $pairwise -log_out $log_out -pad $pad -upad $upad