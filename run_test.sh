src=/Users/joazofeifa/Lab/fast_interval_searcher/src/FIS
name=test
out=/Users/joazofeifa/Lab/fast_interval_searcher/test_output_directory/
IN=/Users/joazofeifa/Lab/fast_interval_searcher/test_input_directory/
query=/Users/joazofeifa/Lab/fast_interval_searcher/DMSO2_3-1_bidirectional_hits_intervals.bed
log_out=/Users/joazofeifa/Lab/fast_interval_searcher/test_output_directory/
rebuild=1
$src -N $name -i $IN -q $query -o $out  -rebuild $rebuild -log_out $log_out