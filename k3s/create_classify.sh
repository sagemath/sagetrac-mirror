# USAGE
# ./create_classify.sh p q pure pathtofolder [path to sage (optional); default is 'sage']
#
# The script must be executed inside the sage folder
#
# EXAMPLE
# ./create_classify.sh 2 3 False k3s/results/no0_order2/
#
# The script prints to output, so do
#
# ./create_classify.sh 2 3 False k3s/results/no0_order2/ >> commands_batch
#
# to have a file with one sage run per line
#
# If you use parallel you can then do
# parallel -j $(number of jobs) < commands_batch
# to start the jobs in parallel
#
# One can do
# cat *.result | grep "complete" | wc -l 
# To count the number of finished files

p=$1
q=$2
pure=$3
folder=$4
if [ -z "$5" ]; then
  sagepath="sage"
else
  sagepath=$5
fi

cd $folder;

for i in *; do
  input="$folder$i"
  output="$folder$i.$q.result"
  log="$folder$i.$q.log"
  echo $sagepath -c "'load(\"k3s/K3_aut_classification.sage\"); classify(p=$p,q=$q,file_read=\"$input\",file_write=\"$output\",pure=$pure,verbose=2)' >$log 2>&1"
done;
