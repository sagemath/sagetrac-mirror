#!/usr/bin/env bash
echo "========================somewhatfptUFVS========================="
cd src




make
if [ $? -ne 0 ]; then
    echo >&2 "Error building the somewhatfptUFVS test program fptUFVS"
    exit 1
fi



make dotreader
if [ $? -e 0 ]; then
    echo "dot2neighbour_matrix converter compile. Will now procede with many easy problems."


else
    echo "Could not compile the dot2neighbour_matrix converter. Will now procede with even harder problems. Btw, no Boost Liberary?"
fi




echo "Doing NP-hard work!"
#time cat ../self_check_data_for_unweighted/real/arvid.neighbour_matrix | ./b.out
#time cat ./self_check_data_for_unweighted/real/16/arvid.neighbour_matrix | ./fptUFVS/b.out




make test
#make testthreaded	#OpenMP


#ERRORS=`grep ^##### $LOG`
#if [[ ! -z "$ERRORS" ]]; then
#    echo >&2 "Error there appeare to be at least one Wrong Answer when testing the test program somewhatfptUFVS. Please give feedback about this grave situation."
#    exit 1
#fi


if [ $? -ne 0 ]; then
#    echo >&2 "Error executing test code."
    echo >&2 "Error executing test."
    exit 1
fi

echo "The somewhatfptUFVS test completed successfully."








echo "================================================================"
