#!/bin/bash
max_procs=10
count=0

calculate(){
    #if /programs/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r $1 -A None -U None -o ${2}.pdbqt 2>&1 | grep -q 'WARNING\|Sorry'; then

    #if /home/yao/programs_src/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r $1 -A None -U None -o ${2}.pdbqt 2>&1 | grep -q 'WARNING\|Sorry'; then
        #echo "File wrong"
        #mv ${2}.pdbqt error_${2}.pdbq
    #else
        #echo "File correct"
    #fi 
    /home/yao/programs_src/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r $1 -A None -U None -o ${2}.pdbqt 2>&1

}
for a in 14  15  16  1b  2  21  25  26  27  28  29  30  50;do
    cd $a
    pwd
    for i in *.pdb; do
	if [[ ${i} != *"renamed"* ]]; then
            let "count = count + 1"
            echo "processing $i, file $count"
            x=$(echo $i | sed 's/.pdb//')
            #calculate $i $x &
            calculate $i $x &
            num_procs=$(ps --no-headers | wc -l)
            echo "num procs: $num_procs"
            if [ ${num_procs} -ge ${max_procs} ];then
	        wait -n 
            fi
	fi
    done
    cd ..
done
