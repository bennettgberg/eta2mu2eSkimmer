for letter in C D E F G 
do
    for num in `seq 0 7`
    do  
        echo "letter $letter number $num"
        sed "s/letter=C/letter=$letter/g" sub_template.jdl >temp.jdl
        sed "s/setnum=0/setnum=$num/g" temp.jdl >sub_jobs.jdl
        rm temp.jdl
        condor_submit sub_jobs.jdl
    done
done
