for csv in $(ls ../Data/*.csv|head -n 2); do
    arg1=$(echo $csv |perl -pe 's/\.\.\/Data\/PCs_//g'|cut -f 1 -d "_")
    arg2=$(echo $csv |perl -pe 's/\.\.\/Data\/PCs_//g'|cut -f 1 -d ".")
    n_col=$(head -1 "$csv" | awk -F',' '{print NF}')
    for col in $(seq 3 $n_col); do
        pc=$((col - 2))
        sbatch ./master_script.sh $arg1 bF_S1_PC${pc} $csv $col
    done
done