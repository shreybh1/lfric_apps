# runscript to run gungho model with different layers and nodes. 
export TIME="0:5:00"
export PPN=32
for MESH_SIZE in 128
do 
  export NUM_NODES=6
  for JOBS in 0 # 1 2 
  do
    echo "Num nodes="$NUM_NODES " TPN="$PPN " time="$TIME  " ARRAY="$JOBS 
    sbatch  --nodes=${NUM_NODES} --tasks-per-node=${PPN} --time=${TIME} --array=${JOBS} --wait pawsey.slurm
    wait
  done
done 
