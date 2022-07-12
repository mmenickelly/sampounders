conc=30


tmux new-session -d 

for ((i=1; i<="$((conc))"; i++))
do 
  tmux split-window -h
  tmux send-keys -t $i "matlab -r \"run_experiment3 $((i))\" -nodesktop -nosplash -nodisplay -nojvm" C-m
  tmux select-layout tile
done

tmux -2 attach-session -d 
