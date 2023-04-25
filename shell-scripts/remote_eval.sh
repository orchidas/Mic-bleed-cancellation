 #!/bin/bash

# send the latest matlab script to the server
scp my\ code/{mle_closed_form_batch,map_closed_form_batch,debleed}.m \
orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/my\ code/"
scp my\ code/test/var_dist.m \
orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/my\ code/test/"

# finished sending scripts
echo "Sent latest scripts"

# tmux session name
# session="var_dist"

# login to the matlab machine on the CCRMA server
ssh orchi@cm-matlab.stanford.edu << ENDSSH

# go to the scratch folder where data is saved
cd /scratch/orchi/src/my\ code/test/

echo "Running MATLAB script"
# run matlab without GUI
nohup matlab -nodisplay -nodesktop -nojvm -nosplash < var_dist.m &

# close tmux session once done
# tmux kill-session -t $session

exit

ENDSSH

# print when done
echo "Finished running matlab script"

