squeue -u $USER -h -o "%i" | while read jobid; do
    echo "=============================="
    echo "JOB $jobid"
    scontrol show job "$jobid" | grep -E "JobId=|JobName=|JobState=|WorkDir=|Command="
done
