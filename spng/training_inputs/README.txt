# To run: 
snakemake --software-deployment-method apptainer --cores 12 \
    --apptainer-args=" -B /cvmfs,/home,/nfs,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid" \
    {targets} --directory={work_directory}