# To run: 
snakemake \
    --cores {NCORES} \
    --software-deployment-method apptainer \ #This is needed to run larsoft
    --apptainer-args=" -B /cvmfs,/home,/nfs,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid" \ #This is needed to run larsoft
    {targets} --directory={work_directory}
