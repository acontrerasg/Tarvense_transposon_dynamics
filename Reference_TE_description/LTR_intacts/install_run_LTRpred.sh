# Create and activate conda environment, optional
conda create -n udocker_prod -c defaults python=2.7
conda activate udocker_prod

# Install most recent version from udocker github repo
pip install git+https://github.com/indigo-dc/udocker

# install udocker
udocker install
# Create a symlink to tmp folder so udocker has more dspace to extract and store containers.
ln -s  /tmp/global2/acontreras/.udocker/  ~/

# Prepare container
udocker pull drostlab/ltrpred
udocker create --name=ltrpred drostlab/ltrpred

# Required on some systems to run the container
export PROOT_NO_SECCOMP=1

# Run container  with mounted folder where data is stored as:
udocker run -v /ebio/abt6_projects8/Tarvense_TE/data/Tarvense_LTRpred:/app/ltrpred_data/ ltrpred
# 
cd ltrpred_data/
####### ADD VARIABLES
genome="/app/ltrpred_data/modified.fasta"
#####
Dfam="/app/ltrpred_data//Dfam_v3.1" # Folder containing Dfam database with the name Dfam. and extensions. 
Rscript  ./run_LTRpred.R  ${genome} ${Dfam}
