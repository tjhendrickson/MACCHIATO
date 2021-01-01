Bootstrap: docker

From: ubuntu:18.04

%files
run.py /run.py
tools /tools
workflow /workflow


%environment
export LC_ALL=C
export CARET7DIR=/opt/workbench/bin_rh_linux64
export OS=Linux
export FS_OVERRIDE=0
export FIX_VERTEX_AREA=
export FSF_OUTPUT_FORMAT=nii.gz
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
export PYTHONPATH=""

%post
# Make local folders/files
mkdir /share
mkdir /scratch
mkdir /local-scratch
mkdir /input_dir
mkdir /output_dir
mkdir /fsf_template_dir
touch /parcel_dlabel.nii

# Install basic utilities
apt-get -qq update
apt-get install -yq --no-install-recommends libquadmath0 libglib2.0-0 python wget bc bzip2 ca-certificates libgomp1 perl-modules tar tcsh unzip git libgomp1 perl-modules curl libgl1-mesa-dev libfreetype6 libfreetype6-dev nano


# Install miniconda3 and needed python tools
cd /opt
wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O /opt/Miniconda3.sh
bash /opt/Miniconda3.sh -b -p /opt/Miniconda3
export PATH="/opt/Miniconda3/bin:${PATH}"
. /opt/Miniconda3/etc/profile.d/conda.sh
conda activate 
conda install -y scipy numpy pandas mpi4py
conda install -c conda-forge -y statsmodels
conda install -c conda-forge -y "h5py>=2.9=mpi*"


pip install nibabel cifti PyWavelets nilearn sklearn git+git://github.com/aestrivex/bctpy.git@5f19d5aa9d14bf638ae6baf1b25280cf1222a476 
pip3 install argunparse

# Install the validator 0.26.11, along with pybids 0.6.5
apt-get update
apt-get install -y curl
curl -sL https://deb.nodesource.com/setup_10.x | bash -
apt-get remove -y curl
apt-get install -y nodejs
npm install -g bids-validator@0.26.11
pip install git+https://github.com/INCF/pybids.git@0.6.5

# Install entropy toolbox: https://raphaelvallat.com/entropy/build/html/index.html

pip install git+git://github.com/raphaelvallat/entropy.git@3bf3b6e937687e77965ecb5bfdc69e2a0f05936

# Install Connectome Workbench version 1.3.2
apt-get update
cd /opt
wget http://brainvis.wustl.edu/workbench/workbench-rh_linux64-v1.3.2.zip
unzip workbench-rh_linux64-v1.3.2.zip
export PATH=/opt/workbench/bin_rh_linux64:${PATH}

# Make scripts executable
chmod +x -R /run.py /tools /workflow

%runscript
. /opt/Miniconda3/etc/profile.d/conda.sh
conda activate
exec /run.py "$@"

