# Madagascar_Propagations

These propagations were developed during a short Madagascar Course. FDTD were used for wave propagation over the 
models Marmousi (2D) and Overthrust (3D).


# ---------------------------------------------------------------------
# You will need Madagascar for this:

# UBUNTU:

sudo apt-get install scons libc6-dev gcc g++ gfortran libxaw7-dev \
freeglut3-dev libnetpbm10-dev libtiff4-dev libgd-dev libplplot-dev \
libavcodec-dev libcairo2-dev libjpeg-dev libblas-dev swig python-dev \
python-numpy libopenmpi-dev libfftw3-dev libsuitesparse-dev python-epydoc

# FEDORA-CENTOS:

yum -y install binutils gcc glibc-headers scons texlive-latex subversion \
gcc-c++ gcc-gfortran numpy python swig octave octave-devel libgomp openmpi \
openmpi-devel blas-devel atlas atlas-devel units gifsicle ffmpeg-devel \
libtiff-devel libjpeg-devel plplot-devel mesa-libGL-devel freeglut-devel \
libXaw-devel netpbm-devel

# Download from here:
 http://sourceforge.net/projects/rsf/files/

# From terminal type

 $ tar xvzf madagascar-1.7.tar.gz
 $ mv -u madagascar-1.7 $HOME/rsfsrc
 $ cd $HOME/rsfsrc
 $ export RSFROOT=$HOME/rsf
 $ export RSFSRC=$HOME/rsfsrc
 $ export DATAPATH=$HOME/Data/
 $ mkdir $HOME/Data
 $ ./configure
 $ make; make install
 $ cat $RSFSRC/env.sh >> $HOME/.bashrc
# ---------------------------------------------------------------------

# First of all untar the 2D and 3D models (from terminal type):

 $ tar -xzvf flujo_3/marmousi2/marmo2.tar.gz
 $ cat flujo_4/overthrust/over.tar.gz.part-a* | tar xz

# ---------------------------------------------------------------------

# Now execute in terminal:

# For 2D Propagation:

 $ cd flujo_3
 $ sftour scons
 $ scons
 $ scons view


# For 3D Propagation:

 $ cd flujo_4
 $ sftour scons
 $ scons
 $ scons view
# ---------------------------------------------------------------------


# ENJOY

