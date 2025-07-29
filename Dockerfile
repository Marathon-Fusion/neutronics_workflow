FROM python:3.11-bookworm

# Install system dependencies
RUN apt-get update && \
    apt-get install -y build-essential \
    g++ \
    cmake\
    libhdf5-dev\
    libpng-dev\ 
    git\
    libeigen3-dev \
    wget \
    autotools-dev\
    autoconf \
    # libnetcdf-dev is needed to allow NETCDF on MOAB which helps with tet meshes in OpenMC
    libnetcdf-dev \
    # libtbb-dev required for DAGMC
    libtbb-dev \
    # libglfw3-dev required for DAGMC
    libglfw3-dev \
    # needed for CadQuery functionality
    libgl1-mesa-glx \
    # needed for CadQuery functionality
    libgl1-mesa-dev \
    # needed for CadQuery functionality
    libglu1-mesa-dev \
    libglu1-mesa\
    # needed for CadQuery functionality
    freeglut3-dev \
    # needed for CadQuery functionality
    libosmesa6 \
    # needed for CadQuery functionality
    libosmesa6-dev \
    # needed for CadQuery functionality
    libgles2-mesa-dev \
    # needed for Gmsh functionality
    libxft2 \
    # needed for gmsh
    libxcursor-dev \
    # needed for gmsh
    libxinerama-dev \
    libexodusii-dev \
    libfmt-dev \
    libxi-dev \
    libxmu-dev



# Set up directories
RUN mkdir /opt/marathon
RUN mkdir /opt/marathon/dagmc_bld
WORKDIR /opt/marathon


# Install Python requirements
RUN python -m ensurepip --upgrade
COPY requirements.txt /opt/marathon/
RUN pip install -r requirements.txt


# Clone and build 
RUN mkdir -p /opt/marathon/dagmc_bld/MOAB/bld
WORKDIR /opt/marathon/dagmc_bld/MOAB
RUN git clone  --single-branch --branch 5.3.1 --depth 1 https://bitbucket.org/fathomteam/moab.git  
WORKDIR /opt/marathon/dagmc_bld/MOAB/bld
RUN cmake ../moab -DENABLE_HDF5=ON \
-DENABLE_NETCDF=ON \
-DENABLE_FORTRAN=OFF \
-DENABLE_BLASLAPACK=OFF \
-DBUILD_SHARED_LIBS=ON \
-DENABLE_PYMOAB=ON \
-DCMAKE_INSTALL_PREFIX=/opt/marathon/dagmc_bld/MOAB/ 
RUN make -j4
RUN make install


# Clone and build DAGMC
RUN mkdir /opt/marathon/dagmc_bld/DAG
WORKDIR /opt/marathon/dagmc_bld/DAG
RUN git clone https://github.com/svalinn/DAGMC.git
RUN mkdir /opt/marathon/dagmc_bld/DAG/bld
WORKDIR /opt/marathon/dagmc_bld/DAG/bld
RUN cmake ../DAGMC -DMOAB_DIR=/opt/marathon/dagmc_bld/MOAB/ \
-DBUILD_TALLY=ON \
-DCMAKE_INSTALL_PREFIX=/opt/marathon/dagmc_bld/DAG
RUN make
RUN make install


# Clone OpenMC repository
WORKDIR /opt/marathon
RUN git clone --recurse-submodules https://github.com/openmc-dev/openmc.git

# Build OpenMC with DAGMC support
RUN mkdir /opt/marathon/openmc/build
WORKDIR /opt/marathon/openmc/build
RUN cmake .. -DOPENMC_USE_DAGMC=on -DCMAKE_PREFIX_PATH=/opt/marathon/dagmc_bld/DAG/
RUN make
RUN make install

# Upgrade pip and install Python dependencies
RUN pip install --upgrade pip
WORKDIR /opt/marathon/openmc
RUN python -m pip install .


# Set up nuclear data
RUN mkdir /opt/marathon/openmc/nucleardata
WORKDIR /opt/marathon/openmc/nucleardata
RUN wget https://anl.box.com/shared/static/9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz && \
    tar -xvf 9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz -C /opt/marathon/openmc/nucleardata/ && \
    rm 9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz

# Set environment variable for cross sections
ENV OPENMC_CROSS_SECTIONS=/opt/marathon/openmc/nucleardata/endfb-vii.1-hdf5/cross_sections.xml


#Build pymoab because somehow this wasnt done already????
WORKDIR /opt/marathon/dagmc_bld/MOAB/bld/pymoab
RUN python3 setup.py install --user

WORKDIR /opt/marathon

#Install yaml-cpp for gmsh2exo
RUN git clone https://github.com/jbeder/yaml-cpp.git
RUN mkdir /opt/marathon/yaml-cpp/build
WORKDIR /opt/marathon/yaml-cpp/build
RUN cmake -DBUILD_SHARED_LIBS=ON ..
RUN make
RUN make install

WORKDIR /opt/marathon

#Install gmsh2exo
RUN git clone https://github.com/andrsd/gmsh2exo.git
RUN mkdir /opt/marathon/gmsh2exo/build
WORKDIR /opt/marathon/gmsh2exo/build

RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/bin ..
RUN make
RUN make install

RUN cp /opt/marathon/gmsh2exo/build/gmsh2exo /usr/local/bin

# Create workspace directory
RUN mkdir /opt/marathon/workspace
WORKDIR /opt/marathon/workspace

ENTRYPOINT ["python", "-u"]