FROM ubuntu:20.04



# install dependency
ENV DEBIAN_FRONTEND noninteractive
ENV TZ 'Asia/Shanghai' 
ENV LANG en_US.UTF-8 
ENV LANGUAGE en_US:en 
ENV LC_ALL en_US.UTF-8
RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list
RUN sed -i s@/security.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list
RUN apt-get update \
    && apt-get install -y gcc g++ git vim make python python3-pip flex bison texinfo \
    && apt install -y automake autoconf libtool pkg-config libgmp3-dev libyaml-dev libclang-dev llvm clang cmake
RUN pip install pandas numpy matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple



# setup user account 
ARG user=sheen
RUN useradd --create-home --no-log-init --shell /bin/bash ${user} \
    && adduser ${user} sudo   \
    && echo "${user}:1" | chpasswd
# RUN usermod -u 1000 ${user} && usermod -G 1000 ${user}
USER ${user}



# download package
WORKDIR /home/${user}
RUN git clone https://github.com/sheenisme/lnlamp.git
RUN git clone https://github.com/sheenisme/llvm-project.git
RUN git clone https://github.com/sheenisme/TAFFO.git
RUN git clone https://github.com/sheenisme/pluto.git



# install PrecTuner
RUN cd /home/sheen/          \
    && mkdir lnlamp-install  \
    && cd lnlamp/            \
    && ./get_submodules.sh   \
    && ./autogen.sh          \
    && ./configure --prefix=/home/sheen/lnlamp-install  \
    && make  \
    && make install

# install LLVM
RUN cd /home/sheen/llvm-project/  \
    && mkdir llvm-install         \
    && git checkout origin/release/12.x  \
    && cmake -S ./llvm -B llvm-build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE="Debug" -DLLVM_VERSION_MAJOR="12"  -DLLVM_TARGETS_TO_BUILD="X86" -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;compiler-rt;openmp;polly" -DLLVM_BUILD_LLVM_DYLIB=ON -DLLVM_LINK_LLVM_DYLIB=ON -DCMAKE_INSTALL_PREFIX=/home/sheen/llvm-project/llvm-install  \
    && cd llvm-build/  \
    && make -j2  \
    && make install

# install LUIS
RUN cd /home/sheen/TAFFO/   \
    && mkdir taffo-install  \
    && export LLVM_DIR=/home/sheen/llvm-project/llvm-install  \
    && cmake -S . -B build -DTAFFO_BUILD_ORTOOLS=ON -DCMAKE_INSTALL_PREFIX=/home/sheen/TAFFO/taffo-install  \
    && cd build/  \
    && make -j2   \
    && make install

# install Pluto
RUN cd /home/sheen/pluto/    \
    && mkdir pluto-install   \
    && git submodule init    \
    && git submodule update  \
    && ./autogen.sh  \
    && ./configure --prefix=/home/sheen/pluto/pluto-install  \
    && make  \
    && make install


# set up paths
ENV PATH /home/sheen/pluto/pluto-install/bin:$PATH
ENV PATH /home/sheen/TAFFO/taffo-install/bin:$PATH
ENV PATH /home/sheen/llvm-project/llvm-install/bin:$PATH
ENV PATH /home/sheen/lnlamp-install/bin:$PATH

ENV LLVM_DIR /home/sheen/llvm-project/llvm-install

ENV LD_LIBRARY_PATH /home/sheen/TAFFO/taffo-install/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /home/sheen/llvm-project/llvm-install/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /home/sheen/lnlamp-install/lib:$LD_LIBRARY_PATH

ENV CPATH /home/sheen/lnlamp/polybench_benchmark/utilities:$CPATH
ENV CPATH /usr/local/include:$CPATH



# launch shell
CMD ["/bin/bash"]
