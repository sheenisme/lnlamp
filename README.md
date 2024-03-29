We presents a holistic approach–`PrecTuner`–by closely coupling the code generator and the autotuner via only one parameter $r$ . Initialized by automatically sampled values, $r$ is first used to generate several code variants in the polyhedral model, combining this optimization with loop transformations. These code variants are then used to solve a performance model expressed in terms of $r$ , possibly under a quality degradation budget. The value of $r$ that produces the best-performing mixed-precision code is finally predicted without evaluating all code variants. 

**The usage of `PrecTuner` is as follows:**
---

```
Usage: lnlamp file [options] ...

Options:
  -s|--specify-output   specify the target file
  -r|--rate-array       specify the rate-array of profiling (default is:0 12 25 37 50 62 75 87 99)
  -e|--error-threshold  specify the error-threshold (default is:99999999999999999999)
  -v|--version          lnlamp version
  -t|--tile             lnlamp specify the tile-sizes(default none) in PPCG
  -o|--openmp           lnlamp specify uses the --openmp in PPCG, and multi-threads(default is:4) running.
  -a|--sched-algo       lnlamp specify the schedule-algorithm to isl|feautrier|no in PPCG, default is isl
  -i|--isl-flags        lnlamp specify the flags in isl library of PPCG,default is null

Examples:
  lnlamp hello.c -s hello_lnlamp.c

Template:
  lnlamp <source-file> [Options] <target-file>
```

It is worth mentioning that all the experimental data in this article were obtained using this tool. In order to better test our method, we embedded the experimental test sub-module (polybench_benchmark) in this tool and encapsulated the testing process into a corresponding script so that our peers can quickly reproduce our experimental process.

**Install the `PrecTuner` separately** :
---

```shell
# prepare
sudo apt update && sudo apt upgrade -y
sudo apt-get install gcc g++ git vim make bc python python3-pip
sudo apt install automake autoconf libtool pkg-config libgmp3-dev libyaml-dev libclang-dev llvm clang
pip install pandas numpy matplotlib

cd /home/sheen/
mkdir lnlamp-install
git clone https://github.com/sheenisme/lnlamp.git

cd lnlamp/
./get_submodules.sh 
./autogen.sh 
./configure --prefix=/home/sheen/lnlamp-install

make
make install
```

For further installation help please refer to: https://repo.or.cz/ppcg.git.

**`Docker`** :
---
We provide a \texttt{dockerfile}(\url{https://github.com/sheenisme/lnlamp/blob/master/Dockerfile}) file for quick installation, but this installation is not recommended given the instability of performance testing in virtual machines.

It is important to note that we \textbf{highly recommend using the Dockerfile to install and use our environment}, even if you want to experiment on a local machine(the operating system we recommend is Ubuntu 20.04.5 LTS - Linux 5.15.0-67-generic) to get stable performance data, the installation steps are recommended to follow the steps described in the dockerfile, and if you want to use a custom path, it is currently recommended to use a soft link(`\texttt{ln -s /your\_dir /home/sheen}`) for the purpose.
We provide a `Dockerfile`([https://github.com/sheenisme/lnlamp/blob/master/Dockerfile](https://github.com/sheenisme/lnlamp/blob/master/Dockerfile)) file for quick installation.

**Reproducing our experiments** :
---

System configuration:
```
  Cpu   : Intel i5 processor (recommended)
  Gpu   : NVIDIA RTX 2070    (recommended)
  Os    : Ubuntu 20.04.5 LTS (recommended)
  Gcc   : 9.4.0  (recommended)
  Clang : 12.0.1 (recommended)
  Memory: Recommended 32G and above
  Disk  : Greater than 256G (artifact need about 124G)
```
As for the computational artifact of this work([https://github.com/sheenisme/lnlamp.git](https://github.com/sheenisme/lnlamp.git)), it implements automatic code generation based on PPCG([http://ppcg.gforge.inria.fr/](http://ppcg.gforge.inria.fr/)), and on Python and Shell scripts for automatic tuning, which closely coupling the code generator and the autotuner via only one parameter $r$ .

So, To reproduce the data in the paper, see the code repository-polybench_benchmark(https://github.com/sheenisme/polybench_benchmark.git)

And, the Artifact is publicly available on 10.5281/zenodo.10275776: https://zenodo.org/records/10275776. 

Note: The default installation path for all tools in this version of artifact is under '/home/sheen' and there is hard code for this path in the script code, we recommend soft linking your custom folder to '/home/sheen' so you can use it quickly and correctly. We also recommend that all users try installing with a dockerfile first, so that even if you have to install locally, you can refer to the docker image for help with some issues.

---
---
---

**PPCG descriptions**:
---

*The following is the ppcg author's description:*
Requirements:

- automake, autoconf, libtool
	(not needed when compiling a release)
- pkg-config (http://www.freedesktop.org/wiki/Software/pkg-config)
	(not needed when compiling a release using the included isl and pet)
- gmp (http://gmplib.org/)
- libyaml (http://pyyaml.org/wiki/LibYAML)
	(only needed if you want to compile the pet executable)
- LLVM/clang libraries, 2.9 or higher (http://clang.llvm.org/get_started.html)
	Unless you have some other reasons for wanting to use the svn version,
	it is best to install the latest supported release.
	For more details, including the latest supported release,
	see pet/README.

If you are installing on Ubuntu, then you can install the following packages:

automake autoconf libtool pkg-config libgmp3-dev libyaml-dev libclang-dev llvm

Note that you need at least version 3.2 of libclang-dev (ubuntu raring).
Older versions of this package did not include the required libraries.
If you are using an older version of ubuntu, then you need to compile and
install LLVM/clang from source.


Preparing:

Grab the latest release and extract it or get the source from
the git repository as follows.  This process requires autoconf,
automake, libtool and pkg-config.

	git clone git://repo.or.cz/ppcg.git
	cd ppcg
	./get_submodules.sh
	./autogen.sh


Compilation:

	./configure
	make
	make check

If you have installed any of the required libraries in a non-standard
location, then you may need to use the --with-gmp-prefix,
--with-libyaml-prefix and/or --with-clang-prefix options
when calling "./configure".


Using PPCG to generate CUDA or OpenCL code

To convert a fragment of a C program to CUDA, insert a line containing

	#pragma scop

before the fragment and add a line containing

	#pragma endscop

after the fragment.  To generate CUDA code run
	
	ppcg --target=cuda file.c

where file.c is the file containing the fragment.  The generated
code is stored in file_host.cu and file_kernel.cu.

To generate OpenCL code run

	ppcg --target=opencl file.c

where file.c is the file containing the fragment.  The generated code
is stored in file_host.c and file_kernel.cl.


Specifying tile, grid and block sizes

The iterations space tile size, grid size and block size can
be specified using the --sizes option.  The argument is a union map
in isl notation mapping kernels identified by their sequence number
in a "kernel" space to singleton sets in the "tile", "grid" and "block"
spaces.  The sizes are specified outermost to innermost.

The dimension of the "tile" space indicates the (maximal) number of loop
dimensions to tile.  The elements of the single integer tuple
specify the tile sizes in each dimension.
In case of hybrid tiling, the first element is half the size of
the tile in the time (sequential) dimension.  The second element
specifies the number of elements in the base of the hexagon.
The remaining elements specify the tile sizes in the remaining space
dimensions.

The dimension of the "grid" space indicates the (maximal) number of block
dimensions in the grid.  The elements of the single integer tuple
specify the number of blocks in each dimension.

The dimension of the "block" space indicates the (maximal) number of thread
dimensions in the grid.  The elements of the single integer tuple
specify the number of threads in each dimension.

For example,

    { kernel[0] -> tile[64,64]; kernel[i] -> block[16] : i != 4 }

specifies that in kernel 0, two loops should be tiled with a tile
size of 64 in both dimensions and that all kernels except kernel 4
should be run using a block of 16 threads.

Since PPCG performs some scheduling, it can be difficult to predict
what exactly will end up in a kernel.  If you want to specify
tile, grid or block sizes, you may want to run PPCG first with the defaults,
examine the kernels and then run PPCG again with the desired sizes.
Instead of examining the kernels, you can also specify the option
--dump-sizes on the first run to obtain the effectively used default sizes.


Compiling the generated CUDA code with nvcc

To get optimal performance from nvcc, it is important to choose --arch
according to your target GPU.  Specifically, use the flag "--arch sm_20"
for fermi, "--arch sm_30" for GK10x Kepler and "--arch sm_35" for
GK110 Kepler.  We discourage the use of older cards as we have seen
correctness issues with compilation for older architectures.
Note that in the absence of any --arch flag, nvcc defaults to
"--arch sm_13". This will not only be slower, but can also cause
correctness issues.
If you want to obtain results that are identical to those obtained
by the original code, then you may need to disable some optimizations
by passing the "--fmad=false" option.


Compiling the generated OpenCL code with gcc

To compile the host code you need to link against the file
ocl_utilities.c which contains utility functions used by the generated
OpenCL host code.  To compile the host code with gcc, run

  gcc -std=c99 file_host.c ocl_utilities.c -lOpenCL

Note that we have experienced the generated OpenCL code freezing
on some inputs (e.g., the PolyBench symm benchmark) when using
at least some version of the Nvidia OpenCL library, while the
corresponding CUDA code runs fine.
We have experienced no such freezes when using AMD, ARM or Intel
OpenCL libraries.

By default, the compiled executable will need the _kernel.cl file at
run time.  Alternatively, the option --opencl-embed-kernel-code may be
given to place the kernel code in a string literal.  The kernel code is
then compiled into the host binary, such that the _kernel.cl file is no
longer needed at run time.  Any kernel include files, in particular
those supplied using --opencl-include-file, will still be required at
run time.


Function calls

Function calls inside the analyzed fragment are reproduced
in the CUDA or OpenCL code, but for now it is left to the user
to make sure that the functions that are being called are
available from the generated kernels.

In the case of OpenCL code, the --opencl-include-file option
may be used to specify one or more files to be #include'd
from the generated code.  These files may then contain
the definitions of the functions being called from the
program fragment.  If the pathnames of the included files
are relative to the current directory, then you may need
to additionally specify the --opencl-compiler-options=-I.
to make sure that the files can be found by the OpenCL compiler.
The included files may contain definitions of types used by the
generated kernels.  By default, PPCG generates definitions for
types as needed, but these definitions may collide with those in
the included files, as PPCG does not consider the contents of the
included files.  The --no-opencl-print-kernel-types will prevent
PPCG from generating type definitions.


GNU extensions

By default, PPCG may print out macro definitions that involve
GNU extensions such as __typeof__ and statement expressions.
Some compilers may not support these extensions.
In particular, OpenCL 1.2 beignet 1.1.1 (git-6de6918)
has been reported not to support __typeof__.
The use of these extensions can be turned off with the
--no-allow-gnu-extensions option.


Processing PolyBench

When processing a PolyBench/C 3.2 benchmark, you should always specify
-DPOLYBENCH_USE_C99_PROTO on the ppcg command line.  Otherwise, the source
files are inconsistent, having fixed size arrays but parametrically
bounded loops iterating over them.
However, you should not specify this define when compiling
the PPCG generated code using nvcc since CUDA does not support VLAs.


CUDA and function overloading

While CUDA supports function overloading based on the arguments types,
no such function overloading exists in the input language C.  Since PPCG
simply prints out the same function name as in the original code, this
may result in a different function being called based on the types
of the arguments.  For example, if the original code contains a call
to the function sqrt() with a float argument, then the argument will
be promoted to a double and the sqrt() function will be called.
In the transformed (CUDA) code, however, overloading will cause the
function sqrtf() to be called.  Until this issue has been resolved in PPCG,
we recommend that users either explicitly call the function sqrtf() or
explicitly cast the argument to double in the input code.


Contact

For bug reports, feature requests and questions,
contact http://groups.google.com/group/isl-development

Whenever you report a bug, please mention the exact version of PPCG
that you are using (output of "./ppcg --version").  If you are unable
to compile PPCG, then report the git version (output of "git describe")
or the version number included in the name of the tarball.


Citing PPCG

If you use PPCG for your research, you are invited to cite
the following paper.

@article{Verdoolaege2013PPCG,
    author = {Verdoolaege, Sven and Juega, Juan Carlos and Cohen, Albert and
		G\'{o}mez, Jos{\'e} Ignacio and Tenllado, Christian and
		Catthoor, Francky},
    title = {Polyhedral parallel code generation for CUDA},
    journal = {ACM Trans. Archit. Code Optim.},
    issue_date = {January 2013},
    volume = {9},
    number = {4},
    month = jan,
    year = {2013},
    issn = {1544-3566},
    pages = {54:1--54:23},
    doi = {10.1145/2400682.2400713},
    acmid = {2400713},
    publisher = {ACM},
    address = {New York, NY, USA},
}
