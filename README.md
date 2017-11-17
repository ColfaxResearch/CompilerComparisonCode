# CompilerComparisonCode
C++ Computational Kernels for "A PERFORMANCE-BASED COMPARISON OF C/C++ COMPILERS"


Build
-------
Build instructions are provided for Linux machines. This project uses the SCons build system (http://scons.org/). Please install SCons 3 or higher following the instructions on the SCons homepage.

To build a given compute kernel `k', you may issue the following commands - 

  `bash-prompt$ cd c++/k` to cd into the kernel folder.
  
  `bash-prompt$ scons archirecture=skl algorithm=ALG compiler=COMP` to build the kernel. Here `ALG` is the name of the algorithm and can be found in the primary source file for the kernel.
  
  For example, for the Structure Function kernel, set `ALG` to `OBLKIO_OPT`.
  
  Finally, `COMP` is the compiler that you wish to compile with. To build with the Intel compiler, set `COMP` to `Intel` etc...
