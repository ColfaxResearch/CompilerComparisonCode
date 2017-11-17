#PBS -N pbsResult
#PBS -l nodes=1:plat8160
#PBS -l walltime=72:00:00
ARCH=skl
ARCH_THREADS=192
ARCH_THREADS_PGI=48
cd $HOME/code/intelShootout/python
python c++Runner.py sfComp OBLKIO_OPT $ARCH --compilers Intel GNU LLVM PGI AOCC ZAPCC --sizes 65536 --maxthreads 1 --postfix OBLKIO_VEC.$ARCH -r 3 -t 10 -p 10
python c++Runner.py sfComp OBLKIO_OPT $ARCH --compilers Intel GNU LLVM AOCC ZAPCC --sizes 65536 --maxthreads $ARCH_THREADS --postfix OBLKIO_OPT.$ARCH -r 3 -t 10 -p 10
python c++Runner.py sfComp OBLKIO_OPT $ARCH --compilers PGI --sizes 65536 --maxthreads $ARCH_THREADS_PGI --postfix OBLKIO_OPT.$ARCH -r 3 -t 10 -p 10
cd $HOME/code/intelShootout/python
python c++Runner.py luDecomp KIJ_OPT $ARCH --compilers Intel GNU LLVM AOCC ZAPCC PGI --sizes 256 --maxthreads 1 --postfix KIJ_VEC.$ARCH -r 3 -t 10 -p 100
python c++Runner.py luDecomp KIJ_OPT $ARCH --compilers Intel GNU LLVM AOCC ZAPCC --sizes 1024 --maxthreads $ARCH_THREADS --postfix KIJ_OPT.$ARCH -r 3 -t 10 -p 10
python c++Runner.py luDecomp KIJ_OPT $ARCH --compilers PGI --sizes 1024 --maxthreads $ARCH_THREADS_PGI --postfix KIJ_OPT.$ARCH -r 3 -t 10 -p 10
cd $HOME/code/intelShootout/python
python c++Runner.py jacobiSolve SOLVE_OPT $ARCH --compilers Intel GNU LLVM AOCC ZAPCC PGI --sizes 64 --maxthreads 1 --postfix SOLVE_VEC.$ARCH -r 3 -t 10 -p 1
python c++Runner.py jacobiSolve SOLVE_OPT $ARCH --compilers Intel GNU LLVM AOCC ZAPCC --sizes 256 --maxthreads $ARCH_THREADS --postfix SOLVE_OPT.$ARCH -r 3 -t 10 -p 1
python c++Runner.py jacobiSolve SOLVE_OPT $ARCH --compilers PGI --sizes 256 --maxthreads $ARCH_THREADS_PGI --postfix SOLVE_OPT.$ARCH -r 3 -t 10 -p 1
