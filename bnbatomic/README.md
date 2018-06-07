#Multithreaded BnB interval solver

Coded by *M. Posypkin*

The program use C++ 17 threads to parallelize the branch-and-bound methods based on simple interval bounds.

**Usage:**

    bnbatomic.exe "name_of_bench" eps max_steps virtual_procs_number parallel_steps_limit



Generic options:

Parameter | Description
------------ | -------------
 "name_of_bench" | Name of the benchmark
eps | Precision
max_steps | The maximal number of steps to perform

Algorithm specific options:

Parameter | Description
------------ | -------------
virtual_procs_number | number of virtual processors
parallel_steps_limit | the minimal number of steps to start parallelization


