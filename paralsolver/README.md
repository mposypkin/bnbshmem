# Multithreaded BnB interval solver

Coded by *Y. Yamchenko*

The program use C++ 17 threads to parallelize the branch-and-bound methods based on simple interval bounds with "Manager-Slave" paradigm.

**Usage:**

    parasol.exe <'name_of_bench'> <eps> <max_steps> <virtual_procs_number> <split_steps_limit> <split_subs_limit> <split_steps_coeff> <split_subs_coeff>

**or (to run all tests):**

    parasol.exe

**or (to list available test functions):**

    parasol.exe list

**or (to see help message):**

    parasol.exe --help

Generic options:

Parameter | Description
------------ | -------------
"name_of_bench" | Name of the benchmark
eps | Precision (default: 0,1)
max_steps | The maximal number of steps to perform (default: 10^6)

Algorithm specific options:

Parameter | Description
------------ | -------------
virtual_procs_number | number of virtual processors (default: 4)
split_steps_limit | limit of steps to split the problem (default: 1000)
split_subs_limit | limit of subproblems in bag to split the problem (default: 10)
split_steps_coeff | percent of steps to deliver to new working thread while splitting (default: 0,5)
split_subs_coeff | percent of subproblems to deliver to new working thread while splitting (default: 0,5)
