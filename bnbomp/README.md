# OpenMP BnB interval solver

Coded by *A. Gorchakov*

The program use OpenMP to parallelize branch-and-bound methods based on simple interval bounds

**Usage:**

    bnbomp.exe name_of_bench knrec|unknrec eps max_steps omp_thread_number

**or (to list available test functions):**
    
    bnbomp.exe  list

Generic options:

Parameter | Description
------------ | -------------
"name_of_bench" | Name of the benchmark
knrec | preset the value of the record equal to the known optimum
uknrec | don't preset the value of the record 
eps | Precision
max_steps | The maximal number of steps to perform

Algorithm specific options:

Parameter | Description
------------ | -------------
omp_thread_number | number of OpenMP threads


