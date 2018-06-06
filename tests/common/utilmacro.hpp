#ifndef __UTILMACRO_HPP__
#define __UTILMACRO_HPP__ 

/**
  * Useful macros
  */

#define CP(rank, str, del) {char name[128]; FILE* f; sprintf(name, "bnb.%d.log", rank);f = fopen(name, "a"); fprintf(f, "%s\n", str); fflush(f); fclose(f); sleep(del);}

#define SGMAX(a,b) (((a)>(b))?(a):(b))

#define SGMIN(a,b) (((a)<(b))?(a):(b))

#define SGABS(a) SGMAX((a), (-(a)))

#define SGSQR(a) ((a) * (a))

#define SGSGN(a) (((a) >= 0) ? 1 : -1)

#define SGBOOLFLOOR(a, T) ((a < 0) ? std::numeric_limits<T>::min() : ((a < 1) ? 0 : 1))

#define SGBOOLCEIL(a, T) ((a > 1) ? std::numeric_limits<T>::max() : ((a > 0) ? 1 : 0))

#endif
