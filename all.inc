#General make file settings
CC = cc
CPP = c++
LIBTOOL = ar -r

#platform-specific options
#Compaq compiler C++ options:
#PCPPOPTS=-D__USE_STD_IOSTREAM -tweak

#Pthread linking support options
PTHREADLD=-pthread

#OpenMP specific compiler options
OPENMPOPT=-fopenmp

#OpenMP specific linker options
OPENMPLD=-fopenmp

#optimization options
#OPTIMOPT=-g
OPTIMOPT=-O


#Static linking option (if needed)
#STATIC=-static

#modern C++ support on
STD_OPT = -std=c++17

TESTLIB = $(ROOT)/tests

EXTRA_INC = -I$(TESTLIB) 
 
#common options
COMOPTS = $(STD_OPT)\
	   -I$(ROOT)\
           $(EXTRA_INC)\
       $(OPENMPOPT)\
       $(OPTIMOPT)\
       $(POLLOPT)

#C options
COPTS=$(COMOPTS) $(ADDIT_INCS) $(PCOPTS)

#C++ options
CPPOPTS= $(COMOPTS) $(ADDIT_INCS) $(PCPPOPTS)

#Linker options for static linking
LDFLAGS = -lm -ldl $(OPENMPLD) $(PTHREADLD) $(STATIC)


MFLAGS = MPICC=$(MPICC) COPTS="$(COPTS)" CPPOPTS="$(CPPOPTS) LIBTOOL="$(LIBTOOL)"

all: $(BINS) $(LIBS)

tests:
	@for i in $(TESTS); do ./$$i; continue; done

# Uncomment this if you want testing finishes after first finished test
#	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; return 1; fi done
	

$(BINS): $(LIBS) $(EXTLIBS)

$(LIBS): $(OBJS)
	$(LIBTOOL) $@ $(OBJS)

indent: 
	indent *.cpp ; indent *.hpp
clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CC) $(COPTS) -MM -c *.c > deps.inc;\
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true

.o.exe: $(GTEST_LIB)/*.a
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe
