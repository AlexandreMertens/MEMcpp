LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -std=c++11 -O3 -fcx-fortran-rules -fcx-limited-range
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

all: cppmem_pp_ttx_epmum.o

cppmem_pp_ttx_epmum.o: cppmem_pp_ttx_epmum.cc $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LIBFLAGS) -I$(INCDIR)

$(LIBDIR)/lib$(MODELLIB).a:
	pushd ../../src && make && popd

.PHONY: clean

clean:
	rm -f *.o
	pushd ../../src && make clean && popd
