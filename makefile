CXXFLAGS := -O2 -g -Wall $(shell root-config --cflags) $(shell lhapdf-config --cflags)
LDFLAGS := -lm $(shell root-config --libs --glibs) $(shell lhapdf-config --ldflags) -lcuba -lDelphes
CXX := g++

sourcedir := examples/


ww_PROC_DIR :=/home/fynu/amertens/scratch/MatrixElement/MG5_aMC_v2_2_3/uu_ww_1d_cpp/
ww_MGproc := $(ww_PROC_DIR)/SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.o $(ww_PROC_DIR)/src/*.o


all: ww

ww: $(sourcedir)/ME_ww_CUBA.o $(sourcedir)/utils.o $(sourcedir)/jacobianF.o
	$(CXX) -o $(sourcedir)/ME_ww_CUBA.exe $(ww_MGproc) $^ $(LDFLAGS)

$(sourcedir)/ME_ww_CUBA.o: $(sourcedir)/ME_ww_CUBA.cpp
	$(CXX) $(CXXFLAGS) -I$(ww_PROC_DIR) -c $< -o $@

$(sourcedir)/jacobianF.o: $(sourcedir)/jacobianF.cpp $(sourcedir)/jacobianF.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(sourcedir)/utils.o: $(sourcedir)/utils.cpp $(sourcedir)/utils.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	-rm $(sourcedir)/*.exe $(sourcedir)/*.o
