#### Common variables

source_dir := src/
include_dir := interface/
objs_dir := objs/

CXXFLAGS := -std=c++14 -O3 -g -Wall $(shell root-config --cflags) $(shell lhapdf-config --cflags) -I$(include_dir)
LDFLAGS := $(shell root-config --libs) -lGenVector $(shell lhapdf-config --ldflags) -lcuba -lDelphes
CXX := g++

_common_objs := binnedTF.o jacobianD.o MEEvent.o MEWeight.o transferFunction.o utils.o
common_objs := $(patsubst %,$(objs_dir)/%,$(_common_objs))
_common_deps := binnedTF.h jacobianD.h MEEvent.h MEWeight.h transferFunction.h utils.h
common_deps := $(patsubst %,$(include_dir)/%,$(common_deps))

#### TTbar specific variables

ttbar_proc_dir := MatrixElements/cppmem_pp_ttx_epmum/
ttbar_proc_obj := $(ttbar_proc_dir)/SubProcesses/P1_Sigma_sm_gg_epvebmumvmxbx/cppmem_pp_ttx_epmum.o

ttbar_dir := ttbar/
ttbar_exec := $(ttbar_dir)/ME_ttbar
_ttbar_objs := Integrand_TTbar.o ME_ttbar_main.o
ttbar_objs := $(patsubst %,$(ttbar_dir)/%,$(_ttbar_objs))

##### Common targets

all: ttbar

$(objs_dir)/%.o: $(source_dir)/%.cpp $(common_deps) | $(objs_dir)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(objs_dir):
	if [ ! -d $(objs_dir) ]; then mkdir $(objs_dir); fi

#### TTbar specific targets

ttbar: $(ttbar_exec)

$(ttbar_dir)/%.o: $(ttbar_dir)/%.cpp $(common_deps)
	$(CXX) -c $< -o $@ $(CXXFLAGS) -I$(ttbar_proc_dir)

$(ttbar_exec): $(common_objs) $(ttbar_objs) $(ttbar_proc_obj)
	$(CXX) -o $(ttbar_exec) $^ $(LDFLAGS) -L$(ttbar_proc_dir)/lib/ -lmodel_sm

#### Clean targets

.PHONY: clean clean_ttbar

clean: clean_ttbar
	-rm $(objs_dir)/*.o
	-if [ -d $(objs_dir) -a ! "$(ls -A $(objs_dir))" ]; then rmdir $(objs_dir); fi

clean_ttbar:
	-rm $(ttbar_dir)/*.o
	-if [ -e $(ttbar_exec) ]; then rm $(ttbar_exec); fi
