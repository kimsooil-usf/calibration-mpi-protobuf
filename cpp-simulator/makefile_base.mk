#this should be 'yes' for the status of the git repository to be checked during compiling

enable_proto=yes

#Set this to -DTIMING to enable timing output
timing = -DTIMING

#Set this to -DDEBUG to enable debugging output
debug =

#Set this to -DMERSENNE_TWISTER to use the Mersenne twister 19937 Random number generator
random = -DMERSENNE_TWISTER


ifeq ($(enable_proto), yes)
#set proto flags
LDLIBS = -lprotobuf
proto_def =-DENABLE_PROTO
obj = cohorts.o train_loader.o agents_store.pb.o initializers.o models.o interventions.o intervention_primitives.o updates.o simulator.o testing.o outputs.o drive_simulator.o
else
LDLIBS =
proto_def =
obj = cohorts.o train_loader.o initializers.o models.o interventions.o intervention_primitives.o updates.o simulator.o testing.o outputs.o drive_simulator.o
endif

include_paths = -Ilibs/ -Ilibs/cxxopts-2.2.0/include/
DEPFLAGS = -MMD -MP -MF $*.d
CXX = g++
CPPFLAGS = -Wall --std=c++14 -Ofast $(DEPFLAGS) $(include_paths) $(parallel) $(proto_def) $(timing) $(debug) $(random)

all: drive_simulator check

drive_simulator: $(obj)
	$(CXX) $(CPPFLAGS) $^ -o $@ $(LDLIBS)

%.o : $.cc %.d
	$(CXX) $(CPPFLAGS) -c $<

$(DEPDIR):
	mkdir -p $@


DEPFILES := $(obj:%.o=%.d)
$(DEPFILES):


.PHONY: clean
clean:
	rm drive_simulator *.o

.PHONY: check
check:
ifeq ($(check_git),yes)
ifeq ($(GIT_TREE_STATE),dirty)
	$(error git state is not clean)
endif
endif

include $(wildcard $(DEPFILES))