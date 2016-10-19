RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-help help
include ../Makefile

simulation: janus.h5

%.h5: %.parameters $(RMPCDMD)/build/single_janus_pbc
	$(RMPCDMD)/build/rmpcdmd run single_janus_pbc $< $@ auto

$(SIMDIR):
	@mkdir $(SIMDIR)

