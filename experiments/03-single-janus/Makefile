RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-help help
include ../Makefile

simulation: janus.h5

%.h5: %.parameters $(RMPCDMD)/build/single_body
	$(RMPCDMD)/build/rmpcdmd run single_body $< $@ auto

$(SIMDIR):
	@mkdir $(SIMDIR)

