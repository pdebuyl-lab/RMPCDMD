RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-help help
include ../Makefile

simulation: $(SIMDIR)
	(cd $(SIMDIR) ; $(RMPCDMD)/build/single_janus_pbc $(EXPDIR)/janus.parameters)

ruckner-kapral.h5: ruckner-kapral.parameters $(RMPCDMD)/build/single_dimer_pbc
	$(RMPCDMD)/build/single_dimer_pbc $<

$(SIMDIR):
	@mkdir $(SIMDIR)

