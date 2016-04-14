RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-help help
include ../Makefile

simulation: $(SIMDIR)
	(cd $(SIMDIR) ; $(RMPCDMD)/build/single_dimer_channel_flow $(EXPDIR)/channel_flow.parameters)

$(SIMDIR):
	@mkdir $(SIMDIR)

