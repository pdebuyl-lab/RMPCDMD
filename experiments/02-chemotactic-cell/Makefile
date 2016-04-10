RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-intro help
include ../Makefile


simulation: $(SIMDIR)
	(cd $(SIMDIR) ; $(RMPCDMD)/build/single_dimer_chemotactic_cell $(EXPDIR)/chemo_cell.parameters)

$(SIMDIR):
	@mkdir $(SIMDIR)

