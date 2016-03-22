RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: $(SIMDIR)
	(cd $(SIMDIR) ; $(RMPCDMD)/build/single_dimer_chemotactic_cell $(EXPDIR)/chemo_cell.parameters)

code:
	make --no-print-directory -C $(RMPCDMD)/build

$(SIMDIR):
	@mkdir $(SIMDIR)
