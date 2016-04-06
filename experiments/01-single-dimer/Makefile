RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: help 01-help

include ../Makefile

simulation: $(SIMDIR)
	(cd $(SIMDIR) ; $(RMPCDMD)/build/single_dimer_pbc $(EXPDIR)/dimer.parameters)

code:
	make --no-print-directory -C $(RMPCDMD)/build

ruckner-kapral.h5: ruckner-kapral.parameters $(RMPCDMD)/build/single_dimer_pbc
	$(RMPCDMD)/build/single_dimer_pbc $<

$(SIMDIR):
	@mkdir $(SIMDIR)

01-help:
	@echo "  ruckner-kapral.h5 : Perform a simulation with the parameters of the Ruckner-Kapral 2007 paper"
