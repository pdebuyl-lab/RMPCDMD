RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)
SIMDIR:=$(PWD)

all: exp-intro help 01-help

include ../Makefile

simulation: dimer.h5

%.h5: %.parameters $(RMPCDMD)/build/single_dimer_pbc
	$(RMPCDMD)/build/rmpcdmd run single_dimer_pbc $< $@ auto

01-help:
	@echo "  ruckner-kapral.h5 : Perform a simulation with the parameters of the Ruckner-Kapral 2007 paper"
