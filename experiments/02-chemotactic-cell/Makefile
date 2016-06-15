RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)

all: exp-intro help
include ../Makefile

simulation: chemo_cell.h5

%.h5: %.parameters $(RMPCDMD)/build/single_dimer_chemotactic_cell
	rmpcdmd run single_dimer_chemotactic_cell $< $@ auto

