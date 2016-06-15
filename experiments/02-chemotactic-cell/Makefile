RMPCDMD:=$(PWD)/../..

EXPDIR:=$(PWD)

all: exp-intro help
include ../Makefile

simulation: chemo_cell.h5

%.h5: %.parameters $(RMPCDMD)/build/chemotactic_cell
	$(RMPCDMD)/build/rmpcdmd run chemotactic_cell $< $@ auto

