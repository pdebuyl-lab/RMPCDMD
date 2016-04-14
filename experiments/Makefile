RMPCDMD:=$(PWD)/../..

.PHONY: help code

exp-help: exp-intro
	@echo "Visit the subdirectories to access the simulation examples"

exp-intro:
	@echo "-- RMPCDMD software --"
	@echo "General information on the experiments available in RMPCDMD"
	@echo "is available in the file experiments/README.md"
	@echo

help:
	@echo "Available options for the make command"
	@echo "  simulation        : perform the example simulation"
	@echo "  code              : rebuild the source code"

code:
	make --no-print-directory -C $(RMPCDMD)/build

