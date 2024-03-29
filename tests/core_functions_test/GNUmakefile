# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)
# GRCHOMBO_SOURCE := Set locally (e.g. export GRCHOMBO_SOURCE=path/to/GRChombo/Source in bash)

GRMILIJUN_SOURCE = $(shell pwd)/../../Source

ebase := Main_Tests

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRMILIJUN_SOURCE)/Background  \
            $(GRMILIJUN_SOURCE)/Matter \
			$(GRMILIJUN_SOURCE)/Diagnostics \
			$(GRMILIJUN_SOURCE)/Core \
            $(GRMILIJUN_SOURCE)/Lagrangians

include $(CHOMBO_HOME)/mk/Make.test