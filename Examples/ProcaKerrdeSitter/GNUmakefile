# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)
# GRCHOMBO_SOURCE := Also set locally (e.g. export GRCHOMBO_SOURCE=path/to/GRChombo/Source in bash)     

GRBoondi_SOURCE = ../../Source

ebase := Main_ProcaKdS

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRDZHADZHA_SOURCE)/Background \
            $(GRDZHADZHA_SOURCE)/Matter \
            $(GRBoondi_SOURCE)/Background  \
            $(GRBoondi_SOURCE)/Matter \
			$(GRBoondi_SOURCE)/Diagnostics \
			$(GRBoondi_SOURCE)/Core \
            $(GRBoondi_SOURCE)/Lagrangians \
            $(GRBoondi_SOURCE)/utils

include $(CHOMBO_HOME)/mk/Make.test