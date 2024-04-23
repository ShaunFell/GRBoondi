define SEARCH_FOR_MAKE
$(wildcard $1/*/GNUmakefile) $(wildcard $1/*/makefile) $(wildcard $1/*/Makefile)
endef

TestDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, Tests)
TestDirs := $(dir $(TestDirsWithGNUmakefile))
RunTestDirs := $(TestDirs:%=run-%)
ExampleDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, Examples)
CleanTestDirs := $(TestDirs:%=clean-%)
RealCleanTestDirs := $(TestDirs:%=realclean-%)

.PHONY: all run examples clean realclean $(TestDirs) $(ExampleDirs)

export GRCHOMBO_SOURCE = $(shell pwd)/Source

ECHO?=@ # set this to null on the command line to increase verbosity

test: $(TestDirs)

run: test $(RunTestDirs)

clean: $(CleanTestDirs) 

realclean: $(RealCleanTestDirs) 

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(ECHO)$(MAKE) -C $@ --no-print-directory all

$(RunTestDirs):
	$(info ################# Running test $@ #################)
	$(ECHO)$(MAKE) -C $(@:run-%=%) --no-print-directory run

$(CleanTestDirs):
	$(ECHO)$(MAKE) -C $(@:clean-%=%) --no-print-directory clean NODEPENDS=TRUE

$(RealCleanTestDirs):
	$(ECHO)$(MAKE) -C $(@:realclean-%=%) --no-print-directory realclean NODEPENDS=TRUE

