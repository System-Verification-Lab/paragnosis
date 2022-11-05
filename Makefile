
DEFAULTBUILDDIR := build
DEBUGBUILDDIR	:= debug
RELEASEBUILDDIR	:= release
CURRENTDIR 		:= $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SOURCEDIR 		:= $(CURRENTDIR)

export GCC_COLORS='error=01;31:warning=01;38;5;172:caret=01;38;5;172:locus=01:quote=:note=01;32:range1=:range2=:fixit-insert=:fixit-delete=:diff-filename=:diff-hunk=:diff-delete=:diff-insert=:type-diff='
export CMAKE_C_COMPILER = $(shell which gcc)
export CMAKE_CXX_COMPILER = $(shell which g++)
export CMAKE_CXX_FLAGS += -fdiagnostics-color=auto
export CMAKE_C_FLAGS += -fdiagnostics-color=auto

MAKEFLAGS += --no-print-directory

.PHONY: build debug release clean help paragnosis all

all: release
	@$(MAKE) pg

pg:
	@cd src/paragnosis && pip3 install .
	cp src/paragnosis/bin/pg bin

build: $(DEFAULTBUILDDIR)/CMakeCache.txt
	@$(MAKE) -Wno-dev -C $(DEFAULTBUILDDIR) install

release: CMAKEFLAGS += -DCMAKE_BUILD_TYPE=Release
release: $(RELEASEBUILDDIR)/CMakeCache.txt
	@$(MAKE) -Wno-dev -C $(RELEASEBUILDDIR) install

debug: CMAKEFLAGS += -DCMAKE_BUILD_TYPE=Debug
debug: $(DEBUGBUILDDIR)/CMakeCache.txt
	@$(MAKE) -Wno-dev -C $(DEBUGBUILDDIR) install

%/CMakeCache.txt:
	@mkdir -p $*
	@cd $* && cmake -Wno-dev $(CMAKEFLAGS) $(SOURCEDIR)

clean:
	@rm -rf $(DEBUGBUILDDIR) $(DEFAULTBUILDDIR) $(RELEASEBUILDDIR)
	@cd src/paragnosis && rm -rf *.egg-info build dist

help:
	@echo "Usage:"
	@echo "    make [option]"
	@echo ""
	@echo "Options:"
	@echo "    build*    : create build directory (cmake) and compile"
	@echo "    release   : same as 'build', but with NDEBUG"
	@echo "    debug     : same as 'build', but with debug symbols"
	@echo "    clean     : remove build directories"
	@echo ""
	@echo "    * = default"
	@echo ""

