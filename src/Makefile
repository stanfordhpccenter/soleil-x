SHELL=bash
# Required paths
ifndef LEGION_DIR
  $(error LEGION_DIR is not set)
endif
ifndef SOLEIL_DIR
  $(error SOLEIL_DIR is not set)
endif

# OS-specific options
ifeq ($(shell uname),Darwin)
  DYNLINK_PATH := DYLD_LIBRARY_PATH
else
  DYNLINK_PATH := LD_LIBRARY_PATH
endif

# CUDA options
USE_CUDA ?= 1

# HDF options
export USE_HDF ?= 1
export HDF_HEADER ?= hdf5.h
HDF_LIBNAME ?= hdf5

#C compiler options
CFLAGS += -ggdb -O0 -Wall -Werror -fno-strict-aliasing -I$(LEGION_DIR)/runtime -I$(LEGION_DIR)/runtime/legion -I$(IMAGE_COMPOSITOR_DIR)/include -I$(LEGION_DIR)/runtime/realm -march=native
CXXFLAGS += -std=c++11 -ggdb -O0 -Wall -Werror -fno-strict-aliasing -I$(LEGION_DIR)/runtime -I$(LEGION_DIR)/runtime/legion -I$(IMAGE_COMPOSITOR_DIR)/include -I$(LEGION_DIR)/runtime/realm -march=native 

CFLAGS += -DBOUNDS_CHECKS
CXXFLAGS += -DBOUNDS_CHECKS

# Regent options
export INCLUDE_PATH := .
ifdef HDF_ROOT
  export INCLUDE_PATH := $(INCLUDE_PATH);$(HDF_ROOT)/include
  export $(DYNLINK_PATH) := $($(DYNLINK_PATH)):$(HDF_ROOT)/lib
endif
REGENT := $(LEGION_DIR)/language/regent.py -g
ifeq ($(USE_CUDA), 1)
  REGENT_FLAGS := -fflow 0 -fopenmp 1 -finner 1 -fcuda 1 -fcuda-offline 1
else
  REGENT_FLAGS := -fflow 0 -fopenmp 1 -finner 1 -fcuda 0
endif

# Link flags
ifdef CRAYPE_VERSION
  LINK_FLAGS += -Bdynamic
  LINK_FLAGS += $(CRAY_UGNI_POST_LINK_OPTS) -lugni
  LINK_FLAGS += $(CRAY_UDREG_POST_LINK_OPTS) -ludreg
endif
LINK_FLAGS += -L$(LEGION_DIR)/bindings/regent -lregent
ifdef HDF_ROOT
  LINK_FLAGS += -L$(HDF_ROOT)/lib
endif
ifeq ($(USE_HDF), 1)
  LINK_FLAGS += -l$(HDF_LIBNAME)
endif
LINK_FLAGS += -lm

# Image compositor linkage
LINK_FLAGS += -L$(IMAGE_COMPOSITOR_DIR)/src -lImageCompositor
export INCLUDE_PATH := $(INCLUDE_PATH);$(IMAGE_COMPOSITOR_DIR)/include;$(LEGION_DIR)/runtime/realm

.PHONY: default all clean

default: soleil.exec

all: soleil.exec dom_host.exec render_standalone

clean:
	$(RM) *.exec *.o *-desugared.rg config_schema.h

%-desugared.rg: %.rg
	./desugar.py $< > $@

HOST = $(shell hostname)

HOST := $(patsubst daint%,pizdaint,$(HOST))
HOST := $(patsubst nid%,pizdaint,$(HOST))
HOST := $(patsubst sapling%,sapling,$(HOST))
HOST := $(patsubst sh-%,sherlock,$(HOST))








ifeq ($(HOST),sapling)
LLVM_LIB := $(SOLEIL_LLVM_DIR)/install/lib/libLLVM-4.0.so
LLVM_BUILD_LIB := $(SOLEIL_LLVM_DIR)/build/lib/libLLVM-4.0.so
LLVM_SRC := $(SOLEIL_LLVM_DIR)/llvm-4.0.0.src
LLVM_SRC_XZ := $(SOLEIL_LLVM_DIR)/llvm-4.0.0.src.tar.xz
LLVM_DOWNLOAD_URL := http://releases.llvm.org/4.0.0/llvm-4.0.0.src.tar.xz
LLVM_BUILD := $(SOLEIL_LLVM_DIR)/build
LLVM_INSTALL := $(SOLEIL_LLVM_DIR)/install
LLVM_INCLUDES=-I${LLVM_INSTALL}/include

${LLVM_SRC_XZ}:	${SOLEIL_LLVM_DIR}
	echo why are we downloading ${LLVM_SRC_XZ} again?
	ls -l ${LLVM_SRC_XZ}
	cd ${SOLEIL_LLVM_DIR} && wget ${LLVM_DOWNLOAD_URL}

${LLVM_SRC}:	${LLVM_SRC_XZ} ${SOLEIL_LLVM_DIR}
	cd ${SOLEIL_LLVM_DIR} && tar xJf ${LLVM_SRC_XZ}

${LLVM_BUILD}:	
	mkdir -p ${LLVM_BUILD}

${LLVM_INSTALL}:
	mkdir -p ${LLVM_INSTALL}

${SOLEIL_LLVM_DIR}:
	mkdir -p ${SOLEIL_LLVM_DIR}


${LLVM_BUILD_LIB}:	${LLVM_SRC} ${LLVM_BUILD} 
ifeq ($(HOST),sapling)
	cd ${LLVM_BUILD} && cmake -DLLVM_BUILD_LLVM_DYLIB=1 -DLLVM_TARGETS_TO_BUILD=X86 ${LLVM_SRC}
else ifeq ($(HOST),pizdaint)
	cd ${LLVM_BUILD} && cmake ${LLVM_SRC}
endif
	cd ${LLVM_BUILD} && make

${LLVM_LIB}: ${LLVM_BUILD} ${LLVM_BUILD_LIB}
	cd ${LLVM_BUILD} && cmake -DCMAKE_INSTALL_PREFIX=$(SOLEIL_LLVM_DIR)/install -P cmake_install.cmake

LLVM_PREREQ := ${LLVM_LIB}


else ifeq ($(HOST),pizdaint)

LLVM_STATIC_LINKAGE := -lLLVMLTO -lLLVMPasses -lLLVMObjCARCOpts -lLLVMSymbolize -lLLVMDebugInfoPDB -lLLVMDebugInfoDWARF -lLLVMCoverage -lLLVMTableGen -lLLVMLineEditor -lLLVMOrcJIT -lLLVMMIRParser -lLLVMObjectYAML -lLLVMLibDriver -lLLVMOption -lLLVMX86Disassembler -lLLVMX86AsmParser -lLLVMX86CodeGen -lLLVMGlobalISel -lLLVMSelectionDAG -lLLVMAsmPrinter -lLLVMDebugInfoCodeView -lLLVMDebugInfoMSF -lLLVMX86Desc -lLLVMMCDisassembler -lLLVMX86Info -lLLVMX86AsmPrinter -lLLVMX86Utils -lLLVMMCJIT -lLLVMInterpreter -lLLVMExecutionEngine -lLLVMRuntimeDyld -lLLVMCodeGen -lLLVMTarget -lLLVMCoroutines -lLLVMipo -lLLVMInstrumentation -lLLVMVectorize -lLLVMScalarOpts -lLLVMLinker -lLLVMIRReader -lLLVMAsmParser -lLLVMInstCombine -lLLVMTransformUtils -lLLVMBitWriter -lLLVMAnalysis -lLLVMObject -lLLVMMCParser -lLLVMMC -lLLVMBitReader -lLLVMProfileData -lLLVMCore -lLLVMSupport -lLLVMDemangle
PIZDAINT_MESA_DIR := /apps/daint/UES/jenkins/6.0.UP07/mc/easybuild/software/Mesa/17.2.8-CrayGNU-18.08/
MESA_LIB := ${PIZDAINT_MESA_DIR}/lib/libOSMesa.so
MESA_INC := -I${PIZDAINT_MESA_DIR}/include

PIZDAINT_EASYBUILD_LLVM_DIR := $(HOME)/easybuild/daint/broadwell/modules/all/LLVM/
${PIZDAINT_EASYBUILD_LLVM_DIR}:
	@echo ========================================================
	@echo ========================================================
	@echo Did not detect a local esybuild of LLVM 4.0.0
	@echo Please issue these commands and then try to make again:
	@echo module load daint-mc
	@echo module load EasyBuild-custom/cscs
	@echo eb ../scripts/LLVM-4.0.0-CrayGNU-18.08.eb 
	@echo module load LLVM-4.0.0-CrayGNU-18.08
	@echo ========================================================
	@echo ========================================================
	@exit -1
LLVM_PREREQ := ${PIZDAINT_EASYBUILD_LLVM_DIR}

endif



soleil.exec: soleil.o soleil_mapper.o config_schema.o json.o render.o renderImage.o renderParticles.o marchingCubes.o 
ifeq ($(HOST),pizdaint)
	@make ${LLVM_PREREQ}
	$(CXX) -o $@ $^ $(LINK_FLAGS) -lGLU ${MESA_LIB} ${LLVM_STATIC_LINKAGE}
else ifeq ($(HOST),sapling)
	@make ${LLVM_PREREQ}
	$(CXX) -o $@ $^ $(LINK_FLAGS) -lGLU -lOSMesa ${LLVM_LIB}
endif

soleil.o: soleil-desugared.rg soleil_mapper.h config_schema.h hdf_helper.rg dom-desugared.rg util.rg
	$(REGENT) soleil-desugared.rg $(REGENT_FLAGS)

dom_host.exec: dom_host.o config_schema.o json.o
	$(CXX) -o $@ $^ $(LINK_FLAGS)

dom_host.o: dom_host.rg config_schema.h dom-desugared.rg util.rg
	$(REGENT) dom_host.rg $(REGENT_FLAGS)

soleil_mapper.o: soleil_mapper.cc soleil_mapper.h config_schema.h
	$(CXX) $(CXXFLAGS) -c -o  $@ $<

config_schema.o config_schema.h: process_schema.rg config_schema.lua json.h util.rg
	$(REGENT) process_schema.rg config_schema.lua $(REGENT_FLAGS)

json.o: json.c json.h
	$(CC) $(CFLAGS) -c -o $@ $<

render.o: render.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ ${MESA_INC} $<

marchingCubes.o:	marchingCubes.cc 
	$(CXX) $(CXXFLAGS) -DUSE_SOFTWARE_OPENGL -c -o $@ ${MESA_INC} $<

renderImage.o:	renderImage.cc 
	$(CXX) $(CXXFLAGS) -c -o $@ ${MESA_INC} $<

renderParticles.o:	renderParticles.cc 
	$(CXX) $(CXXFLAGS) -c -o $@ ${MESA_INC} $<

render_standalone.o:	render_standalone.cc 
	$(CXX) $(CXXFLAGS) -c -o $@ ${MESA_INC}  $<

render_standalone:  render_standalone.o  marchingCubes.o renderImage.o renderParticles.o
	$(CXX) -o $@ $^ -lGLU ${MESA_INC}

