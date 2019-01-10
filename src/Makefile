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

PIZDAINT_MESA_DIR := /apps/daint/UES/jenkins/6.0.UP07/mc/easybuild/software/Mesa/17.2.8-CrayGNU-18.08/
PIZDAINT_MESA_LIB := ${PIZDAINT_MESA_DIR}/lib/libOSMesa.so
PIZDAINT_MESA_INC := ${PIZDAINT_MESA_DIR}/include

ifeq ($(HOST),$(filter $(HOST),sapling pizdaint))
LLVM_LIB := $(SOLEIL_LLVM_DIR)/install/lib/libLLVM-4.0.so
LLVM_SRC := $(SOLEIL_LLVM_DIR)/llvm-4.0.0.src
LLVM_SRC_XZ := $(SOLEIL_LLVM_DIR)/llvm-4.0.0.src.tar.xz
LLVM_DOWNLOAD_URL := http://releases.llvm.org/4.0.0/llvm-4.0.0.src.tar.xz
LLVM_BUILD := $(SOLEIL_LLVM_DIR)/build

${LLVM_SRC_XZ}:	${SOLEIL_LLVM_DIR}
	cd ${SOLEIL_LLVM_DIR} && wget ${LLVM_DOWNLOAD_URL}

${LLVM_SRC}:	${LLVM_SRC_XZ} ${SOLEIL_LLVM_DIR}
	cd ${SOLEIL_LLVM_DIR} && tar xJf ${LLVM_SRC_XZ}

${LLVM_BUILD}:	
	mkdir -p ${LLVM_BUILD}

${SOLEIL_LLVM_DIR}:
	mkdir -p ${SOLEIL_LLVM_DIR}

${LLVM_LIB}:	${LLVM_SRC} ${LLVM_BUILD}
	$(shell module load cmake)
ifeq ($(HOST),sapling)
else ifeq($(HOST),pizdaint)
	cd ${LLVM_BUILD} && cmake $(SOLEIL_LLVM_DIR)/llvm-4.0.0.src
endif
	cd ${LLVM_BUILD} && make
	cd ${LLVM_BUILD} && cmake -DCMAKE_INSTALL_PREFIX=$(SOLEIL_LLVM_DIR)/install -P cmake_install.cmake

endif

soleil.exec: soleil.o soleil_mapper.o config_schema.o json.o render.o renderImage.o renderParticles.o marchingCubes.o ${LLVM_LIB}
ifeq ($(HOST),pizdaint)
	$(CXX) -o $@ $^ $(LINK_FLAGS) -L$(LEGION_DIR)/language/llvm/install/lib -lGLU -lLLVM $(PIZDAINT_MESA_LIBS)
else ifeq ($(HOST),sapling)
	$(CXX) -o $@ $^ $(LINK_FLAGS) -L$(LEGION_DIR)/language/llvm/install/lib -lGLU -lOSMesa $(LLVM_LIB)
else
	$(CXX) -o $@ $^ $(LINK_FLAGS) -lGLU -lOSMesa -lLLVM
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
ifeq ($(HOST),pizdaint)
	$(CXX) $(CXXFLAGS) -c -o $@ -I${PIZDAINT_MESA_INC} $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

marchingCubes.o:	marchingCubes.cc
ifeq ($(HOST),pizdaint)
	$(CXX) $(CXXFLAGS) -DUSE_SOFTWARE_OPENGL -c -o $@ -I${PIZDAINT_MESA_INC} $<
else
	$(CXX) $(CXXFLAGS) -DUSE_SOFTWARE_OPENGL -c -o $@ $<
endif

renderImage.o:	renderImage.cc
ifeq ($(HOST),pizdaint)
	$(CXX) $(CXXFLAGS) -c -o $@ -I${PIZDAINT_MESA_INC} $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

renderParticles.o:	renderParticles.cc
ifeq ($(HOST),pizdaint)
	$(CXX) $(CXXFLAGS) -c -o $@ -I${PIZDAINT_MESA_INC} $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

render_standalone.o:	render_standalone.cc
ifeq ($(HOST),pizdaint)
	$(CXX) $(CXXFLAGS) -c -o $@ -I${PIZDAINT_MESA_INC} $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

render_standalone:  render_standalone.o  marchingCubes.o renderImage.o renderParticles.o
	$(CXX) -o $@ $^ -lGLU -lOSMesa

