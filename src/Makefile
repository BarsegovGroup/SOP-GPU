CUDA_INSTALL_PATH	:= /usr/local/cuda
CUDA_ARCH 			:= sm_13

INCLUDES			+= -I$(CUDA_INSTALL_PATH)/include
LIBS				+= -L$(CUDA_INSTALL_PATH)/lib64
CFLAGS				:= -O3 -g $(INCLUDES)
LDFLAGS				:= $(LIBS) -lrt -lm -lcudart

# compilers
NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc -arch $(CUDA_ARCH) --ptxas-options=-v -use_fast_math 
CXX					:= g++
LINKER				:= g++ -fPIC

OBJDIR				:= obj
include common.mk

# files
SOPGPU_SOURCES	:= \
	main.cpp \
	param_initializer.cpp \
	Integrators/ran2.cpp \
	IO/dcdio.cpp \
	IO/pdbio.cpp \
	IO/topio.cpp \
	IO/configreader.cpp \
	Util/wrapper.cpp \
	Util/mystl.cpp \
	Updaters/dcd_manager.cpp \
	gsop.cu

SOPTOP_SOURCES 	:= \
	soptop.cpp \
	IO/configreader.cpp \
	IO/topio.cpp \
	IO/pdbio.cpp \
	Util/wrapper.cpp

all: sop-top sop-gpu

SOPGPU_OBJS			:= $(call objects, $(SOPGPU_SOURCES))
SOPTOP_OBJS			:= $(call objects, $(SOPTOP_SOURCES))
-include $(SOPGPU_OBJS:.o=.d)
-include $(SOPTOP_OBJS:.o=.d)
 
sop-gpu: $(SOPGPU_OBJS)
	$(LINKER) -o $@ $(SOPGPU_OBJS) $(LDFLAGS)

sop-top: $(SOPTOP_OBJS)
	$(LINKER) -o $@ $(SOPTOP_OBJS) $(LDFLAGS)

clean:
	rm -f sop-gpu sop-top 
	rm -rf "$(OBJDIR)"
	
.PHONY: makedirs clean all

