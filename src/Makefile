# Makefile
# Generic Makefile for making cuda programs
#
BIN					:= sop-gpu
# flags
CUDA_INSTALL_PATH	:= /usr/local/cuda
OBJDIR				:= obj
INCLUDES			+= -I$(CUDA_INSTALL_PATH)/include
LIBS				+= -L$(CUDA_INSTALL_PATH)/lib64
CFLAGS				:= -O3 -g $(INCLUDES)
LDFLAGS				:= $(LIBS) -lrt -lm -lcudart
# compilers
#NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_20 --ptxas-options=-v
NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_13 --ptxas-options=-v -use_fast_math 
CXX					:= g++
LINKER				:= g++ -fPIC
# files
CXX_SOURCES			:= main.cpp \
	IO/dcdio.cpp \
	IO/pdbio.cpp \
	IO/topio.cpp \
	IO/configreader.cpp \
	Util/wrapper.cpp \
	Util/mystl.cpp \
	param_initializer.cpp \
	sop.cpp 
CU_SOURCES			:= gsop.cu
CXX_OBJS				:= $(patsubst %.cpp, $(OBJDIR)/%.cpp.o, $(CXX_SOURCES))
CU_OBJS				:= $(patsubst %.cu, $(OBJDIR)/%.cu.o, $(CU_SOURCES))
 
$(BIN): makedirs clean $(CU_OBJS) $(CXX_SOURCES) sop-top
	$(LINKER) -o $(BIN) $(CU_OBJS) $(CXX_SOURCES) $(LDFLAGS) $(CFLAGS) $(CXXFLAGS)
 
$(OBJDIR)/%.cpp.o: $(CXX_SOURCES)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ -c $<
 
$(OBJDIR)/%.cu.o: $(CU_SOURCES)
	$(NVCC) $(INCLUDES) -o $@ -c $<
	
sop-top:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o sop-top soptop.cpp IO/configreader.cpp IO/topio.cpp IO/pdbio.cpp Util/wrapper.cpp -lm

makedirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(OBJDIR)/IO
	mkdir -p $(OBJDIR)/Util

run: $(BIN)
	LD_LIBRARY_PATH=$(CUDA_INSTALL_PATH)/lib ./$(BIN)
 
clean:
	rm -f $(BIN) sop-top 
	find $(OBJDIR) -name '*.o' -delete
	
install:
	cp $(BIN) /usr/bin/$(BIN)
	cp sop-top	/usr/bin/sop-top

all: sop-top sop-gpu

.PHONY: all install run makedirs

