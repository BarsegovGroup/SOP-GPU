# Makefile
# Generic Makefile for making cuda programs
#
BIN					:= sop-gpu
# flags
CUDA_INSTALL_PATH	:= /usr/local/cuda
OBJDIR				:= obj
INCLUDES			+= -I$(CUDA_INSTALL_PATH)/include
LIBS				:= -L$(CUDA_INSTALL_PATH)/lib64
CFLAGS				:= -O3
LDFLAGS				:= -lrt -lm -lcudart
# compilers
#NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive -arch sm_20 --ptxas-options=-v
NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive -arch sm_11 --ptxas-options=-v -use_fast_math 
CC					:= g++
LINKER				:= g++ -fPIC
# files
C_SOURCES			:= main.c \
	dcdio.c \
	pdbio.c \
	topio.c \
	config_reader.c \
	param_initializer.c \
	sop.c 
CU_SOURCES			:= gsop.cu
C_OBJS				:= $(patsubst %.c, $(OBJDIR)/%.c.o, $(C_SOURCES))
CU_OBJS				:= $(patsubst %.cu, $(OBJDIR)/%.cu.o, $(CU_SOURCES))
 
$(BIN): makedirs clean $(CU_OBJS) gsop-top
	$(LINKER) -o $(BIN) $(CU_OBJS) $(C_SOURCES) $(LDFLAGS) $(INCLUDES) $(LIBS)
 
$(OBJDIR)/%.c.o: $(C_SOURCES)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
 
$(OBJDIR)/%.cu.o: $(CU_SOURCES)
	$(NVCC) $(INCLUDES) -o $@ -c $<
	
gsop-top:
	$(CC) -o sop-top soptop.c config_reader.c topio.c
 
makedirs:
	mkdir -p $(OBJDIR)

run: $(BIN)
	LD_LIBRARY_PATH=$(CUDA_INSTALL_PATH)/lib ./$(BIN)
 
clean:
	rm -f $(BIN) sop-top $(OBJDIR)/*.o
	
install:
	cp $(BIN) /usr/bin/$(BIN)
	cp sop-top	/usr/bin/sop-top