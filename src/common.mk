dir_guard=@mkdir -p $(@D)

$(OBJDIR)/%.cpp.o: %.cpp
	$(dir_guard)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -MMD -o $@ -c $<
 
$(OBJDIR)/%.cu.o: %.cu
	$(dir_guard)
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
	$(NVCC) $(CFLAGS) $(CXXFLAGS) $(INCLUDES) -M $< | sed 's#\([a-zA-Z0-9_-]*\)\.o#$(OBJDIR)/\1.cu.o#' > $(@:.o=.d)

objects = $(patsubst %.cu, $(OBJDIR)/%.cu.o, \
		  $(patsubst %.cpp, $(OBJDIR)/%.cpp.o, \
  		  $(1))) 

