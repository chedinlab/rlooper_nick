# Compiler definitions
CC = g++
CR = gcc

# Directories
SRCDIR = src/
BINDIR = bin/

# Compiler flags
CFLAGS = -g -std=c++11
CRFLAGS =

# Object files list
OBJS = $(SRCDIR)gene.o $(SRCDIR)model.o $(SRCDIR)Rloop_model.o $(SRCDIR)simulation.o $(SRCDIR)structure.o $(SRCDIR)windower.o

# Targets
all: $(BINDIR)rlooper_sim

# Rule to create the bin directory if it doesn't exist
$(BINDIR):
	mkdir -p $(BINDIR)

# Build the rlooper_sim executable
$(BINDIR)rlooper_sim: $(BINDIR) $(OBJS) $(SRCDIR)rlooper_sim.o
	@echo "Building rlooper_sim executable."
	$(CC) -o $(BINDIR)rlooper_sim $(CFLAGS) $(OBJS) $(SRCDIR)rlooper_sim.o

# General rule for C++ object files
$(SRCDIR)%.o: $(SRCDIR)%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# General rule for C object files
$(SRCDIR)%.o: $(SRCDIR)%.c
	$(CR) $(CRFLAGS) -c $< -o $@

# Clean up generated files
clean:
	@echo "Cleaning object files."
	rm -f $(SRCDIR)*.o

# Clean everything including the bin folder
clobber: clean
	@echo "Removing bin directory and its contents."
	rm -rf $(BINDIR)
