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
all: $(BINDIR)rlooper

# Rule to create the bin directory if it doesn't exist
$(BINDIR):
	mkdir -p $(BINDIR)

# Build the rlooper executable
$(BINDIR)rlooper: $(BINDIR) $(OBJS) $(SRCDIR)rlooper.o
	@echo "Building rlooper executable."
	$(CC) -o $(BINDIR)rlooper $(CFLAGS) $(OBJS) $(SRCDIR)rlooper.o

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
