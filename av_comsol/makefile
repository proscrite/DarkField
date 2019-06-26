OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library
LIBGFORTRAN = /Users/pabloherrero/garfield_git/Library/


# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
CFLAGS += -g


LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += -L$(LIBGFORTRAN) -lgfortran
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm


av: av_comsol.C
	$(CXX) $(CFLAGS) av_comsol.C
	@echo $(info $$(CXX) is [${(CXX)}])
	$(CXX) -o av_comsol av_comsol.o $(LDFLAGS)
	rm av_comsol.o
