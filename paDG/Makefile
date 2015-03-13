PROCESSOR := $(shell uname -m)
F90			 =gfortran

SRCDIR = _src
OBJDIR = _objects
APPSDIR = _apps
APPLICATION = _apps/acoustics/variableCoeff/

FFLAGS  ?=-g -C -O4 -ffree-form -fbounds-check
OBJFLAGS ?=$(FFLAGS) -J$(OBJDIR)
LDFLAGS ?=-I$(OBJDIR) -I$(APPSDIR)

.PHONY= 1d_test clean

SOURCES := $(SRCDIR)/driver.f90 \
					 $(SRCDIR)/projectL2.f90 \
					 $(SRCDIR)/outputData.f90 \
					 $(SRCDIR)/evaluateExpansion.f90 \
					 $(SRCDIR)/forcingFunction.f90 \
           $(SRCDIR)/RiemannSolveLLF.f90 \
					 $(SRCDIR)/singleTimeStep.f90 \
					 $(SRCDIR)/singleStage.f90 \
					 $(SRCDIR)/updateCoeffs.f90 \
					 $(SRCDIR)/updateBuffer.f90 \

APPS :=    $(APPLICATION)/fluxFunction.f90 \
				   $(APPLICATION)/setWaveSpeed.f90 \
					 $(APPLICATION)/qinit.f90 \

MODULES := $(SRCDIR)/gridData.f90 \
					 $(SRCDIR)/testParameters.f90 \
					 $(SRCDIR)/modalDGmod.f90

INPUTS  := $(APPSDIR)/inputs.nl
OBJECTS :=$(addprefix $(OBJDIR)/,$(notdir $(SOURCES:.f90=.o)))
MODOBJ  :=$(addprefix $(OBJDIR)/,$(notdir $(MODULES:.f90=.o)))
APPSOBJ	:=$(addprefix $(OBJDIR)/,$(notdir $(APPS:.f90=.o)))

# Set the order in which directories are searched when looking for targets
VPATH = $(APPLICATION):$(APPSDIR):$(SRCDIR)

all: $(SOURCES) $(INPUTS) test_padg_1d

1d_test: test_padg_1d $(INPUTS)
	cp $(INPUTS) _output
	./test_padg_1d

test_padg_1d: $(MODOBJ) $(OBJECTS) $(APPSOBJ) execute.f90
	$(F90) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(MODOBJ) $(APPSOBJ) test_padg_1d

$(OBJDIR)/%.o $(OBJDIR)/%.mod: %.f90
	$(F90) $(OBJFLAGS) -c -o $@ $<
