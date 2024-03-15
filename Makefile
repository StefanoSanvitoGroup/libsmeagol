SHELL = /bin/sh

# the home dir is taken from the current directory
ifneq ($(BUILD_STAGE), true)
SMEAGOLHOME    := $(CURDIR)
export SMEAGOLHOME
endif

MAKEFILE       := $(SMEAGOLHOME)/Makefile
LIBDIR         := $(SMEAGOLHOME)/lib
OBJDIR         := $(SMEAGOLHOME)/obj
SRCDIR         := $(SMEAGOLHOME)/src

# Default Target ============================================================
LIBNAME      := smeagol
LIBRARY      := lib$(LIBNAME)
ARCHIVE_EXT  := .a
default_target: dirs $(LIBRARY)

# Declare PHONY targets =====================================================
.PHONY : dirs default_target $(LIBRARY) clean

ifeq ($(BUILD_STAGE),)
# create directories for library and object files and build SMEAGOL in these directories.
# The easiest way to change the build directory is to call make recursively.
$(LIBRARY): dirs
	@echo "SMEAGOL-FLAGS = " $(FCFLAGS)
	@+$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) $(LIBDIR)/$(LIBRARY)$(ARCHIVE_EXT) BUILD_STAGE=true

dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(LIBDIR)
endif

clean:
	rm -rf $(OBJDIR) $(LIBDIR)

ifeq ($(BUILD_STAGE), true)
# Set default values and read configuration and dependency files
MODDEPS = "lower"
include $(SMEAGOLHOME)/arch.make
include $(SMEAGOLHOME)/all.dep

# Discover files and directories
ALL_SRC_DIRS := $(shell find $(SRCDIR) -type d | awk '{printf("%s:",$$1)}')

vpath %.F     $(ALL_SRC_DIRS)
vpath %.f     $(ALL_SRC_DIRS)
vpath %.F90   $(ALL_SRC_DIRS)
vpath %.f90   $(ALL_SRC_DIRS)

# Build targets

%.o: %.F
	$(FC) -c $(FCFLAGS) $(FCFLAGS_FIXEDFORM) -I'$(dir $<)' -I'$(SRCDIR)' $<

%.o: %.f
	$(FC) -c $(FCFLAGS) $(FCFLAGS_FIXEDFORM) -I'$(dir $<)' -I'$(SRCDIR)' $<

%.o: %.F90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_FREEFORM) -I'$(dir $<)' -I'$(SRCDIR)' $<

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_FREEFORM) -I'$(dir $<)' -I'$(SRCDIR)' $<

%.mod:
	@true

$(LIBDIR)/%:
	@echo "Updating archive $@"
	@$(AR) $@ $?
endif
