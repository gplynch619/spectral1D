CPP := mpicxx
CPPFLAGS := -O3 -std=c++11 -MMD -MP

#Here. We are the source directory. 
SRCDIR := $(shell pwd)
BIN_DIR = $(shell pwd)
BUILD_DIR := build
#List of every .cxx file in the directory
SOURCE_FILES := $(wildcard *.cxx) 

#Finds the folder with FFTW3 executables in it 
FFTW3_DIR = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))

#For every file in SOURCE_FILES, makes a list with a corresponding name
#but with .cxx -> .o
OBJECTS := $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

#Includes headers from here
HEADERS := -I$(SRCDIR) 

#Linking flags
LDFLAGS := -lfftw3

#Dependency Options
DEPENDENCY_OPTIONS = -MM

PROJECT := $(BIN_DIR)/sine_test $(BIN_DIR)/solve

###############################################
#	Makefile cheatsheet
#	$< : the name of the first pre-requisite
#	$@ : the name of the target
#	@ : (before a command) don't echo this
#
###############################################

# -MT makes the target the follow string
# -MF specifies the file to write the dependencies to
$(BUILD_DIR)/%.cxx.d : %.cxx | $(BUILD_DIR)
	@$(CPP) $(DEPENDENCY_OPTIONS) $(HEADERS) $< -MT "$(BUILD_DIR)/$*.o $(BUILD_DIR)/$*.d" -MF $(BUILD_DIR)/$*.cxx.d 
	@echo "[ $(CPP) ] $@"

all : $(DEPS) $(PROJECT)

$(BUILD_DIR)/%.cxx.o : %.cxx | $(BUILD_DIR)
	@$(CPP) $(CPP_FLAGS) $(HEADERS) -c $< -o $@

$(BIN_DIR)/sine_test: $(BUILD_DIR)/fft_test.cxx.o
	@$(CPP) $(CPPFLAGS) -o $@ $< $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/solve: $(BUILD_DIR)/solve.cxx.o $(BUILD_DIR)/initread.cxx.o
	@$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"


#create directories if don't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

ifneq "$(strip $(DEPS))" ""
-include $(DEPS)
endif

.PHONY: depclean
depclean:
	rm -f $(DEPS)

.PHONY: clean
clean:
	rm -rf $(PROJECT) $(OBJECTS)

clean-all: clean depclean
