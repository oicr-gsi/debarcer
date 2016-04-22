#----------------------------------------
# Settings
#----------------------------------------

BINARIES = tools/levenshtien_stream
COMPILE = gcc

#----------------------------------------
# Targets
#----------------------------------------

all : $(BINARIES)

tools/levenshtien_stream : tools/levenshtien_stream.c
	$(COMPILE) -o $@ $^
