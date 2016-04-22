BINARIES = tools/levenshtien_stream
COMPILE = gcc

all : $(BINARIES)

tools/levenshtien_stream : tools/levenshtien_stream.c
	$(COMPILE) -o $@ $^
