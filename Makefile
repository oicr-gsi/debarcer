#----------------------------------------
# File groups
#----------------------------------------

BINARIES = bin/levenshtien_stream
DOCS = docs/manual.pdf

#----------------------------------------
# Tools
#----------------------------------------

COMPILE = gcc
LATEX = pdflatex

#----------------------------------------
# Targets
#----------------------------------------

all : $(BINARIES) $(DOCS)

bin/levenshtien_stream : src/levenshtien_stream.c
	$(COMPILE) -o $@ $^

docs/manual.pdf : docs/manual.tex
	$(LATEX) -output-directory docs $^
