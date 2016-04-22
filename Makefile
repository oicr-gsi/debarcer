#----------------------------------------
# File groups
#----------------------------------------

BINARIES = tools/levenshtien_stream
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

tools/levenshtien_stream : tools/levenshtien_stream.c
	$(COMPILE) -o $@ $^

docs/manual.pdf : docs/manual.tex
	$(LATEX) -output-directory docs $^
