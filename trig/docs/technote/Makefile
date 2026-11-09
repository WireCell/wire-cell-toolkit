# Build the SNB trigger tech note.
#
#   make RUNS=/path/to/runs.json OUT=/path/to/build
#
# RUNS describes the run outputs (see tools/make_results.py); OUT receives
# gen/ (generated figures, tables and numbers) and the PDF.

RUNS ?= runs.json
OUT ?= build
PYTHON ?= python3
NAME = snb-trigger-technote

SOURCES = $(NAME).tex preamble.tex $(wildcard sections/*.tex) $(wildcard diagrams/*.tex)

all: $(OUT)/$(NAME).pdf

gen: $(OUT)/gen/numbers.tex

$(OUT)/gen/numbers.tex: $(RUNS) tools/make_results.py ../../pytools/snbta.py ../../pytools/snbpoint.py
	$(PYTHON) tools/make_results.py --runs $(RUNS) --out $(OUT)/gen

$(OUT)/$(NAME).pdf: $(SOURCES) $(OUT)/gen/numbers.tex
	mkdir -p $(OUT)/sections $(OUT)/diagrams
	TEXINPUTS=.:$(OUT):$${TEXINPUTS} latexmk -pdf -interaction=nonstopmode -halt-on-error \
	  -outdir=$(OUT) $(NAME).tex

clean:
	latexmk -C -outdir=$(OUT) $(NAME).tex

.PHONY: all gen clean
