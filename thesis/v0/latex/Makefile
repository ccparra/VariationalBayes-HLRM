.PHONY: all clean

MAIN := EncabezadoTesisMSc.tex
OUTPUT_DIR := ../document
PDF := $(OUTPUT_DIR)/Master_Thesis.pdf

all: $(PDF)

$(PDF): $(MAIN) Kap1/Kap1.tex Kap2/Kap2.tex Kap3/Kap3.tex Kap4/Kap4.tex Anexos/Anexos.tex BibliMSc.bib
	mkdir -p $(OUTPUT_DIR)
	tectonic -X compile $(MAIN) --outdir $(OUTPUT_DIR)
	mv $(OUTPUT_DIR)/EncabezadoTesisMSc.pdf $(PDF)

clean:
	$(RM) $(PDF)
