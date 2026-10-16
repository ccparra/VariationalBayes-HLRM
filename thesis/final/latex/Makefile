.PHONY: all clean

MAIN := EncabezadoTesisMSc.tex
BUILD_DIR := ../document
PDF := $(BUILD_DIR)/Master_Thesis.pdf

all: $(PDF)

$(PDF): $(MAIN) Kap1/Kap1.tex Kap2/Kap2.tex Kap3/Kap3.tex Kap4/Kap4.tex Anexos/Anexos.tex BibliMSc.bib
	mkdir -p $(BUILD_DIR)
	tectonic -X compile $(MAIN) --outdir $(BUILD_DIR)
	mv $(BUILD_DIR)/EncabezadoTesisMSc.pdf $(PDF)

clean:
	$(RM) $(PDF)
