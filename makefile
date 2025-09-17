# Makefile for OuterSpace Paper
# Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.

# Variables
TEX_DIR = tex
BIB_DIR = bib
FIG_DIR = figures
BUILD_DIR = build
MAIN_TEX = $(TEX_DIR)/main.tex
MAIN_PDF = $(BUILD_DIR)/main.pdf
MAIN_DOCX = $(BUILD_DIR)/main.docx

# LaTeX compiler settings
LATEX = pdflatex
BIBTEX = bibtex
LATEX_FLAGS = -interaction=nonstopmode -halt-on-error -output-directory=$(BUILD_DIR)

# Pandoc settings for DOCX conversion
PANDOC = pandoc
PANDOC_FLAGS = --bibliography=$(BIB_DIR)/bibliography.bib --citeproc

# Conda/environment settings
CONDA = conda
CONDA_ENV_FILE = environment.yml

# Default target
.PHONY: all
all: pdf

# Environment setup
.PHONY: venv
venv:
	@echo "Setting up conda environment for paper building..."
	$(CONDA) env create -f $(CONDA_ENV_FILE) -p venv || \
	$(CONDA) env update -f $(CONDA_ENV_FILE) -p venv
	@echo "Environment created/updated successfully."
	@echo ""
	@echo "IMPORTANT: You need to install LaTeX separately for PDF building:"
	@echo "  Ubuntu/Debian: sudo apt-get install texlive-full"
	@echo "  macOS: brew install --cask mactex"
	@echo "  Windows: Install MiKTeX or TeX Live"
	@echo ""
	@echo "To activate: conda activate ./venv"

.PHONY: venv-remove
venv-remove:
	@echo "Removing conda environment..."
	rm -rf venv

.PHONY: install-latex-ubuntu
install-latex-ubuntu:
	@echo "Installing LaTeX on Ubuntu/Debian..."
	sudo apt-get update && sudo apt-get install -y texlive-full
	@echo "LaTeX installation complete."

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build PDF
.PHONY: pdf
pdf: $(MAIN_PDF)

$(MAIN_PDF): $(MAIN_TEX) $(BIB_DIR)/bibliography.bib | $(BUILD_DIR)
	$(LATEX) $(LATEX_FLAGS) $(MAIN_TEX)
	cp $(MAIN_TEX) $(BUILD_DIR)/
	cd $(BUILD_DIR) && $(BIBTEX) main
	$(LATEX) $(LATEX_FLAGS) $(MAIN_TEX)
	$(LATEX) $(LATEX_FLAGS) $(MAIN_TEX)

# Build DOCX using Pandoc
.PHONY: docx
docx: $(MAIN_DOCX)

$(MAIN_DOCX): $(MAIN_TEX) $(BIB_DIR)/bibliography.bib | $(BUILD_DIR)
	$(PANDOC) $(MAIN_TEX) -o $(MAIN_DOCX) $(PANDOC_FLAGS)

# Clean build artifacts
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

# Clean and rebuild
.PHONY: rebuild
rebuild: clean all

# Watch for changes and rebuild (requires inotify-tools)
.PHONY: watch
watch:
	while inotifywait -e modify $(TEX_DIR)/*.tex $(BIB_DIR)/*.bib; do make pdf; done

# Show help
.PHONY: help
help:
	@echo "OuterSpace Paper Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  venv               - Create/update conda environment"
	@echo "  venv-remove        - Remove conda environment"
	@echo "  install-latex-ubuntu - Install LaTeX on Ubuntu/Debian systems"
	@echo "  pdf                - Build PDF version of the paper (default)"
	@echo "  docx               - Build DOCX version using Pandoc"
	@echo "  all                - Build PDF (same as 'pdf')"
	@echo "  clean              - Remove build artifacts"
	@echo "  rebuild            - Clean and rebuild"
	@echo "  watch              - Watch for changes and auto-rebuild PDF"
	@echo "  check-deps         - Check if required dependencies are installed"
	@echo "  help               - Show this help message"
	@echo ""
	@echo "Setup Instructions:"
	@echo "  1. Install conda/mamba if not already installed"
	@echo "  2. Run 'make venv' to set up the environment"
	@echo "  3. Install LaTeX:"
	@echo "     - Ubuntu/Debian: make install-latex-ubuntu (or sudo apt install texlive-full)"
	@echo "     - macOS: brew install --cask mactex"
	@echo "     - Windows: Install MiKTeX or TeX Live manually"
	@echo "  4. Activate environment: conda activate ./venv"
	@echo "  5. Run 'make pdf' to build the paper"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - Conda or Mamba package manager"
	@echo "  - System-level LaTeX installation (see setup instructions above)"

# Check prerequisites
.PHONY: check-deps
check-deps:
	@echo "Checking dependencies..."
	@which $(LATEX) >/dev/null 2>&1 || (echo "ERROR: $(LATEX) not found. Please install a LaTeX distribution." && exit 1)
	@which $(BIBTEX) >/dev/null 2>&1 || (echo "ERROR: $(BIBTEX) not found. Please install a LaTeX distribution." && exit 1)
	@which $(PANDOC) >/dev/null 2>&1 || (echo "WARNING: $(PANDOC) not found. DOCX conversion will not work.")
	@echo "Dependencies check completed."