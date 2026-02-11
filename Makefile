# Makefile for GenotypeSimulator.jl

JULIA := julia +1.12
PROJECT_DIR := .
JULIAC_DIR := juliac
MAIN_SCRIPT := scripts/main.jl
BINARY_NAME := genotype_simulator
JULIAC_SCRIPT := $(shell $(JULIA) -e "print(joinpath(Sys.BINDIR, \"..\", \"share\", \"julia\", \"juliac\", \"juliac.jl\"))")

# Detect OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    BINARY := $(BINARY_NAME)
endif
ifeq ($(UNAME_S),Darwin)
    BINARY := $(BINARY_NAME)
endif
ifeq ($(OS),Windows_NT)
    BINARY := $(BINARY_NAME).exe
endif

.PHONY: all clean compile test install deps help

all: compile

help:
	@echo "GenotypeSimulator.jl Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  compile    - Compile the binary using juliac"
	@echo "  test       - Run the test suite"
	@echo "  deps       - Install dependencies"
	@echo "  clean      - Clean build artifacts"
	@echo "  install    - Install to ~/.local/bin (Unix only)"
	@echo "  help       - Show this help message"

deps:
	$(JULIA) --project=$(PROJECT_DIR) -e "using Pkg; Pkg.instantiate()"

$(JULIAC_DIR):
	mkdir -p $(JULIAC_DIR)

compile: deps | $(JULIAC_DIR)
	$(JULIA) --project=$(PROJECT_DIR) --experimental $(JULIAC_SCRIPT) \
		--output-exe $(JULIAC_DIR)/$(BINARY) \
		--verbose \
		$(MAIN_SCRIPT)
	@echo "Compilation successful!"
	@echo "Binary location: $(JULIAC_DIR)/$(BINARY)"
	@echo "Run with: ./$(JULIAC_DIR)/$(BINARY) --help"

test: deps
	$(JULIA) --project=$(PROJECT_DIR) -e "using Pkg; Pkg.test()"

install: compile
ifeq ($(OS),Windows_NT)
	@echo "Install target not supported on Windows. Copy $(JULIAC_DIR)/$(BINARY) to your PATH manually."
else
	mkdir -p ~/.local/bin
	cp $(JULIAC_DIR)/$(BINARY) ~/.local/bin/$(BINARY_NAME)
	chmod +x ~/.local/bin/$(BINARY_NAME)
	@echo "Installed to ~/.local/bin/$(BINARY_NAME)"
	@echo "Make sure ~/.local/bin is in your PATH"
endif

clean:
	rm -rf $(JULIAC_DIR)
	rm -f *.csv *.vcf *.ped *.map

.SECONDARY: