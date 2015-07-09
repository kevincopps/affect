.PHONY: build test clean doc

BOLD = `tput bold`
NORMAL = `tput sgr0`
BUILD_PRINT = @echo "${BOLD}Building package and Cython extensions.${NORMAL}"
CLEAN_PRINT = @echo "${BOLD}Removing previous build results.${NORMAL}"

build:
	$(BUILD_PRINT)
	@cd affect && ./build.sh

test:
	cd tests && nosetests -v -s

clean:
	$(CLEAN_PRINT)
	@cd affect && rm -rf exodus.c exodus.so build

doc:
	@cd doc && make html
