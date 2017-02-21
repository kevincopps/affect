.PHONY: build test clean docs

BOLD = `tput bold`
NORMAL = `tput sgr0`
BUILD_PRINT = @echo "${BOLD}Building package and Cython extensions.${NORMAL}"
CLEAN_PRINT = @echo "${BOLD}Removing previous build results.${NORMAL}"

build:
	$(BUILD_PRINT)
	@python setup.py build_ext -i 2>&1 | python makewarn.py "warning: \"Using deprecated NumPy API" 4 3

test:
	pytest -s affect/tests

clean:
	$(CLEAN_PRINT)
	@rm -rf build
	@cd affect && rm -rf exodus.cpp connect.cpp *.so __pycache__

docs:
	@cd docs && make html
