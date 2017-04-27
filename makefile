.PHONY: build test clean docs

BOLD = `tput bold`
NORMAL = `tput sgr0`
BUILD_PRINT = @echo "${BOLD}Building package and Cython extensions.${NORMAL}"
CLEAN_PRINT = @echo "${BOLD}Removing previous build results.${NORMAL}"
TEST_PIP_INSTALL_PRINT = @echo "${BOLD}Testing pip install in empty environment.${NORMAL}"

#
# A function to create a new environment, activate it, and install our package, uninstall, and remove environment.
# $(1) is name of new environment
# $(2) is the current directory (root of our project)
# $(3) is a temporary working directory
#
pip_install =                                  \
	conda create --yes --name $(1) pytest;     \
	source `which activate` $(1);              \
	pip install git+file://$(2);               \
	pip uninstall --yes affect;                \
	source `which deactivate` $(1);            \
	conda remove --yes --name $(1) --all


build:
	$(BUILD_PRINT)
	@python setup.py build_ext -i 2>&1 | python makewarn.py "warning: \"Using deprecated NumPy API" 4 3

test:
	pytest -v -s affect/tests/test_database.py

test_pip_install:
	$(TEST_PIP_INSTALL_PRINT)
	$(call pip_install,affect_test_pip_install,$(CURDIR))

clean:
	$(CLEAN_PRINT)
	@rm -rf build
	@rm -rf affect.egg-info
	@cd affect && rm -rf exodus.cpp connect.cpp *.so __pycache__

docs:
	@cd docs && rm -rf _build && make html
