.PHONY: build build_trace clean coverage test docs test_pip_install

BOLD = `tput bold`
NORMAL = `tput sgr0`
BUILD_PRINT = @echo "${BOLD}Building package and Cython extensions.${NORMAL}"
BUILD_TRACE_PRINT = @echo "${BOLD}Rebuilding with tracing enabled.${NORMAL}"
COVERAGE_PRINT = @echo "${BOLD}Analyzing coverage from tests.${NORMAL}"
CLEAN_PRINT = @echo "${BOLD}Removing previous build results.${NORMAL}"
TEST_PIP_INSTALL_PRINT = @echo "${BOLD}Testing pip install in empty environment.${NORMAL}"

PYTEST := $(shell command -v py.test 2> /dev/null)

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
	@python setup.py build_ext -i 2>&1 | python makewarn.py "\"Using deprecated NumPy API" 4 3

build_trace: clean
	$(BUILD_TRACE_PRINT)
	@python setup.py --linetrace build_ext -i 2>&1 | python makewarn.py "\"Using deprecated NumPy API" 4 3

coverage: build_trace
	$(COVERAGE)
	@py.test --cov-report term:skip-covered --cov-report annotate:cov_annotate --cov-report html:cov_html --cov=affect affect/tests/

test:
	pytest -v -s affect/tests

clean:
	$(CLEAN_PRINT)
	@rm -rf build
	@rm -rf affect.egg-info
	@find . -type d -name '__pycache__' -exec rm -rf {} +
	@cd affect && rm -rf util.cpp exodus.c* connect.cpp dynamics.c* *.so

docs:
	@cd docs && rm -rf _build && make html

test_pip_install:
	$(TEST_PIP_INSTALL_PRINT)
	$(call pip_install,affect_test_pip_install,$(CURDIR))

profile_test:
	python -m cProfile -o profile.pstat $(PYTEST) -s affect/tests/test_display.py::test_write_scene
