.PHONY: build test clean docs

BOLD = `tput bold`
NORMAL = `tput sgr0`
BUILD_PRINT = @echo "${BOLD}Building package and Cython extensions.${NORMAL}"
CLEAN_PRINT = @echo "${BOLD}Removing previous build results.${NORMAL}"
TEST_INSTALL_PRINT = @echo "${BOLD}Testing install in new environment.${NORMAL}"

TEST_INSTALL_ENV := affect_test_install

build:
	$(BUILD_PRINT)
	@python setup.py build_ext -i 2>&1 | python makewarn.py "warning: \"Using deprecated NumPy API" 4 3

test:
	pytest -v -s affect/tests/test_database.py

test_install:
	$(TEST_INSTALL_PRINT)
	conda create --yes --name ${TEST_INSTALL_ENV} --file requirements.txt
	source `which activate` ${TEST_INSTALL_ENV}
	pip install git+file://$(CURDIR)
	pip uninstall affect
	source `which deactivate` ${TEST_INSTALL_ENV}
	conda remove --yes --name ${TEST_INSTALL_ENV} --all


clean:
	$(CLEAN_PRINT)
	@rm -rf build
	@rm -rf affect.egg-info
	@cd affect && rm -rf exodus.cpp connect.cpp *.so __pycache__

docs:
	@cd docs && make html
