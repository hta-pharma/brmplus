check: build
	R-devel CMD check --as-cran brmplus_1.1.0.tar.gz

build:
	R-devel CMD build .
