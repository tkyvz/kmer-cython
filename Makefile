.PHONY: build clean

build:
	python setup.py build_ext --inplace

clean:
	rm -rf build/
	rm -f *.so
	rm -f resources/*.c
