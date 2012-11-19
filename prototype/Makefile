

all: cell_model_opt.so

cell_model_opt.so: cell_model_opt.c
	 gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(HOME)/.local/lib/python2.6/site-packages/numpy/core/include -I/usr/include/python2.6 -o cell_model_opt.so cell_model_opt.c

cell_model_opt.c: cell_model_opt.pyx
	cython -a cell_model_opt.pyx
