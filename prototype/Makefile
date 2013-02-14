
PYTHONDIR=`python-config --includes`
NPDIR=$(HOME)/.local/lib/python2.6/site-packages/numpy/core/include

all: cell_model_opt.so

cell_model_opt.so: cell_model_opt.c
	 gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(NPDIR) -I$(PYTHONDIR) -o cell_model_opt.so cell_model_opt.c

cell_model_opt.c: cell_model_opt.pyx
	cython -a cell_model_opt.pyx
