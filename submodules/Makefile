# order is important, libraries first
modules = sonLib jobTree bwa last cPecan graphmap
.PHONY: all %.all clean %.clean

all : ${modules:%=all.%} 

all.%:
	cd $* && make

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean
