all:
	rootcint -f Dict_maketree.cpp -c `root-config --cflags` -p LinkDef.h
	g++ -o maketree maketree.cc Aux.cc `root-config --cflags --glibs`  -rpath /Users/cperez/external/root_install/lib ###On MAC
	#g++ -o maketree maketree.cc Aux.cc `root-config --cflags --glibs`  ###On Linux

