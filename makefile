all:
	gcc -fPIC -lgomp -shared -o chromagram.so chromagramret.c

