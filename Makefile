all: mpcc.c
	gcc -g -o mpcc mpcc.c

clean:
	rm -rf mpcc mpcc.dSYM