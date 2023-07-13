all: VFS

VFS: main.o EvaluateQAP.o list.o VFS.o RandomPerm.o
	gcc -std=c99 -O2 main.o EvaluateQAP.o list.o VFS.o RandomPerm.o -o VFS -lm

main.o: main.c
	gcc -std=c99 -c main.c

EvaluateQAP.o: EvaluateQAP.c
	gcc -std=c99 -c EvaluateQAP.c

list.o: list.c
	gcc -std=c99 -c list.c

VFS.o: VFS.c
	gcc -std=c99 -c VFS.c

RandomPerm.o: RandomPerm.c
	gcc -std=c99 -c RandomPerm.c

clean:
	rm -rf *o VFS
