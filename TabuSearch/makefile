all: TS

TS: main.o EvaluateQAP.o list.o tabu_search.o RandomPerm.o
	gcc -std=c99 -O2 main.o EvaluateQAP.o list.o tabu_search.o RandomPerm.o -o TS -lm

main.o: main.c
	gcc -std=c99 -c main.c

EvaluateQAP.o: EvaluateQAP.c
	gcc -std=c99 -c EvaluateQAP.c

list.o: list.c
	gcc -std=c99 -c list.c

tabu_search.o: tabu_search.c
	gcc -std=c99 -c tabu_search.c

RandomPerm.o: RandomPerm.c
	gcc -std=c99 -c RandomPerm.c

clean:
	rm -rf *o TS
