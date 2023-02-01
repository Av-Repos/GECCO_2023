#ifndef _LIST_H
#define _LIST_H


struct listitem {
	int i;
	int j;
	struct listitem *next;
};

typedef struct listitem Listitem;

struct list {
	Listitem *head;
};

typedef struct list List;
void initlist (List *ilist);
void insertback(List *ilist, int i, int j);
void popfront(List *ilist);
int inlist(List ilist, int i, int j);
long long int length(List ilist);
void destroy(List * ilist);
#endif
