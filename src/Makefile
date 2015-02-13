LIBS     =   -lm

CFLAGS   =  -O3 -g -pg

CC       =   c++

OBJECTS  =   main.o AncestralProbabilities.o BranchChanges.o CharacterChanges.o Clades.o Data.o MbRandom.o Mcmc.o Model.o Msg.o Node.o Settings.o Tree.o

PROGS    = dolly

all:		$(PROGS)

dolly:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o dolly

main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

AncestralProbabilities.o:	AncestralProbabilities.cpp
		$(CC) $(CFLAGS) -c AncestralProbabilities.cpp

BranchChanges.o:	BranchChanges.cpp
		$(CC) $(CFLAGS) -c BranchChanges.cpp

CharacterChanges.o:	CharacterChanges.cpp
		$(CC) $(CFLAGS) -c CharacterChanges.cpp

Clades.o:	Clades.cpp
		$(CC) $(CFLAGS) -c Clades.cpp

Data.o:	Data.cpp
		$(CC) $(CFLAGS) -c Data.cpp

MbRandom.o:	MbRandom.cpp
		$(CC) $(CFLAGS) -c MbRandom.cpp

Mcmc.o:	Mcmc.cpp
		$(CC) $(CFLAGS) -c Mcmc.cpp

Model.o:	Model.cpp
		$(CC) $(CFLAGS) -c Model.cpp

Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

Settings.o:	Settings.cpp
		$(CC) $(CFLAGS) -c Settings.cpp

Tree.o:	Tree.cpp
		$(CC) $(CFLAGS) -c Tree.cpp

clean:		
		rm -f $(OBJECTS) $(PROGS)
