F: F.o golden.o betai.o nrutil.o gammln.o betacf.o
	$(CC) $(EXTRAS) F.o golden.o betai.o nrutil.o gammln.o betacf.o -o F -lm
F.o: F.c
	$(CC) $(EXTRAS) -c F.c
golden.o: golden.c
	$(CC) $(EXTRAS) -c golden.c
betai.o: betai.c
	$(CC) $(EXTRAS) -c betai.c
nrutil.o: nrutil.c
	$(CC) $(EXTRAS) -c nrutil.c
gammln.o: gammln.c
	$(CC) $(EXTRAS) -c gammln.c
betacf.o: betacf.c
	$(CC) $(EXTRAS) -c betacf.c
