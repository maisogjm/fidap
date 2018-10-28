Dmultivar_multiple_regressionC: Dmultivar_multiple_regressionC.o imgio.o nrutil.o \
	Dlubksb.o Dludcmp.o
	$(CC) $(EXTRAS) Dmultivar_multiple_regressionC.o imgio.o nrutil.o Dlubksb.o \
	Dludcmp.o -o Dmultivar_multiple_regressionC -lm
Dmultivar_multiple_regressionC.o: Dmultivar_multiple_regressionC.c
	$(CC) $(EXTRAS) -c Dmultivar_multiple_regressionC.c
imgio.o: imgio.c
	$(CC) $(EXTRAS) -c imgio.c
nrutil.o: nrutil.c
	$(CC) $(EXTRAS) -c nrutil.c
Dlubksb.o: Dlubksb.c
	$(CC) $(EXTRAS) -c Dlubksb.c
Dludcmp.o: Dludcmp.c
	$(CC) $(EXTRAS) -c Dludcmp.c
