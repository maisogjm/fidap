#include <stdio.h>

main()
{
    float p;
    char c;

    while(scanf("%f",&p)!=EOF){
	printf("%g",p);
	if((c=getchar())!=EOF)
	      printf("%c",c);
    }
}
