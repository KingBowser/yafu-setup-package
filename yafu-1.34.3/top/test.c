/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "arith.h"
#include "util.h"
#include "qs.h"
#include "factor.h"

double spAcc(int m);
double sqr_acc(int m, int sz);
double subadd_acc(int m, int sz);
double shortmuldiv_acc(int m, int sz);
double muldiv_acc(int m, int sz);
double sqrt_acc(int m, int sz);
double gcd_acc(int m, int sz);

void test_dlp_composites()
{
	FILE *in;
	uint64 *comp, f64;
	uint32 *f1;
	uint32 *f2, totBits,minBits,maxBits;
	double t_time;
	clock_t start, stop;
	int i,j,num,correct;
	z tmp,tmp2,t1,t2,t3;
	mpz_t gmptmp;
	//uint64 f64;
	//int64 queue[100];
	//fact_obj_t *fobj;

	mpz_init(gmptmp);
	comp = (uint64 *)malloc(2000000 * sizeof(uint64));
	f1 = (uint32 *)malloc(2000000 * sizeof(uint32));
	f2 = (uint32 *)malloc(2000000 * sizeof(uint32));
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&t2);
	zInit(&t1);
	zInit(&t3);

	in = fopen("pseudoprimes.dat","r");

	start = clock();
	i=0;
	totBits = 0;
	minBits = 999;
	maxBits = 0;
	//read in everything
	while (!feof(in))
	{
		fscanf(in,"%" PRIu64 ",%u,%u",comp+i,f1+i,f2+i);
		sp642z(comp[i],&tmp);
		j = zBits(&tmp);
		totBits += j;
		if ((uint32)j > maxBits)
			maxBits = j;
		if ((uint32)j < minBits && j != 0)
			minBits = j;
		i++;
	}
	num = i;
	num = 10000;
	fclose(in);
	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("data read in %2.4f sec\n",t_time);
	printf("average bits of input numbers = %.2f\n",(double)totBits/(double)i);
	printf("minimum bits of input numbers = %d\n",minBits);
	printf("maximum bits of input numbers = %d\n",maxBits);


	start = clock();

	correct = 0;
	for (i=0;i<num;i++)
	{
		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}

		//sp642z(comp[i],&tmp);
		mpz_set_64(gmptmp, comp[i]);
		f64 = sp_shanks_loop(gmptmp, NULL);

		if ( ((uint32)f64 == f1[i]) || ((uint32)f64 == f2[i]))
			correct++;
	}

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("squfof got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);

	/*
	start = clock();

	correct = 0;
	for (i=0;i<num;i++)
	{
		
		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}
		sp642z(comp[i],&tmp);
		//smallmpqs(&tmp,&tmp2,&t2,&t3);
		//f64 = tmp2.val[0];

		//f64 = sp_shanks_loop(&tmp);
		f64 = (uint64)squfof_jp(&tmp);

		if ( ((uint32)f64 == f1[i]) || ((uint32)f64 == f2[i]))
			correct++;
	}

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("squfof got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);

	start = clock();

	correct = 0;
	for (i=0;i<num;i++)
	{

		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}
		sp642z(comp[i],&tmp);
		f64 = (uint64)SQUFOF_alpertron((int64)z264(&tmp),queue);

		if ( ((uint32)f64 == f1[i]) || ((uint32)f64 == f2[i]))
			correct++;
	}

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("squfof got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);

	start = clock();

	correct = 0;
	for (i=0;i<num;i++)
	{
		int fact1, fact2;
		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}
		sp642z(comp[i],&tmp);
		////sp642z(2128691,&tmp);
		squfof_rds((int64)z264(&tmp),&fact1, &fact2);
		//printf("%" PRIu64 " = %u * %u\n",z264(&tmp), fact1, fact2);


		if ( ((uint32)fact1 == f1[i]) || ((uint32)fact1 == f2[i]))
			correct++;
	}

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("squfof got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);	

	*/

	/*
	fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj);


	start = clock();

	correct = 0;
	num=10000;
	for (i=0;i<num;i++)
	{
		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}
		sp642z(comp[i],&tmp);
		zCopy(&tmp,&fobj->qs_obj.n);
		fobj->qs_obj.flags = 12345;
		pQS(fobj);

		for (j=0; j<fobj->qs_obj.num_factors; j++)
		{
			uint32 fact = (uint32)fobj->qs_obj.factors[j].val[0];
			if ( (fact == f1[i]) || (fact == f2[i]))
			{
				correct++;
				break;
			}
			zFree(&fobj->qs_obj.factors[j]);
		}
		for (; j < fobj->qs_obj.num_factors; j++)
			zFree(&fobj->qs_obj.factors[j]);

	}
	fobj->qs_obj.num_factors = 0;
	free_factobj(fobj);

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("pQS got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);	
	

	fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj);

	start = clock();

	correct = 0;
	for (i=0;i<num;i++)
	{
		if (i%1000 == 0)
		{
			printf("done with %d, %d correct, after %2.2f sec\n",i,
				correct,(double)(clock() - start)/(double)CLOCKS_PER_SEC);
		}

		sp642z(comp[i],&tmp);
		zCopy(&tmp,&fobj->qs_obj.n);
		smallmpqs(fobj);

		for (j=0; j<fobj->qs_obj.num_factors; j++)
		{
			uint32 fact = (uint32)fobj->qs_obj.factors[j].val[0];
			if ( (fact == f1[i]) || (fact == f2[i]))
			{
				correct++;
				break;
			}
			zFree(&fobj->qs_obj.factors[j]);
		}
		for (; j < fobj->qs_obj.num_factors; j++)
			zFree(&fobj->qs_obj.factors[j]);
		fobj->qs_obj.num_factors = 0;

	}
	fobj->qs_obj.num_factors = 0;
	free_factobj(fobj);

	stop = clock();
	t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("smallmpqs got %d of %d correct in %2.2f sec\n",correct,num,t_time);
	printf("percent correct = %.2f\n",100.0*(double)correct/(double)num);
	printf("average time per input = %.2f ms\n",1000*t_time/(double)num);	
	*/

	free(f1);
	free(f2);
	free(comp);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&t2);
	zFree(&t1);
	zFree(&t3);
	mpz_clear(gmptmp);
	return;
}

void test_qsort(void)
{
	//test the speed of qsort in  sorting a few million lists
	//of a few thousand random integers

	uint32 *list;
	int i,j,k;
	const int listsz = 512;
	double t_time;
	clock_t start, stop;

	list = (uint32 *)malloc(listsz * sizeof(uint32));

	for (k=1000;k<1000000;k*=10)
	{
		start = clock();

		for (j=0; j<k; j++)
		{
			for (i=0; i<listsz; i++)
				list[i] = (uint32)spRand(1000000,50000000);
		}

		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("baseline elapsed time for %d lists = %2.2f\n",k,t_time);

		start = clock();

		for (j=0; j<k; j++)
		{
			for (i=0; i<listsz; i++)
				list[i] = (uint32)spRand(1000000,50000000);
			qsort(list,listsz,4,&qcomp_uint32);
		}

		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time for %d sorts = %2.2f\n",k,t_time);
	}

	free(list);
	return;
}
	
void arith_timing(int num)
{

	int i,sz;
	double t_time;
	clock_t start, stop;
	z a,b,c,d,e;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&e);
	
	//spAcc(num);
	
	//sqrt_acc(num, int sz);
	//gcd_acc(num, int sz);
	goto muldiv;

	for (sz=50; sz<=500; sz += 50)
	{
		printf("baseline: generating %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			zRand(&a,sz);
			zRand(&b,sz);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
	}

	for (sz=50; sz<=500; sz += 50)
	{
		subadd_acc(num, sz);
		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("adding %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zAdd(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

	for (sz=50; sz<=500; sz += 50)
	{
		shortmuldiv_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("subtracting %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zSub(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

muldiv:
	//for (sz=50; sz<=500; sz += 50)
	for (sz=5; sz<=50; sz += 1)
	{
		muldiv_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("Multiplying %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zMul(&a,&b,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

	goto sqrt_test;

//sqr_test:
	for (sz=50; sz<=500; sz += 50)
	{
		sqr_acc(num, sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz);
		printf("Squaring %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz);
			zSqr(&a,&c);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}

sqrt_test:	
	for (sz=5; sz<=50; sz += 1)
	//for (sz=50; sz<=150; sz += 10)
	{
		//gcd_acc(num/100,sz);
		sqrt_acc(num,sz);

		/*
		zRand(&a,sz);
		zRand(&b,sz/2);
		printf("Dividing %d random numbers with %d digits: ",num,sz);
		start = clock();
		for (i=0;i<num;i++)
		{
			//zRand(&a,sz);
			//zRand(&b,sz/2);
			zCopy(&a,&e);
			zDiv(&a,&b,&c,&d);
			zCopy(&e,&a);
		}
		stop = clock();
		t_time = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed time = %2.3f\n",t_time);
		*/
	}
	

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&e);

	return;
}

double spAcc(int m)
{
	z a,b,c,d;
	unsigned int shift;
	fp_digit bitmask,k,v;
	fp_digit tmpd[2];
	int j;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test of single precision stuff, %d iterations\n",m);
	start = clock();
	for (j=1; j<m; ++j)
	{
		a.size=2;
		b.size=1;
		a.val[0] = rand()*rand();
		a.val[1] = rand()*rand();
		b.val[0] = 0;
		while (!(b.val[0]))
			b.val[0] = rand()*rand();

		bitmask = HIBITMASK;
		for (shift = 0; shift < BITS_PER_DIGIT; shift++)
		{
			if (b.val[0] & bitmask)
				break;
			bitmask >>= 1;
		}

		b.val[0] <<= shift;
		zShiftLeft(&a,&a,shift);

		tmpd[1] = a.val[1]; tmpd[0] = a.val[0];
		c.val[2] = spDivide(&c.val[1],&c.val[0],tmpd,b.val[0]);

		spMultiply(c.val[1],b.val[0],&d.val[0],&d.val[1]);
		if (c.val[2])
			d.val[1] += b.val[0];
		spAdd(d.val[0],c.val[0],&k,&v);
		spAdd(d.val[1],v,&d.val[1],&v);
		if ((k != a.val[0]) || (d.val[1] != a.val[1]))
		{
			printf("error at %d\na: ",j);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);

	return t;
}

double sqr_acc(int m, int sz)
{
	int i,j;
	clock_t start, stop;
	double t=0.0;

	z a,b,c,d;
	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test, square and verify by mul:\n");
	start = clock();
	for (j=1; j<m; ++j)
	{
		zRand(&a,sz);
		
		if (rand() > 16384)
			a.size *= -1;

		zSqr(&a,&b);
		//zNroot(&b,&c,2);
		zMul(&a,&a,&c);

		//if (a.size < 0)
		//	c.size *= -1;

		i = zCompare(&c,&b);

		if (!(i==0))
		{
			printf("failed at %d\n",j);
			printf("a = %s\nb = %s\nc = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&c,&gstr3));
			exit(-1);
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);

	return t;
}

double subadd_acc(int m, int sz)
{

	int i,j;
	z a,b,c,d;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);

	printf("\nAccuracy test, subtract and verify by addition:\n");
	start = clock();
	for (j=1; j<m; ++j)
	{
		zRand(&a,sz);
		zRand(&b,sz);
		
		if (rand() > 16384)
			a.size *= -1;

		if (rand() > 16384)
			b.size *= -1;

		zAdd(&a,&b,&d);
		zSub(&d,&b,&c);

		i = zCompare(&c,&a);

		if (!(i==0))
		{
			printf("failed at %d\n",j);
			printf("a+b=d\nd-b=c\nassert a == c failed\na: ");
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	return t;
}

double shortmuldiv_acc(int m, int sz)
{

	int i,j;
	z a,b,c,d,q;
	fp_digit v;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&q);

	printf("\nAccuracy test, short divide and verify by short multiplication and addition:\n");
	start = clock();
	for (j=1;j<m;++j)
	{
		zRand(&a,sz);
		
		//make a non-zero digit 'v'
		v = spRand(1,MAX_DIGIT);

		b.val[0] = zShortDiv(&a,v,&q); 
		b.size=1;

		zShortMul(&q,v,&c);
		zAdd(&c,&b,&d);
		zShortAdd(&d,b.val[0],&c);
		zShortSub(&c,b.val[0],&d);
		i = zCompare(&a,&d);
		
		if (!(i==0)) 
		{
			printf("failed at %d\na: ",j);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&q);
	return t;
}

double muldiv_acc(int m, int sz)
{
	int i,j;
	z a,b,c,d,q,rem,tmp;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&q);
	zInit(&rem);
	zInit(&tmp);

	printf("\nAccuracy test, divide and verify by multiplication and addition:\n");
	start = clock();
	for (j=1;j<m;++j)
	{	
		zRand(&a,sz*2);

		do 
		{
			zRandb(&b,(int)((double)sz * 3.33));
		} while (zCompare(&b,&zZero) == 0);

		zCopy(&a,&d);	//because a will be overwritten
		zDiv(&a,&b,&q,&rem);
		zMul(&q,&b,&c);
		zAdd(&c,&rem,&a);

		i = zCompare(&a,&d);
		
		if (i != 0) 
		{
			printf("\nfailed at %d\n",j);
			printf("a = %s\nb = %s\nq = %s\n",
				z2decstr(&d,&gstr1),z2decstr(&b,&gstr2),z2decstr(&q,&gstr3));
			printf("r = %s\nqb = %s\nqb+r = %s\n",
				z2decstr(&rem,&gstr1),z2decstr(&c,&gstr2),z2decstr(&a,&gstr3));
			exit(-1);
			break;
		}
	}
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", m,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&q);
	zFree(&rem);
	zFree(&tmp);
	return t;
}

double sqrt_acc(int m, int sz)
{
	int i,j,k;
	z a,b,c,d,tmp;

	clock_t start, stop;
	double t=0.0;

	zInit(&a);
	zInit(&b);
	zInit(&c);
	zInit(&d);
	zInit(&tmp);

	printf("\nAccuracy test, Newton square root:\n");
	start = clock();
	for (j=1;j<m;++j)
	{
		zRand(&a,sz);
		k = zNroot(&a,&b,2);
		
		if (k >= 10000)
		{
			printf("error, max iterations at %d\na:       ",j);
			break;
		}
		
		zMul(&b,&b,&c);
		i = zCompare(&c,&a);
		if (!((i==0) || (i==MAX_DIGIT)))
		{
			printf("error, not <= 0 at %d, iterations: %d\na:        ",j,k);			
			printf("a = %s\nb = %s\nc = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&c,&gstr3));

			break;
		}
		zShortAdd(&b,1,&c);
		zMul(&c,&c,&d);
		i = zCompare(&d,&a);
		if (i<=0)
		{
			printf("error, not > 0 at %d, iterations: %d\na:       ",j,k);
			printf("a = %s\nb = %s\nd = %s\n",
				z2decstr(&a,&gstr1),z2decstr(&b,&gstr2),z2decstr(&d,&gstr3));
			break;
		}
	}
	stop = clock();

	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("%d numbers verified.  Elapsed time = %6.4f seconds.\n", j,t);

	zFree(&a);
	zFree(&b);
	zFree(&c);
	zFree(&d);
	zFree(&tmp);
	return t;
}

