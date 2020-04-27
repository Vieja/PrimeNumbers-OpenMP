#include<math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <omp.h>
using namespace std;

#define M 2
#define N 100000000
#define THREADS 8
#define WYPISZ_PIERWSZE 1

void wypisz(string tekst, bool* czyPierwsza, int m, int n) {
	int ile = 0;
	for (int i = m; i < n + 1; i++) {
		if (czyPierwsza[i] == true) {
			ile++;
			if (ile <= 100) {
				printf("%d ", i);
				if (ile % 10 == 0) printf("\n");
			}
		}
	}
	cout << tekst << " znalazl " << ile << " liczb pierwszych\n\n";
}

void naiveSequentialOdd(int m, int n) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	czyPierwsza[2] = true;
	int mm = m % 2 == 0 ? m + 1 : m;
	for (int sprawdzany = mm; sprawdzany <= n; sprawdzany += 2)
	{
		bool liczbaPierwsza = true;
		for (int i = 3; i <= sqrt(sprawdzany); i += 2)
		{
			if (sprawdzany % i == 0)
			{
				liczbaPierwsza = false;
				break;
			}
		}
		czyPierwsza[sprawdzany] = liczbaPierwsza;
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("naiveSequentialOdd", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void naiveSequential(int m, int n) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	for (int sprawdzany = m; sprawdzany <= n; sprawdzany++)
	{
		bool liczbaPierwsza = true;
		for (int i = 2; i <= sqrt(sprawdzany); i++)
		{
			if (sprawdzany % i == 0)
			{
				liczbaPierwsza = false;
				break;
			}
		}
		czyPierwsza[sprawdzany] = liczbaPierwsza;
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("naiveSequential", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void naiveParallel(int m, int n, int threads) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	czyPierwsza[2] = true;
	int mm = m % 2 == 0 ? m + 1 : m;
	omp_set_num_threads(threads);
#pragma omp parallel for
	for (int sprawdzany = mm; sprawdzany <= n; sprawdzany += 2)
	{
		bool liczbaPierwsza = true;
		for (int i = 3; i <= sqrt(sprawdzany); i += 2)
		{
			if (sprawdzany % i == 0)
			{
				liczbaPierwsza = false;
				break;
			}
		}
		czyPierwsza[sprawdzany] = liczbaPierwsza;
	}

	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("naiveParallel", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void eratosthenesSequential(int m, int n) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	for (int i = 1; i <= n; i++)
		czyPierwsza[i] = true;

	for (int i = 2; i * i <= n; i++)
		if (czyPierwsza[i])
			for (int j = i * i; j <= n; j += i)
				czyPierwsza[j] = false;

	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesSequential", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void eratosthenesOddSequential(int m, int n) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	for (int i = 1; i <= n; i += 2) {
		czyPierwsza[i - 1] = false;
		czyPierwsza[i] = true;
	}
	czyPierwsza[2] = true;

	for (int i = 3; i * i <= n; i += 2)
		if (czyPierwsza[i])
			for (int j = i * i; j <= n; j += 2 * i)
				czyPierwsza[j] = false;

	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesOddSequential", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

//podejście funkcyjne
//proces działa na podzbiorze liczb wykreślających, operuje na całej tablicy wykreśleń
void eratosthenesFunctionParallel(int m, int n, int threads){
	double real_start, real_stop;
	real_start = omp_get_wtime();
	omp_set_num_threads(threads);
	bool* czyPierwsza = new bool[n + 1];
#pragma omp parallel
	{
#pragma omp for
		for (int i = 1; i <= n; i += 2) {
			czyPierwsza[i - 1] = false;
			czyPierwsza[i] = true;
		}
		czyPierwsza[2] = true;
		int n_sqrt = sqrt(n);
#pragma omp for schedule(dynamic, 1) 
		for (int i = 3; i <= n_sqrt; i += 2)
			if (czyPierwsza[i])
				for (int j = i * i; j <= n; j += 2 * i) {
					czyPierwsza[j] = false;
				}
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesFunctionParallel", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}


vector<int> eratosthenesOnlyOdd_returnVector(int m, int n) {
	vector<int> primes = {};
	int nn = (n - m + 1) / 2;
	bool* czyPierwsza = new bool[nn];
	for (int i = 0; i < nn; i++)
		czyPierwsza[i] = true;
	for (int i = 3; i * i <= n; i += 2)
	{
		int start = ((m + i - 1) / i) * i;
		if (start < i * i)
			start = i * i;
		start += (start % 2 == 0) ? i : 0;
		for (int j = start; j <= n; j += 2 * i)
			czyPierwsza[(j - m) / 2] = false;
	}
	for (int i = 0; i < nn; i++) {
		if (czyPierwsza[i] == true)
			primes.push_back(i * 2 + m + 1);
	}
	delete[] czyPierwsza;
	return primes;
}

void eratosthenesDomainParallel(int m, int n, int sliceSize, int threads) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	omp_set_num_threads(threads);
	bool* primes = new bool[n + 1];
	for (int i = 1; i <= n; i++) primes[i] = false;
	primes[2] = true;
#pragma omp parallel
	{
#pragma omp for 
		for (int from = 2; from <= n; from += sliceSize)
		{
			int to = from + sliceSize;
			if (to > n)
				to = n;
			vector<int> results = eratosthenesOnlyOdd_returnVector(from, to);
			for (int j = 0; j < results.size(); j++) {
				primes[results[j]] = true;
			}
		}
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	wypisz("eratosthenesDomainParallel", primes, m, n);
	delete[] primes;
}

int main()
{
	//naiveSequential(M, N);
	//naiveSequentialOdd(M, N);
	//naiveParallel(M, N, THREADS);
	//eratosthenesSequential(M, N);
	//eratosthenesOddSequential(M, N);
	//eratosthenesFunctionParallel(M, N, THREADS);
	eratosthenesDomainParallel(M, N, 128 * 1024, THREADS);
	return 0;
}