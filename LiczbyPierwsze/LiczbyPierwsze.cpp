#include<math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <omp.h>
using namespace std;

#define M 2
#define N 1000000000
#define LICZBA_PROCESOROW_FIZYCZNYCH 4
#define LICZBA_PROCESOROW_LOGICZNYCH 8
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


void naiveSequential(int m, int n) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	czyPierwsza[2] = true;
	int suma = 0;
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
		wypisz("naiveSequential", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void eratosthenesSequential(int m, int n)
{
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

void eratosthenesOnlyOddSequential(int m, int n)
{
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
		wypisz("eratosthenesOnlyOddSequential", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

void naiveParallel(int m, int n, int threads) {
	double real_start, real_stop;
	real_start = omp_get_wtime();
	bool* czyPierwsza = new bool[n + 1];
	czyPierwsza[2] = true;
	int suma = 0;
	int mm = m % 2 == 0 ? m + 1 : m;
	omp_set_num_threads(threads);
#pragma omp parallel for //schedule(dynamic, 1)
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

//podejście funkcyjne
//proces działa na podzbiorze liczb wykreślających, operuje na całej tablicy wykreśleń
//brak na wejściu liczb pierwszych od 2 do sqrt(N) - proces korzysta ze wszystkich liczb nieparzystych z zakresu

void eratosthenesOnlyOddParallel(int m, int n, int threads)
{
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
				for (int j = i * i; j <= n; j += 2 * i)
					czyPierwsza[j] = false;
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesOnlyOddParallel", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

//podejście funkcyjne
//najpierw funkcja znajduje liczby od 2 do sqrt(n), potem na podstawie utworzonej robi skreślenia
vector<int> eratosthenesOnlyOddParallel_returnVector(int m, int n)
{
	bool* czyPierwsza = new bool[n + 1];
	vector<int> primes = {};
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
				for (int j = i * i; j <= n; j += 2 * i)
					czyPierwsza[j] = false;
	}
	for (int i = m; i < n; i++)
		if (czyPierwsza[i] == true)
			primes.push_back(i);
	delete[] czyPierwsza;
	return primes;
}
//c.d.
void eratosthenesOnlyOddParallelWithVectorOfPrimes(int m, int n, int num_threads)
{
	double real_start, real_stop;
	real_start = omp_get_wtime();
	omp_set_num_threads(num_threads);
	int n_sqrt = sqrt(n);
	vector<int> primes = eratosthenesOnlyOddParallel_returnVector(2, n_sqrt);
	bool* czyPierwsza = new bool[n + 1];
#pragma omp parallel
	{
#pragma omp for
		for (int i = 1; i <= n; i += 2) {
			czyPierwsza[i - 1] = false;
			czyPierwsza[i] = true;
		}
		czyPierwsza[2] = true;
		int size = primes.size();
#pragma omp for schedule(dynamic, 1)
		for (int i = 1; i < size; i++) {
			czyPierwsza[primes[i]] = true;
			for (int j = primes[i] * primes[i]; j <= n; j += 2 * primes[i])
			{
				czyPierwsza[j] = false;
				//printf("%d ",j);
			}
		}
	}
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesOnlyOddParallelWithVectorOfPrimes", czyPierwsza, m, n);
	}
	delete[] czyPierwsza;
}

vector<int> eratosthenesOnlyOdd_returnVector(int m, int n)
{
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
		{
			czyPierwsza[(j - m) / 2] = false;
		}
	}
	for (int i = 0; i < nn; i++) {
		if (czyPierwsza[i] == true)
			primes.push_back(i * 2 + m + 1);
	}
	delete[] czyPierwsza;
	return primes;
}

//podejście domenowe 
//podział przeszukiwanego zakresu liczb po
void eratosthenesBlok(int m, int n, int sliceSize)
{
	double real_start, real_stop;
	real_start = omp_get_wtime();
	volatile bool* primes = new bool[n + 1];
	for (int i = 1; i <= n; i++) primes[i] = false;
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
	real_stop = omp_get_wtime();
	printf("Suma czasow przetwarzania wszystkich watkow wynosi %f sekund\n", ((double)(real_stop - real_start)));
	wypisz("eratosthenesBlok", (bool*)primes, m, n);
	delete[] primes;
}

void eratosthenesBlokParallel(int m, int n, int sliceSize, int threads)
{

	double real_start, real_stop;
	real_start = omp_get_wtime();
	omp_set_num_threads(threads);
	volatile bool* primes = new bool[n + 1];
	for (int i = 1; i <= n; i++) primes[i] = false;
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
	wypisz("eratosthenesBlokParallel", (bool*)primes, m, n);
	delete[] primes;
}

int main()
{
	//naiveSequential(M, N);
	//naiveParallel(M, N, LICZBA_PROCESOROW_LOGICZNYCH);
	//eratosthenesSequential(M, N);
	eratosthenesOnlyOddSequential(M, N);
	//eratosthenesOnlyOddParallel(M, N, LICZBA_PROCESOROW_LOGICZNYCH);
	//eratosthenesOnlyOddParallelWithVectorOfPrimes(M, N, LICZBA_PROCESOROW_LOGICZNYCH);
	//eratosthenesBlok(M, N, 128 * 1024);
	//eratosthenesBlokParallel(M, N, 128 * 1024, LICZBA_PROCESOROW_LOGICZNYCH);
	return 0;
}