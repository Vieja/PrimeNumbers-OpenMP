#include<math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
using namespace std;

#define M 2
#define N 100000000
#define LICZBA_PROCESOROW_FIZYCZNYCH 4
#define LICZBA_PROCESOROW_LOGICZNYCH 8
#define WYPISZ_PIERWSZE 1

void wypisz(string tekst, vector<int> result) {
	cout << tekst << "  znalazl " << result.size() << " liczb pierwszych\n";
	for (unsigned int i = 0; i < 100 && i < result.size(); i++) {
		printf("%d ", result[i]);
		if ((i + 1) % 10 == 0) printf("\n");
	}
	printf("\n");
}

void simpleSequential(int m, int n) {
	vector<int> pierwsze{2};
	vector<int> wynik{};
	if (m == 2) wynik.push_back(2);
	for (int sprawdzany = 3; sprawdzany < n; sprawdzany++)
	{
		bool liczbaZlozona = false;
		for (unsigned int i = 0; pierwsze[i] <= sqrt(sprawdzany) && i < pierwsze.size(); i++)
		{
			if (sprawdzany % pierwsze[i] == 0)
			{
				liczbaZlozona = true;
				break;
			}
		}
		if (!liczbaZlozona)
		{
			pierwsze.push_back(sprawdzany);
			if (sprawdzany >= m)
			{
				wynik.push_back(sprawdzany);
			}
		}
	}
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("simpleSequential",wynik);
	}
} 

void eratosthenesSequential(int m, int n)
{
	vector<int> wynik{};
	bool* czyPierwsza = new bool[n + 1];
	for (int i = 0; i <= n; i++)
		czyPierwsza[i] = true;

	for (int i = 2; i * i <= n; i++)
		if (czyPierwsza[i])
			for (int j = i * i; j <= n; j += i)
				czyPierwsza[j] = false;

	for (int i = m; i <= n; i++)
		if (czyPierwsza[i] == true)
			wynik.push_back(i);

	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesSequential", wynik);
	}
}

void eratosthenesOnlyOddSequential(int m, int n)
{
	vector<int> wynik{};
	if (m == 2) wynik.push_back(2);
	bool* czyPierwsza = new bool[n + 1];
	for (int i = 0; i <= n; i++)
		czyPierwsza[i] = true;

	for (int i = 3; i * i <= n; i+= 2)
		if (czyPierwsza[i])
			for (int j = i * i; j <= n; j += 2*i)
				czyPierwsza[j] = false;

	int start = (m % 2 == 0) ? m+1 : m;
	for (int i = start; i <= n; i += 2)
		if (czyPierwsza[i] == true)
			wynik.push_back(i);
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesOnlyOddSequential", wynik);
	}
}

void eratosthenesOnlyOddSequential2(int m, int n)
{
	vector<int> wynik{};
	int memorySize = (n - 1) / 2;
	if (m == 2) wynik.push_back(2);
	bool* czyPierwsza = new bool[memorySize + 1];
	for (int i = 0; i <= memorySize; i++)
		czyPierwsza[i] = true;

	for (int i = 3; i * i <= n; i += 2)
		if (czyPierwsza[i/2])
			for (int j = i * i; j <= n; j += 2 * i)
				czyPierwsza[j/2] = false;

	//int start = (m % 2 == 0) ? m + 1 : m;
	for (int i = 1; i <= memorySize; i++)
		if (czyPierwsza[i] == true)
			wynik.push_back((i*2)+1);
	if (WYPISZ_PIERWSZE == 1) {
		wypisz("eratosthenesOnlyOddSequential2", wynik);
	}
}


int main()
{
	//simpleSequential(M,N);
	eratosthenesSequential(M, N);
	eratosthenesOnlyOddSequential(M, N);
	eratosthenesOnlyOddSequential2(M, N);
	//sprawdzic co szybsze
	return 0;
}