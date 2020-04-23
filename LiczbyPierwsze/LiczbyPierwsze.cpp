#include<math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
using namespace std;

#define M 2
#define N 10000000
#define LICZBA_PROCESOROW_FIZYCZNYCH 4
#define LICZBA_PROCESOROW_LOGICZNYCH 8
#define WYPISZ_PIERWSZE 1


void simpleSequential(int,int);
void eratosthenesSequential(int,int);
void wypisz(string, vector<int>);


int main()
{
	simpleSequential(M,N);
	eratosthenesSequential(M, N);
	
	return 0;
}

void simpleSequential(int m, int n) {
	vector<int> pierwsze{2};
	vector<int> wynik{};
	if (m == 2) wynik.push_back(2);
	for (int sprawdzany = 3; sprawdzany < n; sprawdzany++)
	{
		bool liczbaZlozona = false;
		for (int i = 0; pierwsze[i] <= sqrt(sprawdzany) && i < pierwsze.size(); i++)
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

vector<int> eratosthenesOnlyOddSequential(int m, int n)
{
	int onlyOdd = (n - 1) / 2;
	vector<int> result;
	bool* czyPierwsza = new bool[onlyOdd + 1];
	for (int i = 0; i <= onlyOdd; i++)
		czyPierwsza[i] = true;

	for (int i = 2; i * i <= n; i++)
		if (czyPierwsza[i])
			for (int j = i * i; j <= n; j += i)
				czyPierwsza[j] = false;

	for (int i = m; i <= onlyOdd; i++)
		if (czyPierwsza[i] == true)
			result.push_back(i);
	return result;
}

void wypisz(string tekst, vector<int> result) {
	cout << tekst << "  znalazl " << result.size() << " liczb pierwszych\n";
	for (int i = 0; i < 100 && i < result.size(); i++) {
		printf("%d ", result[i]);
		if ((i + 1) % 10 == 0) printf("\n");
	}
	printf("\n");
}