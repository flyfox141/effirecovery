#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#define MAXS 500
#define MAXN 200000

int g[ MAXN ][ MAXN ];
int map[ MAXN ];
int now;
int n, m;

void init(void){
	char s[ MAXS ];
	char T, T2[ MAXS ], T3[ MAXS ];
	int i, j, k;
	
	freopen("input.txt", "r", stdin);
	memset(map, -1, sizeof(map));
	now = 0;

	while(!feof(stdin)){
		fgets(s, MAXS, stdin);
		sscanf(s, "%c %s %s", &T, T2, T3);
		if (T == 'D'){
			j = 0;
			for(i=0; i<strlen(T2); i++){
				if (T2[i] >= '0' && T2[i] <= '9')
					j = j * 10 + T2[i] - '0';
				else
					break;
			}
			if (map[j] == -1)
				map[j] = now++;
			
			k = 0;
			for(i=0; i<strlen(T3); i++){
				if (T3[i] >= '0' && T3[i] <= '9')
					k = k * 10 + T3[i] - '0';
				else
					break;
			}
			if (map[k] == -1)
				map[k] = now++;

			if (g[map[j]][map[k]] == 0)
				m++;
			g[map[j]][map[k]] = 1;
			g[map[k]][map[j]] = 1;
		}
	}
	
	n = now;
	return;
}

void output(void){
	int i, j, k;
	freopen("graph_t.txt", "w", stdout);

	printf("%d %d\n", n, m);
	for(i=0; i<n; i++)
		for(j=i+1; j<n; j++){
			if (g[i][j] != 0){
				printf("%d %d %d\n", i, j, g[i][j]);
			}
		}
	
	return;
}


int main(void){
	init();
	output();
	return 0;
}
