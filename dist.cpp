#include <cstdio>
#include <cstdlib>

using namespace std;

int main( void ){
	int i, j, now;
	double frac, last;

	freopen("length_cv.txt", "r", stdin);
	freopen("length_cv_dist.txt", "w", stdout);
	now = -1; last = 0;
	while(!feof(stdin)){
		scanf("%lf %d\n", &frac, &i);
		if (i != now){
			if (now != -1)
				printf("%lf %d\n", frac, now);
			now = i;
			last = frac;
		}			
	}
	if (now < 214748364)
		printf("%lf %d\n", frac, now);

	freopen("length_v.txt", "r", stdin);
	freopen("length_v_dist.txt", "w", stdout);
	now = -1; last = 0;
	while(!feof(stdin)){
		scanf("%lf %d\n", &frac, &i);
		if (i != now){
			if (now != -1)
				printf("%lf %d\n", frac, now);
			now = i;
			last = frac;
		}			
	}
	if (now < 214748364)
	printf("%lf %d\n", frac, now);

	freopen("length_sp.txt", "r", stdin);
	freopen("length_sp_dist.txt", "w", stdout);
	now = -1; last = 0;
	while(!feof(stdin)){
		scanf("%lf %d\n", &frac, &i);
		if (i != now){
			if (now != -1)
				printf("%lf %d\n", frac, now);
			now = i;
			last = frac;
		}			
	}
	if (now < 214748364)
	printf("%lf %d\n", frac, now);

	freopen("length_sp_2.txt", "r", stdin);
	freopen("length_sp_2_dist.txt", "w", stdout);
	now = -1; last = 0;
	while(!feof(stdin)){
		scanf("%lf %d\n", &frac, &i);
		if (i != now){
			if (now != -1)
				printf("%lf %d\n", frac, now);
			now = i;
			last = frac;
		}			
	}
	if (now < 214748364)
	printf("%lf %d\n", frac, now);
	
	return 0;
}

