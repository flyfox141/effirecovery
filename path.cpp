#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <time.h>
#include <vector>
#include <map>

using namespace std;

#define MAXTEST 1000
#define MAXT 10000
#define MAXN 4100
#define MAXPATH 4
#define MAXLABEL 1000
#define MILLION 1000000
#define U_short unsigned short
#define MAXLABEL2 4000
#define MAXNUM 214748364
#define MAXTEST2 1
#define TEST_QUERY 1

FILE *fp1, *fp2;


int n, m;

typedef struct{
	int from, to;
	int v, reverse;
	int next;
}Edge;

struct LABEL_1{
	int v, p;
	int f_hop, l_hop;
	vector <U_short> path;
	int lp;
	bool operator < (const LABEL_1 &a) const{
		return a.v < v;
	}
	LABEL_1 & operator = (const LABEL_1& a ){
		path = a.path;
//		memcpy(path, a.path, sizeof(U_short)*MAXN);
//		memcpy(seq, a.seq, sizeof(U_short)*MAXN);
		v = a.v; p = a.p;
		f_hop = a.f_hop; l_hop = a.l_hop;
		lp = a.lp;
		return *this;
	}
	/*
	LABEL_1(void){
		if ( path == NULL )
			exit( -1 );
		seq = new U_short[ MAXN ];
		if ( seq == NULL )
			exit( -1 );
			}*/
	~LABEL_1( void ){
//		delete path;
//		delete seq;
	}
};


struct LABEL_2{
	int v, p, cv, last;
	bool operator < (const LABEL_2 &a) const{
		return a.v < v;
	}
};

struct Dist_3{
	int v, p, e;
	bool operator < (const Dist_3 &a) const{
		return a.v < v;
	}
};

struct neighbor_4{
	int e, depth, ed;
	bool operator < (const neighbor_4 &a) const{
		if (a.depth > depth)
			return true;
		else if (a.depth < depth)
			return false;
		else if (a.ed < ed)
			return true;
		else
			return false;
	}
};

struct Dist_5{
	int v, p;
	bool operator < (const Dist_5 &a) const{
		return a.v > v;
	}
};

Edge T[ MAXT ];
int pT[ MAXN ];
int lt;
int n_label[ MAXN ];
int t_query[ MAXN ], lt_q = 0;;
int n_prune[ MAXN ];
int all_query = 0;
int times, times2;
int last_v, min_cost, origin;

int n2_label[ MAXN ];
int n2_prune[ MAXN ];


int path[ MAXN ][ MAXN ];
int l_path[ MAXN ];
bool children[ MAXN ][ MAXN ];
bool path2[ MAXN ][ MAXN ];
int depth[ MAXN ][ MAXN ];
neighbor_4 neighbor[ MAXN ][ MAXN ];
int l_neighbor[ MAXN ];
int t_label1, t_label2;

double minp = MAXNUM, maxp = -1;

int q_4[ MAXN ][ MAXN ];
int maxlabel = -1;
int origincv;
int ttt[MAXN][MAXN];


/*
LABEL_1 L_1[ MAXN ][ MAXLABEL ];
LABEL_2 L_2[ MAXN ][ MAXLABEL ];
*/
/*
LABEL_1 *L_1[ MAXN ];
LABEL_2 *L_2[ MAXN ];
*/
vector<LABEL_1> L_1[MAXN];
vector<LABEL_2> L_2[ MAXN ];


int lL_1[ MAXN ], lL_2[ MAXN ];
/*
map<U_short,U_short> ruv[ MAXN ][ MAXN ];
map<U_short,U_short> e_ruv[ MAXN ][ MAXN ];
*/


int dij_w[ MAXN ][ MAXN ];

void insert( int a, int b, int v ){
	T[lt].from = a;
	T[lt].to = b;
	T[lt].v = v;
	T[lt].next = pT[a];
	pT[a] = lt;
	lt++;
	
	T[lt].from = b;
	T[lt].to = a;
	T[lt].v = v;
	T[lt].next = pT[b];
	pT[b] = lt;

	T[lt - 1].reverse = lt;
	T[lt].reverse = lt - 1;
	lt++;
	
	return;
}

void init( void ){
	freopen("graph.in", "r", stdin);

	int i, j;
	int a, b, v;
	memset(pT, -1, sizeof(pT));
	scanf("%d %d\n", &n, &m);
	
	for(i=0; i<m; i++){
		scanf("%d %d %d\n", &a, &b, &v);
		insert(a, b, v);
	}

	srand((unsigned) time(NULL));

	fp1 = fopen("timeused.txt", "w");
	
	return;
}

bool find_uv(LABEL_1 t, int u, int v){
	vector<U_short>::iterator iter;

	for(iter=t.path.begin(); iter!=t.path.end(); iter++){
		if (v == -1){
			if (T[*iter].from == u || T[*iter].to == u)
				return true;
	        }
		else{
			if ((T[*iter].from == u && T[*iter].to == v)||
			    (T[*iter].from == v && T[*iter].to == u))
				return true;
		}
	}
	return false;
}

bool query_1( int v, int u, int length, U_short head, U_short tail ){
	int i, j, k;
	int t1 = 0, t2 = 0;
	int exist1[ MAXN ], le1 = 0, exist2[ MAXN ], le2 = 0;
	bool flag, flag1, flag2;
	int lp1, lp2;
	
	i = 0; j = 0;
	flag1 = true; flag2 = true; times = 0;
	while( i < lL_1[ v ] && j < lL_1[ u ] ){
		times++;
		if (L_1[ v ][ i ].p == L_1[ u ][ j ].p){
			if (L_1[ v ][ i ].v + L_1[ u ][ j ].v < length ) {

				lp1 = L_1[ v ][ i ].lp; lp2 = L_1[ u ][ j ].lp;
				
				if (le1 < 2 && flag1)
				{
					flag = true;
					for (k=0; k<le1; k++){
						if (v != L_1[ v ][ i ].p){
							if (exist1[ k ] ==
							    L_1[ v ][ i ].path.back() || T[exist1[ k ]].reverse == L_1[v][i].path.back())
								flag = false;
						}
						else
							if (exist1[ k ] ==
							    L_1[ u ][ j ].path.front() || T[exist1[k]].reverse == L_1[u][j].path.front())
								flag = false;
					}
					if (flag && !find_uv(L_1[u][j], v, -1)){
						if (v != L_1[ v ][ i ].p)
							exist1[ le1++ ] = L_1[ v ][ i ].path.back();
						else
							exist1[ le1++ ] = L_1[ u ][ j ].path.front();
					}
					if (v != L_1[ v ][ i ].p){
						if (head == L_1[ v ][ i ].path.back() || T[head].reverse == L_1[v][i].path.back())
							flag1 = false;
					}
					else
						if (head == L_1[ u ][ j ].path.front() || T[head].reverse == L_1[u][j].path.front())
							flag1 = false;
				}
				else
					flag1 = false;
				if (le1 >= 2) flag1 = false;

				if (le2 < 2 && flag2)
				{
					flag = true;
					for (k=0; k<le2; k++){
						if (exist2[ k ] ==
						    L_1[ u ][ j ].path.back() || T[exist2[k]].reverse == L_1[u][j].path.back())
							flag = false;	
					}
					if (flag && !find_uv(L_1[v][i],u,-1))
						exist2[ le2++ ] = L_1[ u ][ j ].path.back();
					if (tail == L_1[ u ][ j ].path.back() || T[tail].reverse == L_1[u][j].path.back())
						flag2 = false;
				}
				else
					flag2 = false;
				if (le2 >= 2) flag2 = false;

				if (!flag1 && !flag2)
					return false;
				

				}
			if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1)
			{
				if (L_1[ v ][ i ].v < L_1[ u ][ j ].v)
					i++;
				else
					j++;
			}
			else if (i == lL_1[ v ] - 1)
				j++;
			else
				i++;
		}
		else if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1){
			if (L_1[ v ][ i ].p < L_1[ u ][ j ].p)
				i++;
			else
				j++;
		}
		else if (i == lL_1[ v ] - 1)
			j++;
		else
			i++;
		
	}
	return true;
}

bool quick_prune(int v, int u, U_short head, U_short tail){
	int i, j, k;
	int exist1[ MAXN ], le1 = 0, exist2[ MAXN ], le2 = 0;
	bool flag, flag1, flag2;
	int lp1;

	
	i = 0; flag1 = true; flag2 = true;
	for(i=0; i<lL_1[v]; i++){
		if (L_1[v][i].p == u){
			lp1 = L_1[v][i].lp;
			if (lp1 <= 0) continue;
			if (le1 < 2 && flag1){
				flag = true;
				for(k=0; k<le1; k++){
					if (exist1[k] == L_1[v][i].path.front() || T[exist1[k]].reverse == L_1[v][i].path.front())
						flag = false;
					}
				if (flag)
					exist1[ le1++ ] = L_1[v][i].path.front();
				if (head == L_1[v][i].path.front() || T[head].reverse == L_1[v][i].path.front())
					flag1 = false;
			}
			else
				flag1 = false;
			if (le1 >= 2) flag1 = false;

			if (le2 < 2 && flag2){
				flag = true;
				for(k=0; k<le2; k++){
					if (exist2[k] == L_1[v][i].path.back() || T[exist2[k]].reverse == L_1[v][i].path.back())
						flag = false;
					}
				if (flag)
					exist2[ le2++ ] = L_1[v][i].path.back();
				if (head == L_1[v][i].path.back() || T[head].reverse == L_1[v][i].path.back())
					flag2 = false;
			}
			else
				flag2 = false;
			if (le2 >= 2) flag2 = false;
			if (!flag1 && !flag2)
				return false;
		}
	}
	return true;
}

void Pruned_Dij( int v ){

	int i, j, k;
//	p_queue pq;
	priority_queue <LABEL_1> pq;
	LABEL_1 temp, temp2;
	bool flag;

	temp.p = v;
	temp.v = 0;
	temp.f_hop = -1;
	temp.l_hop = v;
	temp.lp = 0;
	pq.push(temp);
	
	while (!pq.empty()){
		temp = pq.top();
		pq.pop();


		flag = true;

		
		if (temp.lp > 0 && !query_1(v, temp.p, temp.v, temp.path.front(), temp.path.back())){
			n_prune[ v ]++;
			continue;

		}
		n_label[ v ]++;
//		t_query[lt_q++] = times;
//		all_query += times;
		
		if (temp.lp>0 && !quick_prune(temp.p, v, temp.path.front(), temp.path.back()))
			continue;

		lt_q++;

		L_1[temp.p].push_back(temp);
		
//		L_1[ temp.p ][ lL_1[ temp.p ]] = temp;
		L_1[ temp.p ][ lL_1[ temp.p ]].p = v;
		lL_1[ temp.p ]++;
		if (lL_1[temp.p] > maxlabel) maxlabel = lL_1[temp.p];
		
		for (i=pT[ temp.p ]; i!=-1; i=T[i].next){
			if (T[ i ].to <= v)
				continue;

			if (find_uv(temp,T[ i ].to,-1) > 0)
				continue;

			if (T[i].to == 743)
				k = 3245;

			temp2.path = temp.path;
			temp2.p = T[ i ].to;
			temp2.v = temp.v + T[ i ].v;
			temp2.l_hop = temp.p;
			temp2.lp = temp.lp;
			if (temp.p == v)
				temp2.f_hop = temp2.p;
			else
				temp2.f_hop = temp.f_hop;
//			temp2.seq[ T[ i ].to ] = temp2.seq[ temp.p ] + 1;
//			temp2.path[ temp2.lp++ ] = i;
			temp2.path.push_back(i);
			temp2.lp = temp2.path.size();
			if (i == 6700)
				k = 234234;
			
			if (quick_prune(temp2.p, v, temp2.path.front(), temp2.path.back()))
				pq.push( temp2 );
		}
	}
        return;
}

void label_1( void ){
	int i, j;

	for (i=0; i<n; i++){
		Pruned_Dij(i);
	}
	return;
}


int query_2( int v, int u, int head ){
	int i, j, k;
	int lp1, lp2;
	bool flag = false, flag2 = false;

	int result = MAXNUM;
	times = 0;
	
	i = 0; j = 0; last_v = -1; min_cost = result, origin = result;
	while(i < lL_1[ v ] && j < lL_1[ u ] ){
		times++;
		if (L_1[ v ][ i ].p == L_1[ u ][ j ].p){
			if (L_1[ v ][ i ].v + L_1[ u ][ j ].v < min_cost)
				min_cost = L_1[ v ][ i ].v + L_1[ u ][ j ].v;

			if (L_1[ v ][ i ].v + L_1[ u ][ j ].v < origin){
				lp1 = L_1[ v ][ i ].lp; lp2 = L_1[ u ][ j ].lp;
				if (lp1 != 0 && lp2 != 0){
					if (L_1[ v ][ i ].path.back() ==
					    head ||  L_1[v ][i].path.back() == T[head].reverse){
						origin = L_1[v][i].v + L_1[u][j].v;
					}
				}
				else if (lp1 == 0 && L_1[ v ][ i ].p == v && lp2
					 != 0){
					if (L_1[ u ][ j ].path.front() == head 
					    || L_1[ u ][j].path.front() == T[head].reverse){
						origin = L_1[u][j].v;
					}
				}
			        else if (lp2 == 0 && L_1[ u ][ j ].p == v && lp1
					 != 0){
					if (L_1[ v ][ i ].path.back() == head 
					    || L_1[ v ][i].path.back() == T[head].reverse){
						origin = L_1[v][i].v;
					}
				}
			}

			if (origin > 100000)
				k = 234;
			
			if (L_1[ v ][ i ].v + L_1[ u ][ j ].v < result){
				lp1 = L_1[ v ][ i ].lp; lp2 = L_1[ u ][ j ].lp;
				if (lp1 != 0 && lp2 != 0){
					if (L_1[ v ][ i ].path.back() !=
					    head && L_1[v ][i].path.back() != T[head].reverse){
						if (!find_uv(L_1[ u ][ j ],T[head].to,T[head].from))
							result = L_1[v][i].v +	L_1[u][j].v;
						        last_v = L_1[v][i].path.back();
					}
				}
				else if (lp1 == 0 && L_1[ v ][ i ].p == v && lp2
					 != 0){
					if (L_1[ u ][ j ].path.front() != head &&
					    L_1[u][j].v < result && L_1[ u ][j].path.front() != T[head].reverse){
						result = L_1[u][j].v;
						last_v = L_1[u][j].path.front();
					}
					
				}
				else if (lp2 == 0 && L_1[ u ][ j ].p == v && lp1
					 != 0){
					if (L_1[ v ][ i ].path.back() != head &&
					    L_1[v][i].v < result && L_1[ v ][i].path.back() != T[head].reverse){
						result = L_1[v][i].v;
						last_v = L_1[v][i].path.front();
					}
					
				}
				
			}
			flag = false;
			if (find_uv(L_1[ u ][ j ],T[head].to,T[head].from))
				flag = true;
			
			if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1){
				if (L_1[ v ][ i ].v < L_1[ u ][ j ].v && !flag)
					i++;
				else
					j++;
			}
			else if (i == lL_1[ v ] - 1)
				j++;
			else
				i++;
		}
		else if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1){
			if (L_1[ v ][ i ].p < L_1[ u ][ j ].p)
				i++;
			else
				j++;
		}
		else if (i == lL_1[ v ] - 1)
			j++;
		else
			i++;
	}
	return result;
}

int query_4(int v, int u){
	int i;
	for(i=0; i<lL_1[u]; i++)
		if (L_1[u][i].p == v)
			return L_1[u][i].v;
	return MAXNUM;
}

void get_q_4( void ){
	int i, j;
	for(i=0; i<n; i++){
		for(j=i; j<n; j++){
			q_4[i][j] = query_4(i, j);
			q_4[j][i] = q_4[i][j];
		}
	}
	return;
}

bool query_3( int v, int u, int length ){
	int i, j, k, tt, t1, t2, a, b, c, d;
	
	i = 0; j = 0;
	times2 = 0;
	while( i < lL_2[ v ] && j < lL_2[ u ] ){
		times2++;
		if (L_2[ v ][ i ].p == L_2[ u ][ j ].p){
			a = L_2[v][i].v; if (v == L_2[v][i].p) a = MAXNUM;
			c = L_2[u][j].v; if (u == L_2[u][j].p) c = MAXNUM;
/*
			b = query_4(u, L_2[u][j].p);
			d = query_4(v, L_2[v][i].p);
*/
			b = q_4[u][L_2[u][j].p];
			d = q_4[v][L_2[v][i].p];
			
			t1 = (a + b);
			t2 = (c + d);
			if (a == MAXNUM)
				tt = t2;
			else if (c == MAXNUM)
				tt = t1;
			else
				tt = max(t1, t2);
			
			if (tt <= length || tt >= MAXNUM) 
					return false;
			if (i < lL_2[ v ] - 1 && j < lL_2[ u ] - 1)
			{
				if (L_2[ v ][ i ].v < L_2[ u ][ j ].v)
					i++;
				else
					j++;
			}
			else if (i == lL_2[ v ] - 1)
				j++;
			else
				i++;
		}
		else if (i < lL_2[ v ] - 1 && j < lL_2[ u ] - 1){
			if (L_2[ v ][ i ].p > L_2[ u ][ j ].p)
				i++;
			else
				j++;
		}
		else if (i == lL_2[ v ] - 1)
			j++;
		else
			i++;
		
	}
	return true;
}


void Pruned_Sgp( int v )
{
	int i, j, k;
	priority_queue <LABEL_2> pq;
	LABEL_2 temp, temp2;
	bool flag;
	int dist[ MAXN ], tt;
	int ta, tb, last, ncv;

	memset( dist, -1, sizeof(dist));
	temp.p = v;
	temp.v = 0;
	temp.cv = 0;
	temp.last = -1; 
	pq.push(temp);
	
	while (!pq.empty()){
		temp = pq.top();
		pq.pop();

		if (dist[ temp.p ] != -1)
			continue;
		if (!query_3(temp.p, v, temp.v))
		{
			n2_prune[ v ]++;
			continue;
		}

		n2_label[ v ]++;

		L_2[ temp.p ].push_back(temp);
		L_2[ temp.p ][ lL_2[ temp.p ]].p = v;
		lL_2[ temp.p ]++;

		dist[ temp.p ] = temp.v;
		
		for (i=pT[ temp.p ]; i!=-1; i=T[i].next){
			if (T[ i ].to >= v)
				continue;
			
			ta = query_2(T[i].to, v, i);
			tb = T[i].v + temp.v;
/*
			if (temp.last == i)
			tb = temp.v - T[i].v;
			if (ta < origin || (ta >= MAXNUM && temp.p != v))
				ta = -1;
			if ((tb < temp.cv + T[i].v && temp.p != v) || tb >=
			    MAXNUM || (tb == temp.cv + T[i].v && temp.p == v) )
			    tb = -1;*/
			if (ta > tb){
				tt = ta;
				last = last_v;
				ncv = origin;
			}
			else{
				tt = tb;
				last = i;
				ncv = temp.cv + T[ i ].v;
			}

			if (tt == -1 || tt >= MAXNUM)
				continue;

			if (dist[ T[ i ].to ] == -1 || dist[ T[ i ].to ] > tt )
			{
				temp2 = temp;
				temp2.p = T[ i ].to;
				temp2.v = tt;
				temp2.cv = ncv;
				temp2.last = last;
				pq.push( temp2 );
			}

		}
	}
        return;
}

void label_2( void ){
	int i;
	for(i=n-1; i>=0; i--){
//		printf("%d\n", i);
		Pruned_Sgp(i);
		
	}
	return;
}

int get_z(int x, int y){
	int i, j;
	int f1, e1, m1;
	int f2, e2, m2;

	e1 = l_path[ x ] - 1;
	e2 = l_path[ y ] - 1;
	f1 = 0; f2 = 0;

	if (e1 < e2){
		while(f1 < e1){
			m1 = (f1 + e1) / 2;
			if (!path2[ y ][path[ x ][ m1 ]]){
				f1 = m1 + 1;
			}
			else{
				e1 = m1;
			}
		}
		if (!path2[y][path[x][f1]])
			return -1;
		else
			return path[x][f1];

	}
	else if (e1 >= e2){
		while(f2 < e2){
			m2 = (f2 + e2) / 2;
			if (!path2[ x ][path[ y ][ m2 ]]){
				f2 = m2 + 1;
			}
			else{
				e2 = m2;
			}
		}
		if (!path2[x][path[y][f2]])
			return -1;
		else
			return path[y][f2];
	}
	return -1;
}

int query_5( int v, int u ){
	int i, j, k, tt, t1, t2, a, b, c, d;
	i = 0; j = 0; int cv;
	times2 = 0; cv = 0;
	int result = MAXNUM;
	while( i < lL_2[ v ] && j < lL_2[ u ] ){
		times2++;
		if (L_2[ v ][ i ].p == L_2[ u ][ j ].p){
			a = L_2[v][i].v; //if (v == L_2[v][i].p) a = MAXNUM;
			c = L_2[u][j].v; //if (u == L_2[u][j].p) c = MAXNUM;

			b = q_4[u][L_2[u][j].p];
			d = q_4[v][L_2[v][i].p];

//			b = query_4(u, L_2[u][j].p);
//			d = query_4(v, L_2[v][i].p);
			t1 = (a + b);
			t2 = (c + d);
			if (a == MAXNUM){
				tt = t2;
				cv = L_2[u][j].cv + d;
			}
			else if (c == MAXNUM){
				tt = t1;
				cv = L_2[v][i].cv + b;
			}
			
			else{
				if (t1 > t2){
					tt = t1;
					cv = L_2[v][i].cv + b;
				}
				else{
					tt = t2;
					cv = L_2[u][j].cv + d;
				}

			}
			
			if (tt <= result ){
				result = tt;
				origincv = cv;
			}
			if (i < lL_2[ v ] - 1 && j < lL_2[ u ] - 1)
			{
				if (L_2[ v ][ i ].v < L_2[ u ][ j ].v)
					i++;
				else
					j++;
			}
			else if (i == lL_2[ v ] - 1)
				j++;
			else
				i++;
		}
		else if (i < lL_2[ v ] - 1 && j < lL_2[ u ] - 1){
			if (L_2[ v ][ i ].p > L_2[ u ][ j ].p)
				i++;
			else
				j++;
		}
		else if (i == lL_2[ v ] - 1)
			j++;
		else
			i++;
		
	}
	return result;
}
/*
void dij_memo(int p){
	int i, j, k;
	static	int dist[ MAXN ];
	static int father[ MAXN ];
	priority_queue <Dist_3> pq;
	Dist_3 temp, temp2;
	int tt;

	memset(dist, -1, sizeof(dist));
	memset(father, -1, sizeof(father));
	memset(path2, 0, sizeof(path2));
	memset(depth, 0, sizeof(depth));
	temp.p = p;
	temp.v = 0;
	pq.push(temp);

	while (!pq.empty()){
		temp = pq.top();
		pq.pop();

		if (dist[ temp.p ] != -1)
			continue;

		father[temp.p] = T[temp.e].from;
		if (temp.p == p) father[temp.p] = -1;

		k = temp.p; j = 0;
		while(k != -1){
			path[temp.p][j++] = k;
			path2[ temp.p ][k] = true;
			children[k][temp.p] = true;
			k = father[k];
		}
		l_path[ temp.p ] = j;

		dist[ temp.p ] = temp.v;
		
		for(i=pT[temp.p]; i!=-1; i=T[i].next){
			tt = temp.v + T[i].v;
			if (dist[ T[i].to ] == -1 || tt < dist[ T[i].to ]){
				temp2.p = T[i].to;
				temp2.v = tt;
				temp2.e = i;
				pq.push(temp2);
			}
				
		}
	}

static	int ed[ MAXN ][ MAXN ];
static	bool fg[ MAXN ][ MAXN ];
	priority_queue <neighbor_4> pq2[ MAXN ];
	neighbor_4 temp3;
	Dist_5 temp4;
	priority_queue <Dist_5> pq3;
	int x_temp[ MAXN ];

	memset(fg, 0, sizeof(fg));
	       
	for(i=0; i<lt; i++) {
		if (fg[T[i].from][T[i].to])
			continue;
		ed[T[i].from][T[i].to] = dist[T[i].from] + dist[T[i].to] +
	T[i].v;
		ed[T[i].to][T[i].from] = ed[T[i].from][T[i].to];
		fg[T[i].from][T[i].to] = true;
		fg[T[i].to][T[i].from] = true;
		if (i == 10)
			k  = 234;
		tt = get_z(T[i].from, T[i].to);
		depth[ T[i].from ][ T[i].to ] = dist[ tt ];
		depth[ T[i].to ][T[i].from] = dist[tt];
		temp3.e = i;
		temp3.depth = dist[tt];
		temp3.ed = ed[T[i].to][T[i].from];
		pq2[T[i].from].push(temp3);
		temp3.e = T[i].reverse;
		pq2[T[i].to].push(temp3);
	}

	memset(l_neighbor, 0, sizeof(l_neighbor));
	for(i=0; i<n; i++){
		while(!pq2[i].empty()){
			temp3 = pq2[i].top();
			pq2[i].pop();
			neighbor[i][l_neighbor[i]++] = temp3;
		}
	}

	
		ruv[T[i].from][p][  ] = 0;
		ruv[T[i].to][p][ MAXT ] = 0;
		}

	for(i=0; i<n; i++){
		if (dist[i] != -1){
			temp4.v = dist[i];
			temp4.p = i;
			pq3.push( temp4 );
		}
	}

	int u, v, x, y, e;
	int result = MAXNUM;

	while(!pq3.empty()){
		temp4 = pq3.top();
		pq3.pop();

		u = temp4.p;
		if (u == 0)
			k = 345;

		
		for(i=0; i<n; i++){
			for(x_temp[i]=0; x_temp[i]<l_neighbor[i]; x_temp[i]++){
				if (neighbor[i][x_temp[i]].depth < dist[u])
					break;
			}
		}

		e = neighbor[u][x_temp[u]].e;
		v = T[e].to;
		x_temp[u]++;

		result = MAXNUM;
		
		for(i=0; i<n; i++){
			if (!children[u][i])
				continue;
			if (x_temp[i] < l_neighbor[i]){
				y = T[neighbor[i][x_temp[i]].e].to;
				if (ed[i][y] < result)
					result = ed[i][y];
			}
		}
		ruv[u][p][v] = result - dist[u];
		e_ruv[u][p][v] = e;
	}		
	tt = 0;
	return;
}
*/
int dij_naive(int u, int v, int head){
	int dist[ MAXN ];
	bool flag[ MAXN ];
	priority_queue <Dist_3> pq;
	Dist_3 temp, temp2;
	int i, to, k;

	int result = MAXNUM;

	temp.v = 0; temp.p = u; temp.e = -1;
	memset(dist, -1, sizeof(dist));
	memset(flag, 0, sizeof(flag));
	
	pq.push(temp);

	while(!pq.empty()){
		temp = pq.top();
		pq.pop();

		if (flag[temp.p])
			continue;

		dist[temp.p] = temp.v;
		flag[temp.p] = true;

		for(i=pT[temp.p]; i!=-1; i=T[i].next){
			if (i == head || i == T[head].reverse)
				continue;
			to = T[i].to;
			if (to == 563 && u == 0)
				k = 234;
			if (dist[ to ] == -1 || dist[ to ] > temp.v + T[i].v){
				dist[to] = temp.v + T[i].v;
				temp2.v = temp.v + T[i].v;
				temp2.p = to;
				temp2.e = i;
				pq.push(temp2);
			}
		}
	}

	result = dist[v];

	return result;
}

void dij_worst(int u){
	int dist[ MAXN ];
	bool flag[ MAXN ];
	priority_queue <Dist_3> pq;
	Dist_3 temp, temp2;
	int i, to, w = -1, tt, k;

	int result = MAXNUM;

	temp.v = 0; temp.p = u; temp.e = -1;
	memset(dist, -1, sizeof(dist));
	memset(flag, 0, sizeof(flag));
	
	pq.push(temp);

	while(!pq.empty()){
		temp = pq.top();
		pq.pop();

		if (flag[temp.p])
			continue;

		dist[temp.p] = temp.v;
		flag[temp.p] = true;

		for(i=pT[temp.p]; i!=-1; i=T[i].next){
			to = T[i].to;
			if (dist[ to ] == -1 || dist[ to ] > temp.v + T[i].v){
				dist[to] = temp.v + T[i].v;
				temp2.v = temp.v + T[i].v;
				temp2.p = to;
				temp2.e = i;
				w = query_2(to, u, i);
//				w = dij_naive(to, u, i);
				if (w == -1)
					w = MAXNUM;
				if (u == 13 && to == 0)
					k = 23423;
				if (w > dij_w[to][u]) dij_w[to][u] = w;
				if (dij_w[temp.p][u]+T[i].v > dij_w[to][u])
					dij_w[to][u] = dij_w[temp.p][u]+T[i].v;
				pq.push(temp2);
			}
		}
	}

	return;
}




double test( void ){

	int i, j, r, u, v, w;
	int nr[ MAXN ], l_nr;
	
        clock_t t1, t2;
	double p = 0, p2 = 0;

	minp = MAXNUM; maxp = -1;

	for(i=0; i<MAXTEST; i++){
		u = rand() % n;
		v = (rand() % (n-u)) + u;
		l_nr = 0;
		for(j=pT[u]; j!=-1; j=T[j].next)
			nr[ l_nr++ ] = j;
		w = rand() % l_nr;

		t1 = clock();
		r = query_2(u, v, w);
		t2 = clock();

		p2 = (double)(t2 - t1) / CLOCKS_PER_SEC;
		if (p2 < minp) minp = p2;
		if (p2 > maxp) maxp = p2;

		p += p2;
	}

	p /= MAXTEST;



	return p;
}

int query_6( int v, int u){
	int i, j, k;
	int lp1, lp2;
	bool flag = false;

	int result = MAXNUM;
	times = 0;
	
	i = 0; j = 0; 
	while(i < lL_1[ v ] && j < lL_1[ u ] ){
		if (L_1[ v ][ i ].p == L_1[ u ][ j ].p){
			if (L_1[ v ][ i ].v + L_1[ u ][ j ].v < result)
				result = L_1[ v ][ i ].v + L_1[ u ][ j ].v;
			
			if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1){
				if (L_1[ v ][ i ].v < L_1[ u ][ j ].v)
					i++;
				else
					j++;
			}
			else if (i == lL_1[ v ] - 1)
				j++;
			else
				i++;
		}
		else if (i < lL_1[ v ] - 1 && j < lL_1[ u ] - 1){
			if (L_1[ v ][ i ].p < L_1[ u ][ j ].p)
				i++;
			else
				j++;
		}
		else if (i == lL_1[ v ] - 1)
			j++;
		else
			i++;
	}
	return result;
}


void test2(void){
	
        int t1, t2;
	int i, j;
	int v, cv, sp, sp_w;
	double pp3;
	priority_queue<int, vector<int>, greater<int> > p, q, r, s;

	int total = 0;
	
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++)
		{
			v=query_5(i, j);
			sp = query_6(i, j);
			cv = origincv;
			sp_w = dij_w[i][j];
			if (v < MAXNUM && v != -1 && sp_w < MAXNUM){
				if (v <= sp_w){
					p.push(v);
					q.push(sp);
					r.push(cv);
					s.push(sp_w);
				}
				else{
					p.push(sp_w);
					q.push(sp);
					r.push(sp);
					s.push(sp_w);
				}
				total++;
			}
		}
	}

	FILE *fp1 = fopen("length_v.txt", "w");
	t1 = 0;
	while (!p.empty()){
		v = p.top();
		p.pop();
		t1++;
		fprintf(fp1, "%lf %d\n", (double)t1/(double)total, v);
	}

	FILE *fp2 = fopen("length_sp.txt", "w");
	t1 = 0;
	while (!q.empty()){
		sp = q.top();
		q.pop();
		t1++;
		fprintf(fp2, "%lf %d\n", (double)t1/(double)total, sp);
	}
	fclose(fp2);

	fp2 = fopen("length_cv.txt", "w");
	t1 = 0;
	while (!r.empty()){
		cv = r.top();
		r.pop();
		t1++;
		fprintf(fp2, "%lf %d\n", (double)t1/(double)total, cv);
	}

	fp2 = fopen("length_sp_2.txt", "w");
	t1 = 0; int last = 0;
	while (!s.empty()){
		cv = s.top();
		s.pop();
		t1++;
		last = cv;
		fprintf(fp2, "%lf %d\n", (double)t1/(double)total, cv);
	}
	
	return;
}

void find_exist(void){
	int i, j, p;
	int count = 0;
	vector<U_short>::iterator iter;
	map<int,int>::iterator iter2;

	for(i=0; i<n; i++){
		map<int, int> existtimes[ MAXN ];
		for(j=0; j<lL_1[i]; j++){
			for(iter=L_1[i][j].path.begin(); iter!=L_1[i][j].path.end(); iter++){
				p = L_1[i][j].p;
				iter2 = existtimes[p].find(*iter);
				if (iter2 != existtimes[p].end()){
					existtimes[p][*iter]++;
					if (existtimes[p][*iter] >= 3){
						printf("i=%d, p=%d, from=%d, to=%d\n", i, p, T[*iter].from, T[*iter].to);
						count++;
					}
						

				}
				else
					existtimes[p][*iter] = 1;
			}
		}
	}
	printf("count=%d\n", count);
	return;
}


int main( void ){

	int i, j, k;
        clock_t t1, t2;
	double it1, qt1, pp3, pp4;

	init();

	label_1();

	find_exist();
	memset(dij_w, -1, sizeof(dij_w));
/*
	get_q_4();
	
	label_2();

	for(i=0; i<n; i++){
		printf("%d\n", i);
		dij_worst(i);
	}
		

	test2();
*/	
	return 0;
}
