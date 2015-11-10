#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <omp.h>

double timer(){
	static clock_t start = clock();

	clock_t end = clock();
	double  ret = (double)(end - start)/CLOCKS_PER_SEC;
	
	start = end;
	
	return ret;
}

struct Coord{
	union{
		int row, y, i;
	};
	union{
		int col, x, j;
	};

	Coord(int r, int c): row(r), col(c){}

	Coord(): row(0), col(0){}
};

class Front{
private:
	Coord *front;
	size_t count, size;
public:
	Front(size_t N = 0): front(NULL), count(0), size(0){
		reserve(N);
	}

	void reserve(size_t N){
		if (size >= N) return;
		if (front){
			delete[] front;
		}

		front = new Coord[N];
		count = 0;
		size  = N;
	}

	~Front(){
		if (front){
			delete[] front;
		}
	}

	void push(Coord p){
		front[count++] = p;
	}

	Coord &pop(){
		return (count) ? front[--count] : front[0];
	}

	Coord &operator[](size_t i) const{
		return this->at(i);
	}

	inline Coord &at(size_t i) const{
		return front[i];
	}

	size_t length() const{
		return count;
	}

	bool empty() const{
		return count == 0;
	}

	void clear(){
		count = 0;
	}
};

class Map{
private:
	int *_mark[2];
	bool *_isWall;
	size_t N, M;
public:
	Map(size_t pN, size_t pM): N(pN), M(pM){
		_mark[0] = new int[M*N];
		_mark[1] = new int[M*N];
		_isWall  = new bool[M*N];

		for(size_t i=0; i<N*M; ++i){
			_mark[0][i] = _mark[1][i] = 0;
			_isWall[i] = true;
		}
	}

	~Map(){
		delete[] _mark[0];
		delete[] _mark[1];
		delete[] _isWall;
	}


	inline size_t numRows() const{
		return N;
	}

	inline size_t numCols() const{
		return M;
	}

	inline int &mark(int id, const Coord &c) const{
		return _mark[id][c.i*M + c.j];
	}

	inline bool &isWall(const Coord &c) const{
		return _isWall[c.i*M + c.j];
	}

	inline int &mark(int id, int i, int j) const{
		return _mark[id][i*M + j];
	}

	inline bool &isWall(int i, int j) const{
		return _isWall[i*M + j];
	}

	void clearMarks(){
		for(size_t i=0; i<N*M; ++i){
			_mark[0][i] = _mark[1][i] = 0;
		}
	}

	void toHTML(const char *fname) const{
		FILE *fp = fopen(fname, "w");
		fprintf(fp, "<!DOCTYPE html><html><body>"
			"<style>"
			"table{border:1px solid #eee; border-right:0; border-top: 0} "
			"td{border:1px solid #eee; border-left:0; border-bottom:0; width: 5px; height: 5px} "
			"td.wall{background-color: #f0f0f0} "
			"</style>"
			"<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\">");
		
		int max[2] = {0, 0};
		for(size_t i=0; i<M*N; ++i){
			if (max[0] < _mark[0][i]) max[0] = _mark[0][i];
			if (max[1] < _mark[1][i]) max[1] = _mark[1][i];
		}

		size_t k = 0;
		for(size_t i=0; i<N; ++i){
			fprintf(fp, "<tr>");
			for(size_t j=0; j<M; ++j){
				if (_isWall[k]){
					fprintf(fp, "<td class=\"wall\"></td>");
				}
				else{
					int cl0 = round(ceil(255.0*_mark[0][k]/max[0]));
					int cl1 = round(ceil(255.0*_mark[1][k]/max[1]));
					fprintf(fp, "<td style=\"background-color:rgb(%d, %d, 255)\"></td>", cl0, cl1);
				}
				k++;
			}
			fprintf(fp, "</tr>");
		}
		fprintf(fp, "</table></body></html>");
		fclose(fp);
	}
};

class Wave{
private:
	Front *front, *frontNew;
	int nStep, id;
	void destroy(){
		if (front){
			delete front; front = NULL;
		}
		
		if (frontNew){
			delete frontNew; frontNew = NULL;
		}
	}
public:
	Wave(): nStep(0), front(NULL), frontNew(NULL){}

	~Wave(){
		destroy();
	}

	void init(int pid, Coord c, const Map &map){
		id = pid;

		size_t N = 4*(map.numRows()+map.numCols());
		destroy();

		front    = new Front(N);
		frontNew = new Front(N);

		nStep = 1;
		front->push(c);

		map.mark(id, c) = 1;
	}

	void step(Map &map){
		frontNew->clear();
		nStep++;

		while(!front->empty()){
			Coord c = front->pop();
			if (c.i > 0 && !map.mark(id, c.i-1, c.j) && !map.isWall(c.i-1, c.j)){
				map.mark(id, c.i-1, c.j) = nStep;
				frontNew->push(Coord(c.i-1, c.j));
			}
			if (c.i < map.numRows()-1 && !map.mark(id, c.i+1, c.j) && !map.isWall(c.i+1, c.j)){
				map.mark(id, c.i+1, c.j) = nStep;
				frontNew->push(Coord(c.i+1, c.j));
			}
			if (c.j > 0 && !map.mark(id, c.i, c.j-1) && !map.isWall(c.i, c.j-1)){
				map.mark(id, c.i, c.j-1) = nStep;
				frontNew->push(Coord(c.i, c.j-1));
			}
			if (c.j < map.numCols()-1 && !map.mark(id, c.i, c.j+1) && !map.isWall(c.i, c.j+1)){
				map.mark(id, c.i, c.j+1) = nStep;
				frontNew->push(Coord(c.i, c.j+1));
			}
		}

		Front *f = front;
		front = frontNew;
		frontNew = f;
	}

	bool isEnd(Map &map) const{
		for(size_t i=0; i<front->length(); ++i){
			if (map.mark(0, front->at(i)) && map.mark(1, front->at(i))){
				return true;
			}
		}
		return false;
	}
};

class Warm{
private:
	Coord c;
	int id, nStep;
public:
	Warm():id(0), nStep(0){}

	void init(int pid, Coord pc, const Map &map){
		c    = pc;
		id   = pid;
		nStep = 1;

		if (id < 2) map.mark(id, c) = 1;
		map.isWall(c) = false;
	}

	void step(Map &map){
		int delta = (rand() % 2) ? 1 : -1;
		if (rand()%2){
			if (c.i == 0 && delta == -1) delta = 1;
			else if (c.i == map.numRows() - 1 && delta == 1) delta = -1;

			c.i += delta;
		}
		else{
			if (c.j == 0 && delta == -1) delta = 1;
			else if (c.j == map.numCols() - 1 && delta == 1) delta = -1;

			c.j += delta;
		}
		nStep++;

		if (id < 2 && !map.mark(id, c)){
			map.mark(id, c) = nStep;
		}
		map.isWall(c) = false;
	}

	bool isEnd(Map &map) const{
		return (map.mark(0, c) && map.mark(1, c));
	}
};

bool generate(Map &map, int maxIter){
	srand(time(NULL));
	Warm warms[4];
	warms[0].init(0, Coord(0, 0), map);
	warms[1].init(1, Coord(map.numRows()-1, map.numCols()-1), map);
	warms[2].init(2, Coord(map.numRows()-1, 0), map);
	warms[3].init(3, Coord(0, map.numCols()-1), map);

	int iteration = 0;
	bool isContinue = true;
	timer();
	#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		while(isContinue){
			warms[id].step(map);
			#pragma omp barrier
			#pragma omp master
			{
				iteration++;
				isContinue = iteration < maxIter && !warms[0].isEnd(map) && !warms[1].isEnd(map);
			}
			#pragma omp barrier
		}
	}
	
	double T = timer();
	printf("Generated at %e s; %e s per iteration\n", T, T/iteration);
	return (iteration < maxIter);
}

bool find(Map &map, int maxIter){
	Wave waves[2];
	waves[0].init(0, Coord(0, 0), map);
	waves[1].init(1, Coord(map.numRows()-1, map.numCols()-1), map);

	int iteration = 0;
	bool isContinue[2] = {true, true};
	timer();
	#pragma omp parallel num_threads(2)
	{
		int id = omp_get_thread_num();
		while(iteration < maxIter && isContinue[0] && isContinue[1]){
			waves[id].step(map);
			#pragma omp barrier
			#pragma omp master
			{
				iteration++;
			}
			isContinue[id] = !waves[id].isEnd(map);
			#pragma omp barrier
		}
	}
	
	double T = timer();
	printf("Finded    at %e s; %e s per iteration\n", T, T/iteration);
	return (iteration < maxIter);
}

int main(void){
	const size_t N = 512;
	const size_t M = 512;
	
	Map map(N, M);
	generate(map, M*N);
	map.toHTML("map.html");

	map.clearMarks();
	find(map, M+N);
	map.toHTML("map_waves.html");

	printf("Press enter to exit... ");
	getchar();
	
	return EXIT_SUCCESS;
}