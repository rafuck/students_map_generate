#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <omp.h>

#pragma warning (disable : 4996)

double timer(){
	static clock_t start = clock();

	clock_t end = clock();
	double  ret = (double)(end - start)/CLOCKS_PER_SEC;
	
	start = end;
	
	return ret;
}

typedef enum{
	DirectionFirst = 0,
	Left,
	Right,
	Up,
	Down,
	DirectionLast
} Direction;

class IMap{
public:
	virtual size_t numRows() const = 0;
	virtual size_t numCols() const = 0;
};

struct Coord{
	union{
		int row, y, i;
	};
	union{
		int col, x, j;
	};

	Coord(int r, int c): row(r), col(c){}

	Coord(): row(0), col(0){}

	inline bool canStepTo(Direction direction, const IMap &map) const{
		switch(direction){
			case Left:  return x > 0;
			case Right: return x < map.numCols() - 1;
			case Up:    return y < map.numRows() - 1;
			case Down:  return y > 0;
			default:    return false;
		}
	}

	Coord stepTo(Direction direction) const{
		switch(direction){
			case Left:  return Coord(y, x-1);
			case Right: return Coord(y, x+1);
			case Up:    return Coord(y+1, x);
			case Down:  return Coord(y-1, x);
			default:    return Coord(x, y);
		}
	}
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

	inline Coord &operator[](size_t i) const{
		return this->at(i);
	}

	inline Coord &at(size_t i) const{
		return front[i];
	}

	inline size_t length() const{
		return count;
	}

	inline bool empty() const{
		return count == 0;
	}

	void clear(){
		count = 0;
	}
};

class Map: public IMap{
private:
	int *_mark[2];
	bool *_isWall;
	bool *_isPath;
	size_t N, M;
public:
	Map(size_t pN, size_t pM): N(pN), M(pM){
		_mark[0] = new int[M*N];
		_mark[1] = new int[M*N];
		_isWall  = new bool[M*N];
		_isPath  = new bool[M*N];

		for(size_t i=0; i<N*M; ++i){
			_mark[0][i] = _mark[1][i] = 0;
			_isWall[i] = true;
			_isPath[i] = false;
		}
	}

	~Map(){
		delete[] _mark[0];
		delete[] _mark[1];
		delete[] _isWall;
		delete[] _isPath;
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

	inline bool &isPath(const Coord &c) const{
		return _isPath[c.i*M + c.j];
	}

	void clearMarks(){
		for(size_t i=0; i<N*M; ++i){
			_mark[0][i] = _mark[1][i] = 0;
		}
	}

	void toHTML(const char *fname) const{
		FILE *fp = fopen(fname, "w");
		fprintf(fp, "<!DOCTYPE html><html><body> "
			"<style> "
			"html, body{background-color:#fff; padding:0; margin:0} "
			"table{border:1px solid #eee; border-right:0; border-top: 0} "
			"td{border:1px solid #eee; border-left:0; border-bottom:0; width: 5px; height: 5px} "
			"td.wall{background-color: #f0f0f0} "
			"td.path{background-color: #00ff00} "
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
				else if (_isPath[k]){
					fprintf(fp, "<td class=\"path\"></td>");
				}
				else{
					int cl0 = std::round(std::ceil(255.0*_mark[0][k]/max[0]));
					int cl1 = std::round(std::ceil(255.0*_mark[1][k]/max[1]));
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
			const Coord c = front->pop();
			for(int dir = DirectionFirst+1; dir != DirectionLast; ++dir){
				if (c.canStepTo(static_cast<Direction>(dir), map)){
					Coord cNew = c.stepTo(static_cast<Direction>(dir));
					if (!map.mark(id, cNew) && !map.isWall(cNew)){
						map.mark(id, cNew) = nStep;
						frontNew->push(cNew);
					}
				}	
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

	Coord meetCell(const Map &map) const{
		for(size_t i=0; i<front->length(); ++i){
			if (map.mark(0, front->at(i)) && map.mark(1, front->at(i))){
				return front->at(i);
			}
		}
		return Coord(0, 0);
	}

	void backward(Map &map, Coord c) const{
		for(int k=nStep; k>0; --k){
			for(int dir = DirectionFirst+1; dir != DirectionLast; ++dir){
				if (c.canStepTo(static_cast<Direction>(dir), map)){
					Coord cNew = c.stepTo(static_cast<Direction>(dir));
					if (map.mark(id, cNew) == k){
						map.isPath(cNew) = true;
						c = cNew;
						break;
					}
				}	
			}	
		}
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

		Coord meetCell = waves[(!isContinue[0]) ? 0 : 1].meetCell(map);
		waves[id].backward(map, meetCell);
	}

	double T = timer();
	printf("Finded    at %e s; %e s per iteration\n", T, T/iteration);
	return (iteration < maxIter);
}

int main(void){
	const size_t N = 256;
	const size_t M = 256;

	Map map(N, M);
	generate(map, M*N);
	map.toHTML("map.html");

	map.clearMarks();
	find(map, 4*(M+N));
	map.toHTML("map_waves.html");

	printf("Press enter to exit... ");
	getchar();
	
	return EXIT_SUCCESS;
}