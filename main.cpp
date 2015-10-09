#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <omp.h>

struct Cell{
	int mark[2];
	bool isWall;

	Cell(){
		mark[0] = mark[1] = 0;
		isWall = true;
	}
};

class Map{
private:
	Cell *map;
	size_t N, M;
public:
	Map(size_t pN, size_t pM){
		N = pN;
		M = pM;
		map = new Cell[pN*pM];
	}

	~Map(){
		delete[] map;
	}


	size_t numRows() const{
		return N;
	}

	size_t numCols() const{
		return M;
	}

	Cell* row(size_t i) const{
		return map+M*i;
	}

	Cell *cell(size_t i, size_t j) const{
		return map + M*i + j;
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
		Cell *cell = map;
		int max[2] = {0, 0};
		for(size_t i=0; i<M*N; ++i){
			if (max[0] < cell[i].mark[0]) max[0] = cell[i].mark[0];
			if (max[1] < cell[i].mark[1]) max[1] = cell[i].mark[1];
		}

		for(size_t i=0; i<N; ++i){
			fprintf(fp, "<tr>");
			for(size_t j=0; j<M; ++j){
				if (cell->isWall){
					fprintf(fp, "<td class=\"wall\"></td>");
				}
				else{
					int cl0 = round(ceil(255.0*cell->mark[0]/max[0]));
					int cl1 = round(ceil(255.0*cell->mark[1]/max[1]));
					fprintf(fp, "<td style=\"background-color:rgb(%d, %d, 255)\"></td>", cl0, cl1);
				}

				cell++;
			}
			fprintf(fp, "</tr>");
		}
		fprintf(fp, "</table></body></html>");
		fclose(fp);
	}
};

class Warm{
private:
	Cell *cell;
	int i, j, id, nStep;
public:
	Warm(){
		cell = NULL;
		id = i = j = 0;
		nStep = 0;
	}

	void init(int pid, int pi, int pj, const Map &map){
		i    = pi;
		j    = pj;
		id   = pid;
		nStep = 1;
		cell = map.cell(i, j);
		if (id < 2) cell->mark[id] = 1;
		cell->isWall   = false;
	}

	void step(Map &map){
		int delta = (rand() % 2) ? 1 : -1;
		if (rand()%2){
			if (i == 0 && delta == -1) delta = 1;
			else if (i == map.numRows() - 1 && delta == 1) delta = -1;

			i += delta;
		}
		else{
			if (j == 0 && delta == -1) delta = 1;
			else if (j == map.numCols() - 1 && delta == 1) delta = -1;

			j += delta;
		}
		nStep++;

		printf("%d: %d, %d -> %d\n", id, i, j, nStep);

		cell = map.cell(i, j);
		if (id < 2 && !cell->mark[id]) cell->mark[id] = nStep;
		cell->isWall   = false;
	}

	bool isEnd() const{
		return (cell && cell->mark[0] && cell->mark[1]);
	}
};

bool generate(Map &map, int maxIter){
	srand(time(NULL));
	Warm warms[4];
	warms[0].init(0, 0, 0, map);
	warms[1].init(1, map.numRows()-1, map.numCols()-1, map);
	warms[2].init(2, map.numRows()-1, 0, map);
	warms[3].init(3, 0, map.numCols()-1, map);
	int iteration = 0;
	#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		while(iteration < maxIter && !warms[0].isEnd() && !warms[1].isEnd()){
			warms[id].step(map);
			#pragma omp master
			{
				iteration++;
			}
			#pragma omp barrier
		}
	}
	return (iteration < maxIter);
}

int main(void){
	const size_t N = 256;
	const size_t M = 256;
	
	Map map(N, M);
	generate(map, M*N);
	
	map.toHTML("map.html");
	
	return EXIT_SUCCESS;
}