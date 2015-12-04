#ifndef __BARRIER_H__
#define __BARRIER_H__
#include <atomic>
#include <cstdlib>
#include <omp.h>

struct CoreCacheLine{
	std::atomic<int> signalLeft;
	std::atomic<int> signalRight;
	char padding[64-2*sizeof(std::atomic<int>)];
};
class Barrier{
private:
	CoreCacheLine *messages;
	size_t nThreads;

	inline void sendRight(size_t threadId){
		while(messages[threadId].signalRight.load(std::memory_order_relaxed) != 0);
		messages[threadId].signalRight.store(1, std::memory_order_relaxed);

		while(messages[threadId].signalRight.load(std::memory_order_relaxed) != 2);
		messages[threadId].signalRight.store(0, std::memory_order_relaxed);
	}

	inline void receiveLeft(size_t threadId){
		while(messages[threadId-1].signalRight.load(std::memory_order_relaxed) != 1);
		messages[threadId-1].signalRight.store(2, std::memory_order_relaxed);
	}

	inline void sendLeft(size_t threadId){
		int ex;
		do{
			ex = 0;
		}
		while(!messages[threadId].signalLeft.compare_exchange_weak(ex, 1, std::memory_order_relaxed));

		do{
			ex = 2;
		}
		while(!messages[threadId].signalLeft.compare_exchange_weak(ex, 0, std::memory_order_relaxed));
	}

	inline void receiveRight(size_t threadId){
		int ex;
		do{
			ex = 1;
		}
		while(!messages[threadId+1].signalLeft.compare_exchange_weak(ex, 2, std::memory_order_relaxed));
	}
public:
	Barrier(size_t n = 0){
		nThreads = n;
		messages = new CoreCacheLine[nThreads];
		for(size_t i=0; i<nThreads; ++i){
			messages[i].signalLeft.store(0);
			messages[i].signalRight.store(0);
		}
	}

	void barrier(size_t threadId){
		if (threadId == 0){
			sendRight(threadId);
		}
		else{
			receiveLeft(threadId);
			if (threadId != nThreads - 1){
				sendRight(threadId);
			}
		}

		if (threadId == nThreads-1){
			sendLeft(threadId);
		}
		else{
			receiveRight(threadId);
			if (threadId != 0){
				sendLeft(threadId);
			}
		}
	}

	inline void barrier(){
		barrier(omp_get_thread_num());
	}
};

class BarrierCounter{
private:
	std::atomic<int> entered;
	std::atomic<int> leaved;
	size_t nThreads;
public:
	BarrierCounter(size_t n){
		nThreads = n;
		entered.store(0);
		leaved.store(0);
	}
	
	inline void barrier(size_t threadId){
		while(leaved.load(std::memory_order_relaxed));
		
		entered.fetch_add(1, std::memory_order_relaxed);
		while(entered.load(std::memory_order_relaxed) < nThreads);

		leaved.fetch_add(1, std::memory_order_relaxed);
		while(leaved.load(std::memory_order_relaxed) < nThreads);

		int old = entered.fetch_add(-1, std::memory_order_relaxed);
		if (old == 1){
			leaved.store(0, std::memory_order_relaxed);
		}
	}
};
#endif