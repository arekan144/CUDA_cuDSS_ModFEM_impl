#pragma once
#ifndef CHRONOTIMER
#define CHRONOTIMER
#include <chrono>
class ChronoTimer
{
public:
	ChronoTimer() :start(), end() {};
	inline void setStartTime()
	{
		start = std::chrono::steady_clock::now();
	}
	inline void setEndTime()
	{
		end = std::chrono::steady_clock::now();
	}
	
	inline long long getDiffInS() const
	{
		return std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
	}
	
	inline long long getDiffInMS() const
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	}
	
private:
	std::chrono::steady_clock::time_point start;
	std::chrono::steady_clock::time_point end;
};

#endif // !CHRONOTIMER