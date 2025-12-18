#pragma once
#include <chrono>
#include <vector>
#include <fstream>
#ifdef _DEBUG
#include <iostream>
#endif
class ChronoTimer
{
public:
	ChronoTimer(size_t init_reserve = 5llu){
		start.reserve(init_reserve);
		counted_times.reserve(init_reserve);
	}
	// Save current time point as start
	inline void setStartTime()
	{
#ifdef _DEBUG
		std::cout << "Starting time\n";
#endif
		start.push_back(std::chrono::steady_clock::now());
	}
	// Only duration is saved to array (from start(last ind - a) to steady_clock::now())
	void saveTimeNow(size_t a = 0llu)
	{
#ifdef _DEBUG
		std::cout << "Saving time\n";
#endif
		if(a < start.size())
			counted_times.push_back(std::chrono::steady_clock::now() - *(start.end() - a - 1llu));
	}
	// Get (from start to end) duration in sec
	inline long long getDiffInS(size_t a = 0llu) const
	{
		return std::chrono::duration_cast<std::chrono::seconds>(counted_times.at(a)).count();
	}
	// Get (from start to end) duration in milisec
	inline long long getDiffInMS(size_t a = 0llu) const
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(counted_times.at(a)).count();
	}
	// Get (from start to end) duration in nanosec
	inline long long getDiffInNS(size_t a = 0llu) const
	{
		return std::chrono::duration_cast<std::chrono::nanoseconds>(counted_times.at(a)).count();
	}
	// Will save counted_times in ns
	void saveTimesToFile(std::ofstream& sf, size_t indent = 0llu)
	{
		size_t i = 0llu;
		for (auto it = counted_times.begin(); it != counted_times.end(); it++)
		{
			sf << std::chrono::duration_cast<std::chrono::nanoseconds>(*it).count() << " ";
			if (i++ == indent - 1llu) 
			{
				sf << "\n";
				i = 0llu;
			}
		}
	}

private:
	std::vector<std::chrono::steady_clock::time_point> start;
	std::vector<std::chrono::steady_clock::duration> counted_times;
	
};

