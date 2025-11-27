#pragma once
#include <chrono>
#include <fstream>

class ChronoTimer
{
public:
	ChronoTimer(size_t init_reserve = 5UI64){
		start.reserve(init_reserve);
		counted_times.reserve(init_reserve);
	}
	// Save current time point as start
	inline void setStartTime()
	{
		start.push_back(std::chrono::steady_clock::now());
	}
	// Only duration is saved to array (from start(last ind - a) to steady_clock::now())
	void saveTimeNow(size_t a = 0Ui64)
	{
		if(a < start.size())
			counted_times.push_back(std::chrono::steady_clock::now() - *(start.end()-a-1Ui64));
	}
	// Get (from start to end) duration in sec
	inline long long getDiffInS(size_t a = 0Ui64) const
	{
		return std::chrono::duration_cast<std::chrono::seconds>(counted_times.at(a)).count();
	}
	// Get (from start to end) duration in milisec
	inline long long getDiffInMS(size_t a = 0Ui64) const
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(counted_times.at(a)).count();
	}
	// Get (from start to end) duration in nanosec
	inline long long getDiffInNS(size_t a = 0Ui64) const
	{
		return std::chrono::duration_cast<std::chrono::nanoseconds>(counted_times.at(a)).count();
	}
	// Will save counted_times in ns
	void saveTimesToFile(std::ofstream& sf, size_t indent=0)
	{
		unsigned i = 0Ui64;
		for (auto it = counted_times.begin(); it != counted_times.end(); it++)
		{
			sf << std::chrono::duration_cast<std::chrono::nanoseconds>(*it).count() << " ";
			if (i++ == indent-1Ui64) 
			{
				sf << "\n";
				i = 0Ui64;
			}
		}
	}

private:
	std::vector < std::chrono::steady_clock::time_point > start;
	std::vector<std::chrono::steady_clock::duration> counted_times;
	
};
