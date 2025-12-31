#pragma once
#include <map>
#include <string>
#include <iostream>
/* matrix file -f
matrix file type -l (default "text") [(t)ext, (b)inary]
program mode -t (default "eigen+cudss") [(E)igen+cudss,(C)usolve+cudss,(S)ave txt as binary]
matrix type -m (default "general") [(g)eneral, (s)ymmetric, (h)ermitian, symmetric-positive-defined (spd), hermitian-positive-defined (hpd)]
evaluation -e (default yes) [(y)es, (n)o, (l)oad solution]
solution file -s
*/
class ArgumentParser
{
public: 
	ArgumentParser() = delete;
	ArgumentParser(int& argc, char** argv) { // throws std::string
		for (int i = 1; i < argc; i++)
		{
			if (argv[i][0] == '-' && argv[i][2] == '\0') {
				if (i + 1 == argc || (argv[i + 1][0] == '-' && argv[i + 1][2] == '\0'))
					throw std::string("Error: Expected value after argument: ") + std::string(argv[i]);
				auto pos = inputs.find(argv[i][1]);
				if (pos != inputs.end()) 
					pos->second.first = std::string(argv[i + 1]);
				else
					throw std::string("Error: Wrong argument: ") + std::string(argv[i]);
				// we already checked the following argument/value, thrown an error or saved the value
				i++;
			}
			else
			{
				bool found = false;
				if (i <= 1 && inputs['f'].first.empty()) {
					inputs['f'].first = std::string(argv[i]);
					found = true;
				}
				else if (i > 1 && argv[i][1] == '\0') {
					found = true;
					switch (argv[i][0]) {
					case 'y': case 'n': case 'Y': case 'N':
						inputs['e'].first = std::string(argv[i]);
						break;
					case 't': case 'b': case 'T': case 'B':
						inputs['l'].first = std::string(argv[i]);
						break;
					case 'E': case 'S': case 'C':
						inputs['t'].first = std::string(argv[i]);
						break;
					case 'g': case 'h': case 's':
						inputs['m'].first = std::string(argv[i]);
						break;
					default:
						found = false;
						break;
					}
				}
				if(!found)
					throw std::string("Error: Wrong value: ") + std::string(argv[i]);

			}
		}
		parse();
	}
	
	std::pair<std::string,short> getArgument(const char opt) {
		auto pos = inputs.find(opt);
		if (pos != inputs.end())
			return pos->second;
		return {"", (short)(-1)};
	}

	void print()
	{
		//std::map<const char, std::pair<std::string, short>>::iterator it;
		for (auto it = inputs.begin(); it != inputs.end(); it++) {
			std::cout << it->first << ": "
				<< it->second.first << ' '
				<< it->second.second << '\n';
		}
	}

private:

	void parse() {
		if (inputs['f'].first.empty())
			throw std::string("Error: you need to provide file name.");

		if (inputs['l'].first == "t" || inputs['l'].first == "txt" || inputs['l'].first == "text")
			inputs['l'].second = 0;
		else if (inputs['l'].first == "b" || inputs['l'].first == "bin" || inputs['l'].first == "binary")
			inputs['l'].second = 1;
		else 
			throw std::string("Error: Wrong file type.");
		
		if (inputs['m'].first == "g" || inputs['m'].first == "gen" || inputs['m'].first == "general")
			inputs['m'].second = 0;
		else if (inputs['m'].first == "s" || inputs['m'].first == "sym" || inputs['m'].first == "symmetric")
			inputs['m'].second = 1;
		else if (inputs['m'].first == "h" || inputs['m'].first == "her" || inputs['m'].first == "hermitian")
			inputs['m'].second = 2;
		else if (inputs['m'].first == "spd" || inputs['m'].first == "symmetric-pd" || inputs['m'].first == "symmetric-positive-definite")
			inputs['m'].second = 3;
		else if (inputs['m'].first == "hpd" || inputs['m'].first == "hermitian-pd" || inputs['m'].first == "hermitian-positive-defined")
			inputs['m'].second = 4;
		else
			throw std::string("Error: Wrong matrix type.");

		if (inputs['t'].first == "e" || inputs['t'].first == "E" || inputs['t'].first == "eigen" || inputs['t'].first == "eigen+cudss")
			inputs['t'].second = 0;
		else if (inputs['t'].first == "c" || inputs['t'].first == "C" || inputs['t'].first == "cusolve" || inputs['t'].first == "cusolve+cudss")
			inputs['t'].second = 1;
		else if (inputs['t'].first == "s" || inputs['t'].first == "S" || inputs['t'].first == "save")
			inputs['t'].second = 2;
		else if (inputs['t'].first == "p" || inputs['t'].first == "P" || inputs['t'].first == "Pattern")
			inputs['t'].second = 3;
		else
			throw std::string("Error: Wrong program type.");

		if (inputs['e'].first == "y" || inputs['e'].first == "Y" || inputs['e'].first == "yes") 
			inputs['e'].second = 1;
		else if (inputs['e'].first == "n" || inputs['e'].first == "N" || inputs['e'].first == "no") 
			inputs['e'].second = 0;
		else if (inputs['e'].first == "s" || inputs['e'].first == "S" || inputs['e'].first == "save")
			inputs['e'].second = 2;
		else if (inputs['e'].first == "l" || inputs['e'].first == "L" || inputs['e'].first == "load") {
			inputs['e'].second = 3;
			if(inputs['s'].first.empty())
				throw std::string("Error: evaluation option set to load, but no file was specified.");
		}
		else
			throw std::string("Error: Wrong evaluation option type.");
		
		if (inputs['e'].second != 3 && !inputs['s'].first.empty())
			std::cout << "Warning: Solution filed provieded, with option " << inputs['e'].first << '\n';

	}
	std::map<const char, std::pair<std::string, unsigned short>> inputs{ 
		{'f',{std::string(""),0}},   // matrix file -f
		{'l',{std::string("t"),0}},  // matrix file type -l (default "text") [(t)ext, (b)inary]
		{'t',{std::string("E"),0}},  // program mode -t (default "eigen+cudss") [(E)igen+cudss,(C)usolve+cudss,(S)ave txt as binary]
		{'m',{std::string("g"),0}},  // matrix type -m (default "general") [(g)eneral, (s)ymmetric, (h)ermitian, symmetric-positive-defined (spd), hermitian-positive-defined (hpd)]
		{'e',{std::string("y"),1}},  // evaluation -e (default yes) [(y)es, (n)o, (l)oad solution]
		{'s',{std::string(""),0}},	 // solution file
	};
	
};
