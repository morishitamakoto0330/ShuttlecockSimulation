#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>



int main()
{
	std::ifstream ifs("../python_work/_data_1213.txt");
	std::string str, buf;
	std::vector<int> v;
	
	if(ifs.fail()) {
		std::cerr << "failed to read file." << std::endl;
		return -1;
	}
	
	while(getline(ifs, str)) {
		std::stringstream _str;
		_str << str;
		
		while(getline(_str, buf, ',')) {
			v.push_back(std::stoi(buf));
		}
	}

	for(int i = 0; i < v.size(); i++) {
		std::cout << v[i] << std::endl;
	}
	
	return 0;
}



