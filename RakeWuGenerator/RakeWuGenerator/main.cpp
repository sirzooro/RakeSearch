// Генерация заданий для поиска "граблением"

# include <iostream>
# include <string>
# include "Square.h"
# include "Generator.h"

int main()
{
	int workunitsCount = 1854;

	string startFileName = "start_parameters.txt";
	string resultFileName = "result.txt";
	string checkpointFileName = "checkpoint.txt";
	string tempCheckpointFileName = "checkpoint_new.txt";

 	Generator wuGenerator;
	wuGenerator.GenerateWorkunits(startFileName, resultFileName, checkpointFileName, tempCheckpointFileName, workunitsCount);

	cout << "Press any key to exit ... " << endl;
	cin.get();
}
