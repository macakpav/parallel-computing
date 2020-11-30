/**MPI parallelization of serial game of life code.
 * Pavel Mačák (pavel.macak@fs.cvut.cz)
 *
 * Serial implementation of the game of life for the course Parallel Computing.
 * Emil Loevbak (emil.loevbak@cs.kuleuven.be)
 * First implementation: November 2019
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <cmath>
#include <tuple>

int const globalBufferLength = 50;

int getPID(const int row, const int col, const int pRows, const int pCols)
{
	// Get ID of processor on row row and column col in processor grid (pRows x pCols)

	if (col == pCols)
		return getPID(row, 0, pRows, pCols);
	if (col == -1)
		return getPID(row, pCols - 1, pRows, pCols);
	if (row == pRows)
		return getPID(0, col, pRows, pCols);
	if (row == -1)
		return getPID(pRows - 1, col, pRows, pCols);

	if (row < 0 || col < 0 || row >= pRows || col >= pCols)
	{
		return -1; // means no neighboor
	}
	return row * pRows + col;
}

std::tuple<int, int> getPcoords(const int pID, const int pRows, const int pCols)
{
	// Get row and column of processor pID
	int pRow = 0 + pID / pCols;
	int pCol = 0 + pID - (pID / pCols) * (pCols);
	return {pRow, pCol};
}

std::tuple<int, int> getRowRange(const int pID, const int pRows, const int pCols, const int rowsPerP, const int colsPerP)
{
	// Get row range of processor pID
	auto [pRow, pCol] = getPcoords(pID, pRows, pCols);
	int pFirstRow = pRow * rowsPerP;
	int pLastRow = pFirstRow + rowsPerP - 1;
	return {pFirstRow, pLastRow};
}

std::tuple<int, int> getRowRange(const int pRow, const int pCol, const int pRows, const int pCols, const int rowsPerP, const int colsPerP)
{
	// Get row range of processor at coordinates pRow, pCol

	int pFirstRow = pRow * rowsPerP;
	int pLastRow = pFirstRow + rowsPerP - 1;
	return {pFirstRow, pLastRow};
}

std::tuple<int, int> getColRange(const int pID, const int pRows, const int pCols, const int rowsPerP, const int colsPerP)
{
	// Get column range of processor pID
	auto [pRow, pCol] = getPcoords(pID, pRows, pCols);
	int pFirstCol = pCol * colsPerP;
	int pLastCol = pFirstCol + colsPerP - 1;
	return {pFirstCol, pLastCol};
}

std::tuple<int, int> getColRange(const int pRow, const int pCol, const int pRows, const int pCols, const int rowsPerP, const int colsPerP)
{
	// Get column range of processor at coordinates pRow, pCol

	int pFirstCol = pCol * colsPerP;
	int pLastCol = pFirstCol + colsPerP - 1;
	return {pFirstCol, pLastCol};
}

void initializeBoard(std::vector<std::vector<bool>> &board)
{
	int deadCellMultiplyer = 2;
	srand(time(0));
	for (auto &col : board)
	{
		for (auto element : col)
		{
			element = (rand() % (deadCellMultiplyer + 1) == 0);
		}
	}
}

void updateBoard(std::vector<std::vector<bool>> &board)
{
	const size_t rows = board.size();
	const size_t cols = board[0].size();
	std::vector<std::vector<int>> liveNeighbors(rows, std::vector<int>(cols, 0));

	//Count live neighbors
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			if (board[i][j])
			{
				for (int di = -1; di <= 1; ++di)
				{
					for (int dj = -1; dj <= 1; ++dj)
					{
						//Periodic boundary conditions
						liveNeighbors[(i + di + rows) % rows][(j + dj + cols) % cols]++;
					}
				}
				liveNeighbors[i][j]--; //Correction so that a cell does not concider itself as a live neighbor
			}
		}
	}

	//Update board
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			board[i][j] = ((liveNeighbors[i][j] == 3) || (board[i][j] && liveNeighbors[i][j] == 2));
		}
	}
}

void writeBoardToFile(std::vector<std::vector<bool>> &board, size_t firstRow, size_t lastRow, size_t firstCol, size_t lastCol, std::string fileName, int iteration, uint processID)
{
	//Open file
	std::ofstream outputFile(fileName + "_" + std::to_string(iteration) + "_" + std::to_string(processID) + ".gol");
	//Write metadata
	outputFile << std::to_string(firstRow) << " " << std::to_string(lastRow) << std::endl;
	outputFile << std::to_string(firstCol) << " " << std::to_string(lastCol) << std::endl;
	//Write data
	// std::ostream_iterator<bool> outputIterator(outputFile, "\t");
	for (size_t i = 1; i < board.size()-1; ++i){
		for (size_t j = 1; j < board[0].size()-1; j++)
		{
			outputFile << std::to_string(board[i][j]) << "\t";
		}
		outputFile << std::endl;
	}
	// for (size_t i = 1; i <= board.size(); ++i)
	// {
	// 	copy(++board[i].begin(), --board[i].end(), outputIterator);
	// 	outputFile << std::endl;
	// }

	//Close file
	outputFile.close();
}

std::string setUpProgram(size_t rows, size_t cols, int iteration_gap, int iterations, int processes)
{
	//Generate progam name based on current time, all threads should use the same name!
	time_t rawtime;
	struct tm *timeInfo;
	char buffer[globalBufferLength];
	time(&rawtime);
	timeInfo = localtime(&rawtime);
	strftime(buffer, sizeof(buffer), "%Y-%m-%d-%H-%M-%S", timeInfo);
	std::string programName(buffer);

	//Generate main file
	std::ofstream outputFile(programName + ".gol");
	outputFile << std::to_string(rows) << " " << std::to_string(cols) << " " << std::to_string(iteration_gap) << " " << std::to_string(iterations) << " " << std::to_string(processes) << std::endl;
	outputFile.close();

	return programName;
}

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		std::cout << "This program should be called with four arguments! \nThese should be, the total number of rows; the total number of columns; the gap between saved iterations and the total number of iterations, in that order." << std::endl;
		return 1;
	}
	size_t noRows, noCols;
	int iteration_gap, iterations;
	try
	{
		noRows = atoi(argv[1]);
		noCols = atoi(argv[2]);
		iteration_gap = atoi(argv[3]);
		iterations = atoi(argv[4]);
	}
	catch (std::exception const &exc)
	{
		std::cout << "One or more program arguments are invalid!" << std::endl;
		return 1;
	}
	// Initialize the MPI environment. The two arguments to MPI Init are not
	// currently used by MPI implementations, but are there in case future
	// implementations might need the arguments.
	MPI_Init(NULL, NULL);
	MPI_Status Stat;

	// Get the number of processes
	int processes;
	MPI_Comm_size(MPI_COMM_WORLD, &processes);

	// Get the rank of the process
	int myID;
	MPI_Comm_rank(MPI_COMM_WORLD, &myID);

	// size_t firstRow = 0, lastRow = noRows - 1, firstCol = 0, lastCol = noCols - 1;
	// All threads use the same name
	std::string programName = setUpProgram(noRows, noCols, iteration_gap, iterations, processes);

	int usedProcesses = (int)std::sqrt(processes) * (int)std::sqrt(processes);
	if (myID == 0)
	{
		printf("Game of Life has started...\n");
		printf("INFO: Number of tasks = %d. Only using %d tasks.\n", processes, usedProcesses);
	}
	processes = usedProcesses;

	//2D-square partitioning
	int procCols = (int)std::sqrt(processes);
	int procRows = procCols;

	// assuming noCols and noRows is divisible by procCols and procRows
	int rowsPerProc = noRows / procRows;
	int colsPerProc = noCols / procCols;

	//Build board separately for each thread
	auto [myRow, myCol] = getPcoords(myID, procRows, procCols);
	auto [myFirstRow, myLastRow] = getRowRange(myID, procRows, procCols, rowsPerProc, colsPerProc);
	auto [myFirstCol, myLastCol] = getColRange(myID, procRows, procCols, rowsPerProc, colsPerProc);

	int northNB = getPID(myRow - 1, myCol, procRows, procCols);
	int southNB = getPID(myRow + 1, myCol, procRows, procCols);
	int westNB = getPID(myRow, myCol - 1, procRows, procCols);
	int eastNB = getPID(myRow, myCol + 1, procRows, procCols);
	MPI_Barrier(MPI_COMM_WORLD);
#ifndef NDEBUG
	std::cout << "Processor " << myID << " will initialize its board." << std::endl;
#endif

	std::vector<std::vector<bool>> myBoard((myLastRow - myFirstRow + 1 + 2), std::vector<bool>(myLastCol - myFirstCol + 1 + 2));
	initializeBoard(myBoard);
#ifndef NDEBUG
	std::cout << "Processor " << myID << " finished initializing its board." << std::endl;
#endif
	MPI_Barrier(MPI_COMM_WORLD);

	//Do iteration
	bool westExchange[rowsPerProc];
	bool eastExchange[rowsPerProc];
	bool northExchange[colsPerProc];
	bool southExchange[colsPerProc];
#ifndef NDEBUG
	std::cout << "Processor " << myID << " will write its board to file." << std::endl;
#endif
	writeBoardToFile(myBoard, myFirstRow, myLastRow, myFirstCol, myLastCol, programName, 0, myID);
#ifndef NDEBUG
	std::cout << "Processor " << myID << " finished writing to file." << std::endl;
#endif
	for (int i = 1; i <= iterations; ++i)
	{
#ifndef NDEBUG
		std::cout << "Processor " << myID << " has board of size " << myBoard.size() << "x" << myBoard[0].size() << std::endl;
#endif
		if (myCol % 2 == 0)
		{
			for (int i = 0; i < rowsPerProc; i++)
			{
				westExchange[i] = myBoard[i][colsPerProc];
				eastExchange[i] = myBoard[i][1];
			}
#ifndef NDEBUG
			std::cout << "Processor " << myID << " has his exchange data ready." << std::endl;
#endif
			MPI_Send(westExchange, rowsPerProc, MPI_CXX_BOOL, westNB, 1, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << westNB << std::endl;
#endif
			MPI_Send(eastExchange, rowsPerProc, MPI_CXX_BOOL, eastNB, 2, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << eastNB << std::endl;
#endif
			MPI_Recv(eastExchange, rowsPerProc, MPI_CXX_BOOL, eastNB, 3, MPI_COMM_WORLD, &Stat);
			MPI_Recv(westExchange, rowsPerProc, MPI_CXX_BOOL, westNB, 4, MPI_COMM_WORLD, &Stat);
			for (int i = 0; i < rowsPerProc; i++)
			{
				myBoard[i][colsPerProc + 1] = eastExchange[i];
				myBoard[i][0] = westExchange[i];
			}
		}
		else
		{
#ifndef NDEBUG
			std::cout << "Processor " << myID << " waiting to recieve." << std::endl;
#endif
			MPI_Recv(eastExchange, rowsPerProc, MPI_CXX_BOOL, eastNB, 1, MPI_COMM_WORLD, &Stat);
			MPI_Recv(westExchange, rowsPerProc, MPI_CXX_BOOL, westNB, 2, MPI_COMM_WORLD, &Stat);
			for (int i = 0; i < rowsPerProc; i++)
			{
				myBoard[i][colsPerProc + 1] = eastExchange[i];
				myBoard[i][0] = westExchange[i];

				eastExchange[i] = myBoard[i][colsPerProc];
				westExchange[i] = myBoard[i][1];
			}
			MPI_Send(westExchange, rowsPerProc, MPI_CXX_BOOL, westNB, 3, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << westNB << std::endl;
#endif
			MPI_Send(eastExchange, rowsPerProc, MPI_CXX_BOOL, eastNB, 4, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << eastNB << std::endl;
#endif
		}

		if (myRow % 2 == 0)
		{
			for (int i = 0; i < colsPerProc; i++)
			{
				northExchange[i] = myBoard[1][i];
				southExchange[i] = myBoard[rowsPerProc][i];
			}
			MPI_Send(northExchange, colsPerProc, MPI_CXX_BOOL, northNB, 5, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << northNB << std::endl;
#endif
			MPI_Send(southExchange, colsPerProc, MPI_CXX_BOOL, southNB, 6, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << southNB << std::endl;
#endif
			MPI_Recv(southExchange, colsPerProc, MPI_CXX_BOOL, southNB, 7, MPI_COMM_WORLD, &Stat);
			MPI_Recv(northExchange, colsPerProc, MPI_CXX_BOOL, northNB, 8, MPI_COMM_WORLD, &Stat);
			for (int i = 0; i < colsPerProc; i++)
			{
				myBoard[rowsPerProc + 1][i] = southExchange[i];
				myBoard[0][i] = northExchange[i];
			}
		}
		else
		{
			MPI_Recv(southExchange, colsPerProc, MPI_CXX_BOOL, southNB, 5, MPI_COMM_WORLD, &Stat);
			MPI_Recv(northExchange, colsPerProc, MPI_CXX_BOOL, northNB, 6, MPI_COMM_WORLD, &Stat);
			for (int i = 0; i < colsPerProc; i++)
			{
				myBoard[rowsPerProc + 1][i] = southExchange[i];
				myBoard[0][i] = northExchange[i];

				northExchange[i] = myBoard[1][i];
				southExchange[i] = myBoard[rowsPerProc][i];
			}
			MPI_Send(northExchange, colsPerProc, MPI_CXX_BOOL, northNB, 7, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << northNB << std::endl;
#endif
			MPI_Send(southExchange, colsPerProc, MPI_CXX_BOOL, southNB, 8, MPI_COMM_WORLD);
#ifndef NDEBUG
			std::cout << "Processor " << myID << " send data to processor " << southNB << std::endl;
#endif
		}

		updateBoard(myBoard);
		if (i % iteration_gap == 0)
		{
			writeBoardToFile(myBoard, myFirstRow, myLastRow, myFirstCol, myLastCol, programName, i, myID);
		}
	}

	// Finalize the MPI environment. No more MPI calls can be made after this
	MPI_Finalize();
}