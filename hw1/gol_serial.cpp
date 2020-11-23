/**
 * Serial implementation of the game of life for the course Parallel Computing.
 * Emil Loevbak (emil.loevbak@cs.kuleuven.be)
 * First implementation: November 2019
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>

int const globalBufferLength = 50;

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
    std::ostream_iterator<bool> outputIterator(outputFile, "\t");
    for (size_t i = 0; i < board.size(); ++i)
    {
        copy(board[i].begin(), board[i].end(), outputIterator);
        outputFile << std::endl;
    }
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
    size_t rows, cols;
    int iteration_gap, iterations;
    try
    {
        rows = atoi(argv[1]);
        cols = atoi(argv[2]);
        iteration_gap = atoi(argv[3]);
        iterations = atoi(argv[4]);
    }
    catch (std::exception const &exc)
    {
        std::cout << "One or more program arguments are invalid!" << std::endl;
        return 1;
    }
    int processes = 1, processID = 0;
    size_t firstRow = 0, lastRow = rows - 1, firstCol = 0, lastCol = cols - 1;
    std::string programName = setUpProgram(rows, cols, iteration_gap, iterations, processes);

    //Build board
    std::vector<std::vector<bool>> board((lastRow - firstRow + 1), std::vector<bool>(lastCol - firstCol + 1));
    initializeBoard(board);

    //Do iteration
    writeBoardToFile(board, firstRow, lastRow, firstCol, lastCol, programName, 0, processID);
    for (int i = 1; i <= iterations; ++i)
    {
        updateBoard(board);
        if (i % iteration_gap == 0)
        {
            writeBoardToFile(board, firstRow, lastRow, firstCol, lastCol, programName, i, processID);
        }
    }
}