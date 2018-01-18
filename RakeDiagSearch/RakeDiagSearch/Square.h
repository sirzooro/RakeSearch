// Diagonal Latin square

# if !defined Square_h
# define Square_h

# include <iostream>

using namespace std;

class Square
{
public:
	static const int Rank = 9;			// Rank of the square
	static const int Empty = -1;			// Empty value
	static const char HeadToken = '{';		// Start token for the data about the square in stream
	static const char TailToken = '}';		// End token for the data about the square in stream

	static int OrthoDegree(const Square& a, const Square& b);	// Orthogonality degree for squares a and b

	Square();					// Default constructor
	Square(int source[Rank][Rank]);			// Matrix constructor
	Square(const Square& source);				// Copy constructor

	int operator == (const Square& value) const;
		// Overloading the comparison operator: the matrix components are being compared
	Square& operator = (const Square& value);
		// Overloading the assignment operator
	friend std::ostream& operator << (std::ostream& os, const Square& value);
		// Overloading the square data output operator
	friend std::istream& operator >> (std::istream& is, Square& value);
		// Overloading the square data input operator

	int IsDiagonal() const;	// Checking the square for being diagonal Latin square 
	int IsLatin() const;		// Checking the square for being Latin square
	void Initialize(const int source[Rank][Rank]);	// Initialization of the square components
	void Reset();					// Resetting all values of variables
	void Read(std::istream& is);			// Reading the square from stream
	void Write(std::ostream& os) const;			// Writing the square into stream

	int Matrix[Rank][Rank];				// Matrix of the square

protected:
private:
};

# endif
