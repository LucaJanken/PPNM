# Compiler settings
CXX = g++ 
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Source files
SRCS = main.cc ../matrix/vector.cc spline.cc

# Datafile
LINDATAFILE = lspline.csv
QDATAFILE = qspline.csv
CDATAFILE = cspline.csv

# Object files
OBJS = $(SRCS:.cc=.o)

# Target executable
TARGET = main

# Default targets
all: lspline.svg qspline.svg cspline.svg

# Rule to generate the output file
lspline.csv: $(TARGET)
	./$(TARGET) > $@

# Rule to create the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDLIBS)

# Rule to create the object files
%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Rule to generate the plot
lspline.svg: $(LINDATAFILE)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'lspline.svg'; \
	set title 'Linear Spline'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set yrange [-2:2]; \
	plot \
	'$(LINDATAFILE)' using 1:2 with linespoints title 'f(x) = cos(x)', \
	'$(LINDATAFILE)' using 1:3 with linespoints title 'F(x) = sin(x)';"

# Rule to generate the plot
qspline.svg: $(QDATAFILE)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'qspline.svg'; \
	set title 'Quadratic Spline'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set yrange [-2:2]; \
	plot \
	'$(QDATAFILE)' using 1:2 with linespoints title 'f(x) = sin(x)', \
	'$(QDATAFILE)' using 1:4 with linespoints title 'df(x)/dx = cos(x)', \
	'$(QDATAFILE)' using 1:3 with linespoints title 'F(x) = -cos(x)';"

# Rule to generate the plot
cspline.svg: $(CDATAFILE)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'cspline.svg'; \
	set title 'Cubic Spline'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set yrange [-2:2]; \
	plot \
	'$(CDATAFILE)' using 1:2 with linespoints title 'f(x) = sin(x)', \
	'$(CDATAFILE)' using 1:4 with linespoints title 'df(x)/dx = cos(x)', \
	'$(CDATAFILE)' using 1:3 with linespoints title 'F(x) = -cos(x)';"

# Clean
clean:
	rm -f $(TARGET) $(OBJS) lspline.svg lspline.csv qspline.svg qspline.csv cspline.svg cspline.csv

# Phony
.PHONY: all clean
