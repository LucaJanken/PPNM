CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra
LDFLAGS=-lstdc++ -lm

# Object files
OBJS=main.o quadrature.o

# Executable name
EXEC=main

# Default target
all: $(EXEC) Out.txt erf_plot diff_plot

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $(EXEC) $(LDFLAGS)

# Error function tabulated file in ../../exercises/plots/erf_tab.txt
erf_data = ../../exercises/plots/erf_tab.txt

# Error function approximate values from ../../exercises/plots/erf_data.txt
erf_approx = ../../exercises/plots/erf_data.txt

# Object file dependencies
main.o: main.cc quadrature.h
	$(CXX) $(CXXFLAGS) -c main.cc

quadrature.o: quadrature.cc quadrature.h
	$(CXX) $(CXXFLAGS) -c quadrature.cc

Out.txt: $(EXEC)
	./$(EXEC) > Out.txt
	python main.py >> Out.txt  # Append Python script output directly after $(EXEC) execution

# Plot the error function alongside the tabulated and approximate values
erf_plot: erf_plot.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'erf_plot.svg'; \
	set title 'Error function'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set xrange [-3:3]; \
	set grid; \
	plot 'erf_plot.csv' using 1:2 with lines title 'erf(x)', \
	'"$(erf_data)"' using 1:2 with points pointtype 1 linecolor 'green' title 'erf_{data}(x)', \
	'"$(erf_approx)"' using 1:2 with points pointtype 6 title 'erf_{approx}(x)', \
	'"$(erf_data)"' using (\$$1 * -1):(\$$2 * -1) with points pointtype 1 linecolor 'green' notitle;"
	
diff_plot: diff.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'diff_plot.svg'; \
	set title 'Difference between erf implementations'; \
	set xlabel 'x'; \
	set ylabel 'Difference'; \
	set xrange [0:3]; \
	set grid; \
	plot 'diff.csv' using 1:2 with points pointtype 1 linecolor 'red' title 'Difference (My erf - Tabulated)', \
	'diff.csv' using 1:3 with points pointtype 6 linecolor 'blue' title 'Difference (Approx - Tabulated)';"


# Clean the build
clean:
	rm -f $(EXEC) $(OBJS) Out.txt erf_plot.csv diff.csv erf_plot.svg diff_plot.svg

