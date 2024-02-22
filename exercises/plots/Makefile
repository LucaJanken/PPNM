# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -O2

# Target executable
EXECUTABLE = functions.exe

# Plot files (SVG)
PLOT_ERF = erf_plot.svg
PLOT_GAMMA = gamma_plot.svg
PLOT_LNGAMMA = lngamma_plot.svg

# Tabulated files
TABULATED_ERF = erf_tab.txt
TABULATED_GAMMA = gamma_tab.txt
TABULATED_LNGAMMA = lngamma_tab.txt

# Data files
DATAFILE_ERF = erf_data.txt
DATAFILE_GAMMA = gamma_data.txt
DATAFILE_LNGAMMA = lngamma_data.txt

# Default target
all: $(EXECUTABLE) $(PLOT_ERF) $(PLOT_GAMMA) $(PLOT_LNGAMMA)

# Rule to create the target executable
$(EXECUTABLE): main.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to generate the output file
$(DATAFILE_ERF) $(TABULATED_GAMMA) $(TABULATED_LNGAMMA) $(DATAFILE_GAMMA) $(DATAFILE_LNGAMMA): $(EXECUTABLE)
	./$(EXECUTABLE)

# Rule to plot the error function
$(PLOT_ERF): $(DATAFILE_ERF) $(TABULATED_ERF) Makefile
	echo '\
		set terminal svg background "white"; \
		set output "$@"; \
		set key left; \
		set xlabel "x"; \
		set ylabel "erf(x)"; \
		set tics out; \
		set grid; \
		set title "Error function"; \
		plot [-3:3] \
			 "$(DATAFILE_ERF)" using 1:2 with lines title "erf(x)", \
			 "$(TABULATED_ERF)" using 1:2 with points pointtype 6 title "Tabulated Values", \
			 "$(TABULATED_ERF)" using ($$1 * -1):($$2 * -1) with points pointtype 6 title "Reflected Values"; \
	' | gnuplot

# Rule to plot the gamma function
$(PLOT_GAMMA): $(TABULATED_GAMMA) $(DATAFILE_GAMMA) Makefile
	echo '\
        set terminal svg background "white"; \
        set output "$@"; \
        set key left; \
        set xlabel "x"; \
        set ylabel "Gamma(x)"; \
		set xrange [-2:5]; \
		set yrange [-50:50]; \
        set tics out; \
        set grid; \
        set title "Gamma function"; \
        plot \
                 "$(DATAFILE_GAMMA)" using 1:2 with lines title "Continuous Gamma(x)", \
                 "$(TABULATED_GAMMA)" using 1:2 with points pointtype 6 title "Tabulated Gamma(x)"; \
' | gnuplot

# Rule to plot the natural logarithm of the gamma function
$(PLOT_LNGAMMA): $(TABULATED_LNGAMMA) $(DATAFILE_LNGAMMA) Makefile
	echo '\
        set terminal svg background "white"; \
        set output "$@"; \
        set key left; \
        set xlabel "x"; \
        set ylabel "ln(Gamma(x))"; \
		set xrange [0.1:5]; \
		set yrange [-50:50]; \
        set tics out; \
        set grid; \
        set title "Natural logarithm of the Gamma function"; \
        plot \
                 "$(DATAFILE_LNGAMMA)" using 1:2 with lines title "Continuous ln(Gamma(x))", \
                 "$(TABULATED_LNGAMMA)" using 1:2 with points pointtype 6 title "Tabulated ln(Gamma(x))"; \
' | gnuplot

clean:
	rm -f $(EXECUTABLE) $(PLOT_ERF) $(PLOT_GAMMA) $(PLOT_LNGAMMA) $(DATAFILE_ERF) $(TABULATED_GAMMA) $(TABULATED_LNGAMMA) $(DATAFILE_GAMMA) $(DATAFILE_LNGAMMA)

.PHONY: all clean
