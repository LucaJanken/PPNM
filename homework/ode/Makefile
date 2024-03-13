all: main plot1 plot2 plot4 plot5 plot6 plot7

main: main.cc 
	g++ main.cc rkode.cc ../matrix/vector.cc ../splines/spline.cc -std=c++11 -lstdc++ -lm -o main

shm.csv: main
	./main

shm_friction.csv: main
	./main

ncm.csv: main
	./main

nem.csv: main
	./main

plot1: shm.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'shm.svg'; \
	set title 'Harmonic oscillator'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot 'shm.csv' using 1:2 with lines title 'u(x) = cos(x)', \
	'shm.csv' using 1:3 with lines title 'v(x) = -sin(x)'" \

plot2: shm_friction.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'shm_friction.svg'; \
	set title 'Harmonic oscillator with Friction'; \
	set xlabel 't'; \
	set ylabel 'y'; \
	set grid; \
	plot 'shm_friction.csv' using 1:2 with lines title 'theta(t)', \
	'shm_friction.csv' using 1:3 with lines title 'omega(t)'" \

plot4: ncm.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'ncm.svg'; \
	set title 'Newtonian Circular Motion'; \
	set size ratio -1;\
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot 'ncm.csv' using (1/column(2))*cos(column(1)):(1/column(2))*sin(column(1)) with lines notitle" \

plot5: nem.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'nem.svg'; \
	set title 'Newtonian Elliptical Motion'; \
	set size ratio -1;\
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot 'nem.csv' using (1/column(2))*cos(column(1)):(1/column(2))*sin(column(1)) with lines notitle;"

plot6: rp.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'rp.svg'; \
	set title 'Relativistic Precession'; \
	set size ratio -1;\
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot 'rp.csv' using (1/column(2))*cos(column(1)):(1/column(2))*sin(column(1)) with lines notitle;"

plot7: figure8.csv
	gnuplot -e "\
		set terminal svg size 800,600 background 'white'; \
		set output 'figure8.svg'; \
		set title 'Figure-8 Orbit'; \
		set xlabel 'x'; \
		set ylabel 'y'; \
		set size ratio -1;\
		plot '$<' using 1:2 with lines title 'Planet A', \
		            '$<' using 3:4 with lines title 'Planet B', \
		            '$<' using 5:6 with lines title 'Planet C';"


clean:
	rm -f main *.svg *.csv

.PHONY: clean all plot1 plot2 plot4 plot5 plot6 plot7