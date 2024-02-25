CXX = g++
CXXFLAGS = -std=c++11 -pthread
TARGET = main
SRC = main.cc
OBJ = $(SRC:.cc=.o)
N = 1e8
TIME = /usr/bin/time -p

# Default build: builds the serial version
all: $(TARGET) time_test time_test_openmp time_test_openmp_tls

# Builds the serial version
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiles the object file
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Time tests for the serial version
time_test:
	@echo "Running with 1 thread" > Out.txt
	@$(TIME) ./$(TARGET) -threads:1 -terms:$(N) >> Out.txt 2>&1
	@echo "Running with 2 threads" >> Out.txt
	@$(TIME) ./$(TARGET) -threads:2 -terms:$(N) >> Out.txt 2>&1
	@echo "Running with 3 threads" >> Out.txt
	@$(TIME) ./$(TARGET) -threads:3 -terms:$(N) >> Out.txt 2>&1
	@echo "Running with 4 threads" >> Out.txt
	@$(TIME) ./$(TARGET) -threads:4 -terms:$(N) >> Out.txt 2>&1

# Time tests for the OpenMP version
time_test_openmp:
	@echo "Compiling with OpenMP..."
	@$(CXX) $(CXXFLAGS) -fopenmp -DUSE_OPENMP $(SRC) -o $(TARGET)_openmp
	@echo "When multiple threads simultaneously attempt to update the shared totalSum variable without proper synchronization, some updates are lost." > Out_OpenMP.txt
	@echo "Running OpenMP version" >> Out_OpenMP.txt
	@for thread_count in 1 2 3 4 ; do \
		echo "Running with $$thread_count threads" >> Out_OpenMP.txt; \
		$(TIME) ./$(TARGET)_openmp -threads:$$thread_count -terms:$(N) >> Out_OpenMP.txt 2>&1; \
	done

# Time tests for the OpenMP version with Thread-local Storage
time_test_openmp_tls:
	@echo "Compiling with OpenMP and Thread-local Storage..."
	@$(CXX) $(CXXFLAGS) -fopenmp -DUSE_OPENMP -DUSE_TLS $(SRC) -o $(TARGET)_openmp_tls
	@echo "OpenMP's reduction operation automatically handles the "thread-local storage" aspect for the variable being reduced across threads, so no updates are lost." > Out_OpenMP_TLS.txt
	@echo "Running OpenMP with Thread-local Storage version" >> Out_OpenMP_TLS.txt
	@for thread_count in 1 2 3 4 ; do \
		echo "Running with $$thread_count threads" >> Out_OpenMP_TLS.txt; \
		$(TIME) ./$(TARGET)_openmp_tls -threads:$$thread_count -terms:$(N) >> Out_OpenMP_TLS.txt 2>&1; \
	done


clean:
	rm -f $(OBJ) $(TARGET) $(TARGET)_openmp $(TARGET)_openmp_tls Out.txt Out_OpenMP.txt Out_OpenMP_TLS.txt

.PHONY: all clean time_test time_test_openmp time_test_openmp_tls
