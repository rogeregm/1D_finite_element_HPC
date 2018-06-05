# HPC-Coursework
# Roger Gonzalez
# 26/03/17

CC=g++
PARA=mpicxx
CXXFLAGS=-O0
LDLIBS=-lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas
HDRS=functions
TARGET=compile
OBJS=Q1_main.cpp Common_func.cpp Q1_func.cpp
		# length - Elements - Area - Inertia - Youngs Modulus - time - timesteps - density
PARAMETERS=10 24 0.012 0.0000144 210000000000.0 1.0 10000.0 7850.0

all: $(TARGET)
default: $(TARGET)

$(TARGET): $(OBJS)
	$(PARA)  $(OBJS) $(LDLIBS) -o $(TARGET)

task1: $(TARGET)
	mpiexec -np 1 ./$(TARGET) $(PARAMETERS) 0 0
	# python plotting.py

task2: $(TARGET)
	mpiexec -np 1 ./$(TARGET) $(PARAMETERS) 1 0
	# python plotting2.py

task3: $(TARGET)
	mpiexec -np 1 ./$(TARGET) $(PARAMETERS) 1 1
	# python plotting3.py

task4: $(TARGET)
	mpiexec -np 2 ./$(TARGET) $(PARAMETERS) 2 0
	# python plotting4.py

task5: $(TARGET)
	mpiexec -np 4 ./$(TARGET) $(PARAMETERS) 2 1
	# python Plotting5.py

# Specify clean is not a real target
.PHONY: clean	
	
# Clean up 
clean:		
	rm -f *.o	
