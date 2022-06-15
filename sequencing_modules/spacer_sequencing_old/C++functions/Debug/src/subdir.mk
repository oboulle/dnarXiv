################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/clustering.cpp \
../src/clustering_keys.cpp \
../src/old_clustering.cpp 

OBJS += \
./src/clustering.o \
./src/clustering_keys.o \
./src/old_clustering.o 

CPP_DEPS += \
./src/clustering.d \
./src/clustering_keys.d \
./src/old_clustering.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


