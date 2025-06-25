#pragma once

void galaxySimulationPM(const char* outputDir);

void galaxySimulationPMTiming(const char* outputDir);

void galaxyCollisionPM(const char* outputDir);

void galaxySimulationP3M(const char* outputDir);

void galaxySimulationP3MTiming(const char* outputDir);

void galaxyCollisionP3M(const char* outputDir);

void galaxySimulationBH(const char* outputDir);

void galaxySimulationBHTiming(const char* outputDir);

void anisotropy(const char* outputDir);

void finiteDiffRelError(const char* outputDir);

void assignmentSchemes(const char* outputDir);

void poorLaplace(const char* outputDir);

void pmAccuracy(const char* outputDir, bool ppSoftening);

void pmOptimal(const char* outputDir);

void pmOptimalVaryingDiameter(const char* outputDir);

void pmOptimalAssignments(const char* outputDir);

#ifndef CUDA
void p3mAccuracyAssignments(const char* outputDir);
void p3mAccuracyShapes(const char* outputDir);
#endif

void bhAccuracy(const char* outputDir);

void plummerBH(const char* outputDir);

void plummerPM(const char* outputDir);

#ifndef CUDA
void plummerP3M(const char* outputDir);
#endif

void galaxyCollisionBH(const char* outputDir);

void galaxySimulationPPTiming(const char* outputDir);