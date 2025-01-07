# ParticleSimulation

## Build
```bash
mkdir build; cd build
cmake ..
cmake --build .
```

## Test
```bash
ctest
```

## Animation
```bash
python script/animate.py <output file path>
```
For example (MSVC in Debug mode):
```bash
python script/animate.py build/source/Debug/output.txt
```

## Total energy & momentum plots
```bash
python script/conserved_quantities.py <energy file path> <momentum file path>
```
For example (MSVC in Debug mode):
```bash
python script/conserved_quantities.py build/source/Debug/energy.txt build/source/Debug/momentum.txt
```