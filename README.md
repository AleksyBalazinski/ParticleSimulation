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

## Total energy & momentum plots
```bash
python script/conserved_quantities.py <energy file path> <momentum file path>
```