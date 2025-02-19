# ParticleSimulation
Implementation of the [particle-mesh method](https://www.cs.cmu.edu/afs/cs/academic/class/15850c-s96/www/nbody.html#pm).

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

## Output viewing
### 3D
```bash
python script/animate.py <output file path>
```
For example (MSVC in Debug mode):
```bash
python script/animate.py build/source/Debug/output.txt
```
### 2D
```bash
python script/animate2d.py <output file path>
```
For example (MSVC in Debug mode):
```bash
python script/animate2d.py build/source/Debug/output.txt
```

### Final positions
```bash
python script/show_last.py <output file path>
```
For example (MSVC in Debug mode):
```bash
python script/show_last.py build/source/Debug/output.txt
```

## Total energy & momentum plots
```bash
python script/conserved_quantities.py <energy file path> <momentum file path>
```
For example (MSVC in Debug mode):
```bash
python script/conserved_quantities.py build/source/Debug/energy.txt build/source/Debug/momentum.txt
```
## Results
### Spiral galaxy simulation
![](results/img/spiral_galaxy.png)