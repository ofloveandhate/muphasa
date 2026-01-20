# Muphasa

Muphasa is a software to compute presentations of multi-parameter persistent homology. It was devoloped by Matías R. Bender, Oliver Gäfvert, and Michael Lesnick.

## Instalation

To install Muphasa you need CMake. There are two options:

- If you want to install Muphasa, together with its Python package and interface, you should be able to use the setuptools build system from the root folder.  Simply use `pip install .` from root in this repo.

- If you want to compile only its console version, you can use the following lines, 
   ```
   mkdir build; cd build
   cmake ../src/cpp
   cmake --build .
   ```

## Examples

- If you installed the Python interface, you can try the notebook file [/python/notebooks/examples.ipynb](https://github.com/olivergafvert/muphasa/blob/main/python/notebooks/examples.ipynb) .

- If you want to test the console version, from the build folder, you can run the following lines
  ```
  ./mph ../examples/example1.txt
  ./mph --firep ../examples/example2.firep
  ```
