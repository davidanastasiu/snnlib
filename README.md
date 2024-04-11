### Local install
```bash
conda create --name snnlib python=3.10
conda activate snnlib
conda install -c conda-forge gcc gxx
conda install -c anaconda cmake
python -m pip install jupyter scikit-learn matplotlib  # extras for demos and such
cd ssnlib
cmake -B build -DCMAKE_BUILD_TYPE=RELEASE
make -C build
python -m pip install .
```


### Compile tests and demos
```bash
python -m pip install "pybind11[global]"
cmake -B build -DCMAKE_BUILD_TYPE=RELEASE
make -C build
```
To see more details during make:
```bash
make VERBOSE=1 -C build
```

### Python install
```bash
python -m pip install .
```

To isntall in development mode:
```python
python -m pip install -e .
```

Check out `test.ipynb` for some examples of running KNNG in Python.

### Run a subset of the tests
```bash
build/tests/snnlib_test --gtest_filter=sort.*
```

### Work in progress
This project is a work in progress. More work is needed to finalize nearest neighbor search methods and other driver utilities.
