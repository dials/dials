To set up a DIALS development environment you start off with an empty directory.

On Linux and MacOS you can then run
```bash
$ wget https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py
$ python bootstrap.py
```
which will install a complete development environment into the directory `conda_base` and
download, configure, and compile [cctbx](https://github.com/cctbx/cctbx_project), DIALS,
and [xia2](https://github.com/xia2/xia2). You will find all the source code repositories
in the `modules` directory.

Once built, you need to load the DIALS development environment by running
```bash
$ . dials
```

For more information please go to https://dials.github.io/
