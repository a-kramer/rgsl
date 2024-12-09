# Installation Instructions

Using the `remotes` package:

```R
remotes::install_github("a-kramer/rgsl")
```

Ensure that the [GNU Scientific
Library](https://www.gnu.org/software/gsl/doc/html/index.html) (gsl)
is installed in your system and
[pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/) or
[pkgconf](http://pkgconf.org/) can find it:

```sh
$ pkg-config --libs gsl
```

The above command should print something similar to:

```sh
-lgsl -lgslcblas -lm
```
