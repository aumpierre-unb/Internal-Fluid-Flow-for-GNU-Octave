# Generate hash string sha256 for .yaml file:

## 1) Release version on Github under a tag, e.g.
tag="v0.2.0"

## 2) Get the address of .tar.gz:
url=["https://github.com/aumpierre-unb/" ...
     "Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/" ...
     tag ".tar.gz"]

## 3) Set file name:
file=[tempdir tag ".tar.gz"]

## 4) Generate file from url:
urlwrite(url,file)

## 5) Generate hash string sha256 for .yaml:
hash("sha256",fileread(file))

## 6) Copy and paste hash string under quotes to .yaml file
