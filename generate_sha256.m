
# how to get the hash string sha256 for .yaml file:

# 1) Release version on Github under a tag, e.g.
tag="v0.2.0-alpha"

# 2) Execute this line to check on consistency of the package:

pkg install("https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/v0.1.1.tar.gz")

# 3) Get the address of .tar.gz:
url=["https://github.com/aumpierre-unb/" ...
     "Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/" ...
     tag ...
     ".tar.gz"]

# 4) Set file:
file=[tag ".tar.gz"]

# 5) Generate file from url:
urlwrite(url,file)

# 6) Generate hash string sha256 for .yaml:
hash("sha256",fileread(file))

# 7) copy and paste hash string under quotes to .yaml file



