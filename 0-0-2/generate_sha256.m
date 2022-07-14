#generating sh256 for .yaml file:
#1)release version on Github
#2)get the address of .tar.gz source (e.g. v0.0.1.tar.gz):
url="https://github.com/aumpierre-unb/
Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/v0.0.1.tar.gz"
#3)set destination .tar.gz file (in the scurrent directory):
file="teste.tar.gz"
#4)generate file from url:
urlwrite(url,file)
#5)generate sha256 for file:
hash("sha256",fileread(file))
#copy and paste ans to .yaml file
