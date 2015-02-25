#!/bin/bash

# Install glm
curl -L http://sourceforge.net/projects/ogl-math/files/glm-0.9.6.3/glm-0.9.6.3.zip > ./glm.zip
unzip glm.zip

# Clean up
rm glm.zip
