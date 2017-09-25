Algorithm documentation
=======================
While the tutorial in the main Readme should help you get up and running the documents 
in this folder are meant to give some insight into the internal structure and 
working of the programs.
 

- [Read generation POS and CIGAR algorithm](readgeneration.md)
- [Alignment scoring](alignmentscoring.md)
- [Read model format](readmodelformat.md)


Docker image creation
=====================

### Build
```
cd mitty
docker build -f containerization/dockerfile -t images.sbgenomics.com/kghosesbg/mitty3:latest .
```

### Test
Run the docker image with the demo scripts folder mounted

```
docker run -ti --rm -v /Users/kghose/Code/mitty-demo-data:/data images.sbgenomics.com/kghosesbg/mitty3:latest
```

And then run individual test scripts to check if everything works on the Linux container

### Push
```
docker push images.sbgenomics.com/kghosesbg/mitty3:latest
```


Pushing to PyPi
===============

(via https://packaging.python.org/tutorials/distributing-packages)

```
python setup.py sdist
python setup.py bdist_wheel
# Fill out ~/.pypirc if you haven't already
twine upload dist/*
```
