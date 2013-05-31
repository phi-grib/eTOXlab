PREDICT

Usage is (predict -h shows usage message)

	predict [-f filename.sdf][-v 1]

Can be run using the interpreter

	python ./predict.py -v -1 -f mini.sdf

or just calling the code (#!/usr/bin/env python)

	./predict.py -v -1 -f mini.sdf

BUILD

Usage is (build -h shows usage message)

	build.py [-f filename.sdf][-v 1]

Can be run using the interpreter

	python ./build.py -v -1 -f mini.sdf

or just calling the code (#!/usr/bin/env python)

	./build.py -v -1 -f mini.sdf

For building a customized version edit the build function of model in the imodel.py if the previous version.  


DOCUMENTATION

We can use the program "pydoc" included in the Python programming utilities. For example we can set up a server:

	pydoc -p 8080

And then, all the documentation is accessible using a web browser pointing to 

	localhost:8080

Alternativelly, it can be dumped to HTML using

	pydoc -w predict

... will produce predict.html
