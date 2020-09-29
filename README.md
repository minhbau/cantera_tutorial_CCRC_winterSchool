
# Cantera tutorial CCRC winter school 2019

minhbau.luong@kaust.edu.sa

Refs:
```
https://cantera.org/
https://cantera.org/tutorials/python-tutorial.html
https://cantera.org/examples/jupyter/index.html
```


# Install cantera on Linux and MacOS (Tested with lastest Anacona3)

Create cantera_py3 environment to install cantera and neccesary packages
```
conda create -n cantera_py3  python=3
```

Activate the cantera_py3 environment
```
source activate cantera_py3
```
Install cantera version 2.4 

```
conda install -c cantera cantera=2.4
```

Install some usefull packages
```
conda install -c conda-forge matplotlib  numpy ipython  scipy  spyder jupyter pandas
```



Open your browser an go to
```
https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool
```
to download a zip file of cantera_tutorial_CCRC_winterSchool and unzip it to your local machine
then go the unzip directory, for example 

```
cd $HOME/Downloads/cantera_tutorial_CCRC_winterSchool

```
or from "terminal", type

```
cd $HOME/Downloads/
git clone https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool.git

```

and run a test
```
cd $HOME/Downloads/cantera_tutorial_CCRC_winterSchool
python example.py
```

To logout the the cantera_py3 environment
``` 
conda deactivate
```
Re-login the cantera_py3 environment
```
source activate cantera_py3
```


# Install cantera on Windows (Tested with the lastest Anacona3 conda -v 4.7.12)


NOTE: Open the Anaconda Prompt terminal and type the following command in sequence

All the commands and run the code are done in "Anaconda Prompt terminal" 

Create cantera_py3 environment to install cantera and neccesary package
```
conda create -n cantera_py3  python=3
```
Activate the cantera_py3 environment
```
conda activate cantera_py3
```

Install cantera version 2.4 located in cantera chanel
```
conda install -c cantera cantera=2.4
```

Install some usefull packages
```
conda install -c conda-forge matplotlib  numpy ipython  scipy  spyder jupyter pandas
```

Open your browser an go to

```
https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool
```

Download a zip file of cantera_tutorial_CCRC_winterSchool and unzip it to your local machine
then go the unzip directory, for example in my case
from the "Anaconda Prompt terminal", type
```
cd C:\Users\minhbau\Downloads\cantera_tutorial_CCRC_winterSchool
```
and run the test
```
python example.py
```
To logout the the cantera_py3 environment
```
conda deactivate
```
Re-login the cantera_py3 environment
```
conda activate cantera_py3
```
