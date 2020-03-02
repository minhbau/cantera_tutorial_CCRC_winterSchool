#!###########################!###########################!##########################
# install cantera on Linux and MacOS (TESTED with lastest Anacona3)
#!###########################!###########################!##########################
# create cantera_py3 environment to install cantera and neccesary package
conda create -n cantera_py3  python=3
# activate the cantera_py3 environment
source activate cantera_py3
# Install cantera version 2.4 located in cantera chanel
conda install -c cantera cantera=2.4
# install some usefull packages
conda install -c conda-forge matplotlib  numpy ipython  scipy  spyder jupyter pandas


# -----------------------------------------------------------
# open your browser an go to
https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool
# download a zip file of  cantera_tutorial_CCRC_winterSchool and unzip it to your local machine
# then go the unzip directory, for example in my case
# -----------------------------------------------------------

# from "terminal", type
cd $HOME/Downloads/cantera_tutorial_CCRC_winterSchool
# and run the test
python example.py

# to logout the the cantera_py3 environment
conda deactivate
# re-login the cantera_py3 environment
source activate cantera_py3


#!###########################!###########################!##########################
# install cantera on Windows (TESTED with lastest Anacona3 conda -V 4.7.12)
#!###########################!###########################!##########################
# NOTE: Open the Anaconda Prompt terminal and type the following command in sequence
# all the commands and run the code are done in "Anaconda Prompt terminal" 
# create cantera_py3 environment to install cantera and neccesary package
conda create -n cantera_py3  python=3
# activate the cantera_py3 environment
conda activate cantera_py3
# Install cantera version 2.4 located in cantera chanel
conda install -c cantera cantera=2.4
# install some usefull packages
conda install -c conda-forge matplotlib  numpy ipython  scipy  spyder jupyter pandas

# -----------------------------------------------------------
# open your browser an go to
https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool
# download a zip file of  cantera_tutorial_CCRC_winterSchool and unzip it to your local machine
# then go the unzip directory, for example in my case
# -----------------------------------------------------------
# from the "Anaconda Prompt terminal", type
cd C:\Users\minhbau\Downloads\cantera_tutorial_CCRC_winterSchool
# and run the test
python example.py

# to logout the the cantera_py3 environment
conda deactivate
# re-login the cantera_py3 environment
conda activate cantera_py3
