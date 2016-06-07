StoneX
===========

Description
-----------
*StoneX* implements a temperature dependent coherent magnetization reversal model for size-distributed assemblies of ferromagnetic nanoparticles and ferromagnetic-antiferromagnetic core-shell nanoparticles.
The model is based on coherent magnetization reversal, with thermal fluctuation using the Néel-Arrhenius law.
Non-thermal model are also proposed, like the Stoner-Wohlfarth and Meiklejohn-Bean models.
The program allow to produce angular hysteresis behavior, as well as temperature dependent hysteresis evolution.


Citations and licence
---------------------
If you use StoneX in any work or publication, we kindly ask you to cite our article:

“[Thermal simulation of magnetization reversals for size-distributed assemblies of core-shell exchange biased nanoparticles](https://arxiv.org/abs/1602.07877), J. Richy *et al.*, ArXiv e-prints, **1602**.07877, 2016”


StoneX is open-source software. You are free to modify and distribute the source code under the GPLv3 licence.


Installation / Dependencies
-------------
This program uses Python version 3, and can therefore be used in any platform. See the specific installation instructions
The following packages are needed:
- sys, shutil, time
- numpy
- matplotlib
- scipy

For Mac Os X, using the package manager *macport* and the *ipython* shell:
```
sudo port selfupdate
sudo port install py35-ipython py35-scipy py35-matplotlib
```

For GNU/Linux, using the package manager *apt* and the *ipython* shell:
```
sudo apt update
sudo apt install ipython3 python3-scipy python3-matplotlib
```


How to use
------------
The *StoneX* program is composed of a module (contained in the StoneX folder), and a script for controlling the parameters and launching the simulation.

The *examples* folder contain standard script examples for running the simulation. To indicate the module path, the easiest way is to  modify the variable *module_path* in the script file.
For example,
```
module_path = '/home/Me/Python_Modules/StoneX_project/'
```

The script can also be placed in the same directory than the StoneX folder.
```
├── Project
|   ├── StoneX/
|   └── script.py
```

Contributing
------------
Contributions are gratefully accepted. To contribute code, fork our repo on github and send a pull request.
