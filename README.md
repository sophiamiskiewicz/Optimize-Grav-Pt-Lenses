# Optimize-Grav-Pt-Lenses
This is my final project for Fall 2024 Computational Methods in Astrophysics. I created an optimization script for point lenses, to be used with the lensing simulation software QLens.

In order to execute this code, you will first need to git clone or download as a ZIP file the qlens-beta development branch. This is a gravitational lensing software developed by Quinn Minor, and I have forked that repository and put it on my GitHub here: https://github.com/sophiamiskiewicz/qlens-beta/tree/development. First, let's go over the steps to install the python wrapper for QLens:
1. Either download the ZIP file or clone the repository at the above link.
2. Navigate to the QLens development directory in your terminal and then go into the python directory with the command cd python/.
3. Create a directory called build with the command mkdir build.
4. Stay in the python directory and execute the cmake .. command. There is a fair chance this won't immediately work unless you specify where python is located on your computer. To do this, execute the command 'which python' in your terminal and copy that path. Then execute the command as cmake .. -DPYTHON_EXECUTABLE=insertyourcopiedpathhere. If that doesn't work, you may need to find the path for your python include files and then execute the command as cmake .. -DPYTHON_EXECUTABLE=insertyourcopiedpathhere -DPYTHON_INCLUDE_DIR=insertpathtopythonincludefileshere.
5. Afterwards, you'll want to execute the command make deps, and then the command make python.

From there, QLens should be compiled for you! If something goes awry, please feel free to contact me at smiskiewicz@gradcenter.cuny.edu. Make sure you have the g++ compiler installed for this.

Once QLens is compiled, then download or git clone my script: testalpha_new.py and place it in the python directory. To run the script, go to the python directory and use the command python -i testalpha_new.py. To modify the parameters being guessed in the script, you can add adjust the numbers in lines 77 and 82. To modify the search range, you can change line 123. And finally, you can change the number of iterations of Powell's method in line 191, where it originally says "300". When you initially run the script, it will show you how close your guesses are to the lens images in the figure with blue/green dots which represent our model vs. generated data for the lens. When the script pauses, press Ctrl d to continue the code, and you will see the results of the final optimization.  
