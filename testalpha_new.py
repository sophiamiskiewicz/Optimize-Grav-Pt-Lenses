from qlens_helper import *
import numpy as np
#from qlens import *

# This is the definition for the function for the GoldenRatioSearch
def GoldenRatioSearch(f, params_1, params_4, epsilon):
    # params_1 and params_4 are the outermost values for the search
    f_1 = f(params_1)
    f_4 = f(params_4)

    # params_3 and params_2 should have all the same parameters as params_1 and params_4 to start
    params_3 = params_1
    params_2 = params_1

    # define the "Golden Ratio"
    z = (1 + np.sqrt(5))/2
    dist = np.abs(params_4 - params_1)

    # now we'll use 1/z to find the interior points
    params_3 = params_1 + dist*(1/z)
    params_2 = params_4 - dist*(1/z)

    f_3 = f(params_3)
    f_2 = f(params_2)

    while params_4 - params_1 > epsilon:

        # if f(params_2) is less than f(params_3) then we will move in the direction of params_2
        if f_2 < f_3:
            params_4 = params_3
            params_3 = params_2
            dist = np.abs(params_4 - params_1)
            params_2 = params_4 - dist*(1/z)

            f_2 = f(params_2)
            f_3 = f(params_3)

        # if f(params_3) is less than f(params_2) then we will move in the direction of params_2
        if f_2 > f_3:
            params_1 = params_2
            params_2 = params_3
            dist = np.abs(params_4 - params_1)
            params_3 = params_1 + dist*(1/z)

            f_2 = f(params_2)
            f_3 = f(params_3)

        # if f_2 = f_3, then we will try to continue narrowing down the search and set the outermost points to be params_2 and params_3
        if f_2 == f_3:
            print("uh-oh! looks like a degeneracy!")
            params_1 = params_2
            params_4 = params_3

            dist = np.abs(params_4 - params_1)

            # now we'll use 1/z to find the interior points
            params_3 = params_1 + dist*(1/z)
            params_2 = params_4 - dist*(1/z)

            f_3 = f(params_3)
            f_2 = f(params_2)

    # the final solution will be somewhere in the middle between params_2 and params_3
    sol = (params_2 + params_3)/2

    return sol

q = QLens()

q.sci_notation = False
# This is the input file, which gives us the image locations of the gravitational lensing system that we are fitting
q.imgdata_read("alphafit.dat")
q.imgdata_display()

# Edit the numbers after the parameter names in order to vary your initial parameter guess
# How the code is currently set, your guesses for alpha and s will be ignored
q.add_lens(SPLE_Lens({"b": 2.5, "alpha": 1, "s": 0.0, "q": 0.5, "theta": 0, "xc": 0.0, "yc": 0.2}))
# This line tells us what parameters we will be varying
q.lens[0].setvary([1,0,0,1,1,1,1])

# Edit the numbers after the parameter names in order to vary initial parameter guess for the shear
q.add_lens(Shear({"shear": 0.02, "theta": -20, "xc": 0.0, "yc": 0.2}))
# This line anchors the shear to the center of the gravitational lens
q.lens[1].anchor_center(0)
q.lens[1].setvary([1,1])

# These are various settings included in QLens
q.central_image = False
q.imgplane_chisq = True
q.flux_chisq = True
q.chisqtol = 1e-6
# QLens will automatically guess the best point for locating the source
q.analytic_bestfit_src = True
#q.set_sourcepts_auto()
q.nrepeat = 2
print("Fit model:")
q.fitmodel()

print(q.fitparams())
# The following two lines must be run before evaluating the likelihood
q.setup_fitparams(True) # This sets up the fit parameters
initial_params = q.fitparams()  # This will return a python list giving the fit parameters
q.init_fitmodel() # This makes a copy of the model (lenses, sources, etc.) that can be varied during a fit (this is called the "fitmodel" object)
q.LogLike(q.fitparams())
fit_plotimg(q)
pause()  # must run the script with -i to run in interactive mode
# after pausing, press control-d to run the script

# this is where I define the function to run Powell's method
# used ChatGPT & Numerical Recipes for some various parts of this function
def Powells_method(f, params_initial, iterations):
    # we'll want to start by initializing a matrix that gives us the initial directions we will search in
    # the number of dimensions of the matrix will be dependent upon how many parameters we are optimizing
    # in this case we start with 7 parameters, and create a 7x7 identity matrix to begin the optimization
    dimension = len(params_initial)
    params = np.copy(params_initial)
    differences = []
    direction_matrix = np.eye(dimension)
    params = np.array(params)

    # these are the search ranges for the parameters that we are optimizing
    # ex: for parameter 0, we are searching in a range of +/- 2.5 from the initial parameter guess
    search_range = [2.5, 2.0, 115, 1.0, 1.0, 1.0, 90]

    # shows us the initial direction matrix
    print(direction_matrix)

    # this is where we begin the implementation of Powell's method
    for j in range(iterations):

        # the iteration number represents how many times we have done a complete optimization across all the directions
        print("this is the iteration: ", j)
        params_init = np.copy(params)

        for i in range(dimension):
            # we will iterate over i to search over each direction (which is represented by a single row of the direction_matrix)
            print("this is the index: ", i)
            params_old = np.copy(params)
            direction_vector = direction_matrix[i]

            # this is our line minimization function -- we minimize the parameter t, which will be added to the guess/updated guess in the appropriate direction we are optimizing in
            def line_minimization(t):
                return f(list(params + direction_vector*t))

            # for our first iteration, we will use the simple search range that is defined in line 123
            if j == 0:
                minimum_t = GoldenRatioSearch(line_minimization, -1*search_range[i], search_range[i], 1e-6)

            # for all the following iterations, we will alter the search range slightly based on the new direction we are traveling in (if we are traveling in a new direction)
            if j > 0:
                newsearch = search_range*direction_vector
                newrange = np.linalg.norm(newsearch)
                minimum_t = GoldenRatioSearch(line_minimization, -1*newrange, newrange, 1e-6)

            # the new parameters after finding the ideal t value
            params = params_old + direction_vector*minimum_t
            # we keep track of the change between each parameter in each direction
            differences.append(np.abs(params[i] - params_old[i]))

        # this is the sharpest decrease/increase between all the directions
        maximum = max(differences)

        f_E = f(list(2*params - params_init))
        f_0 = f(list(params_init))
        f_N = f(list(params))

        # this is a condition from Numerical Recipes for if we change directions in the next iteration over the variable j, we will cut out out the direction that had the sharpest change in parameter value(s)
        if f_E < f_0 or 2*(f_0 - 2*f_N + f_E)*((f_0-f_N)-maximum)**2 < ((f_0 - f_E)**2)*maximum:
            new_direction = params - params_init
            direction_of_change = differences.index(maximum)
            print("this is direction of change: ", direction_of_change)
            #print(new_direction)
            direction_matrix[direction_of_change] = new_direction
            print(direction_matrix)

        else:
            print("no direction change")

        differences = []

        print("new params: ", params)

        # optional statement if you want to end the loop when the parameters are barely shifting anymore
        #if np.linalg.norm(params - params_init) < 1e-6:
            #print('reached sufficiently small interval')
            #break

    return params

# you can alter the final value in line 191 to change the number of iterations
solution = Powells_method(q.LogLike, initial_params, 300)
txt = f"b: {solution[0]}, q: {solution[1]}, theta: {solution[2]}, xc: {solution[3]}, yc: {solution[4]}, shear: {solution[5]}, theta: {solution[6]}"
print(txt)

new_params = list(solution)

# QLens destroys the previous fitmodel object, so we have to re-initialize it
q.init_fitmodel()
# the new LogLikelihood value
print(q.LogLike(new_params))
# adopt the new parameters for the fit
q.adopt_model(new_params)
# show what the final optimized fit looks like
fit_plotimg(q)

plt.show()
