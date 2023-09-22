
def make_recipe_size_distribution(stru_table,weights, dat_path):
    """
    Creates and returns a Fit Recipe object for a size distribution

    Parameters
    ----------
    stru_table : array of string, each string being the full path to the XYZ files to load.
    weigths : array of weights for each structure xyz file
    (stru_table and weights must have the same length)
    dat_path :  string, The full path to the PDF data to be fit.

    Returns
    ----------
    fitrecipe : The initialized Fit Recipe object using the datname and structure
                provided.
    """
    # 5: build array of diffpy Structure objects
    stru_array=[]
    
    for structure in stru_table:
        
        stru=Structure(filename = DPATH / structure)
        stru_array.append(stru)
    
    # 6: Create a Profile object for the experimental dataset and
    # tell this profile the range and mesh of points in r-space.
    profile = Profile()
    parser = PDFParser()
    parser.parseFile(DPATH / dat_path)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)    
    
    # 7: Create a Debye PDF Generator object for each discrete structure model.
    # build an array of generators as done for structure objects
    if len(stru_table)==len(weights):
        generator_cluster_array=[]
        contribution = FitContribution("cluster")
        # initialize index for iterative naming of variables inside the loop
        index=0
        for stru in stru_array:
            generator_cluster = DebyePDFGenerator("G%d"%index)
            generator_cluster.setStructure(stru, periodic=False)
            contribution.addProfileGenerator(generator_cluster) 
            generator_cluster.qdamp.value = QDAMP_I
            generator_cluster.qbroad.value = QBROAD_I
            generator_cluster.setQmax(QMAX)
            generator_cluster.setQmin(QMIN) 
            if RUN_PARALLEL:
                try:
                    import psutil
                    import multiprocessing
                    from multiprocessing import Pool
                except ImportError:
                    print("\nYou don't appear to have the necessary packages for parallelization")
                syst_cores = multiprocessing.cpu_count()
                cpu_percent = psutil.cpu_percent()
                avail_cores = np.floor((100 - cpu_percent) / (100.0 / syst_cores))
                ncpu = int(np.max([1, avail_cores]))
                pool = Pool(processes=ncpu)
                generator_cluster.parallel(ncpu=ncpu, mapfunc=pool.map)   
            # store the generator in an array for further use (variable definitions)
            generator_cluster_array.append(generator_cluster)
            index+=1            
        
        # 8: Create a contribution and Set an equation, based on your PDF generators.
        # the equation is: scale_factor*(sum(weigths[i]*G[i]) (here defined by string)
        contribution.setProfile(profile, xname="r")
        distribution="("
        for index in range(len(weights)):
            distribution+="%f"%weights[index]+"*G%d"%index+"+"
            
        distribution=distribution[:-1]+")"
        equation="s1*"+distribution
        print(equation+"\n")
        contribution.setEquation(equation) 
        
        # Create the recipe and add variables
        recipe = FitRecipe()
        recipe.addContribution(contribution)
        recipe.addVar(contribution.s1, SCALE_I, tag="scale")
        print("scale factor variable added to refinement \n")
        i=0
        for generator_cluster in generator_cluster_array:
            phase_cluster = generator_cluster.phase
            lattice = phase_cluster.getLattice()
            recipe.newVar("zoomscale_%d"%i, ZOOMSCALE_I, tag="lat")
            print("zoomscale_%d"%i+" variable added to refinement\n")
            recipe.constrain(lattice.a, 'zoomscale_%d'%i)
            recipe.constrain(lattice.b, 'zoomscale_%d'%i)
            recipe.constrain(lattice.c, 'zoomscale_%d'%i) 
            atoms=phase_cluster.getScatterers()
            recipe.newVar("Au_Uiso_%d"%i, UISO_Au_I, tag="adp")
            print("Au_UISO_%d"%i+" variable added to refinement\n")
            for atom in atoms:
                if atom.element.title() == "Au":
                    recipe.constrain(atom.Uiso, "Au_Uiso_%d"%i)
            recipe.addVar(generator_cluster.delta2, name="Au_Delta2_%d"%i, value=DELTA2_I, tag="d2")
            print("Au_Delta2_%d"%i+" variable added to refinement\n")
            i+=1
    else:
        print("Tables of weight must be equal to the number of models")
    return recipe