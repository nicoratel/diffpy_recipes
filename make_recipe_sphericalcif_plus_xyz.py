def make_recipe_sphericalcif_plus_xyz(cif_path1, stru_path, dat_path, anis_adp_Flag):
    """
    Creates and returns a Fit Recipe object with two phases.

    Parameters
    ----------
    cif_path1 : string, The full path to the structure CIF file to load, for the first phase.
    
    dat_path :  string, The full path to the PDF data to be fit.

    Returns
    ----------
    recipe :    The initialized Fit Recipe object using the datname and structure path
                provided.
    """
    # 9: Create two CIF file parsing objects, parse and load the structures, and
    # grab the space group names.
    p_cif1 = getParser('cif')
    
    stru1 = p_cif1.parseFile(cif_path1)
    
    sg1 = p_cif1.spacegroup.short_name
    
    stru2 = Structure(filename=str(stru_path))

    # 10: Create a Profile object for the experimental dataset and
    # tell this profile the range and mesh of points in r-space.
    profile = Profile()
    parser = PDFParser()
    parser.parseFile(dat_path)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)

    # 11a: Create a PDF Generator object for a periodic structure model
    # of phase 1.
    generator_crystal = DebyePDFGenerator("G1")
    generator_crystal.setStructure(stru1, periodic=True)
    generator_cluster = DebyePDFGenerator("G2")
    generator_cluster.setStructure(stru2, periodic=False)
    
    

    # 12: Create a Fit Contribution object. This is new, as we
    # need to tell the Fit Contribution about BOTH the phase
    # represented by 'generator_crystal1' AND the phase represented
    # by 'generator_crystal2'.
    contribution = FitContribution("crystal")
    contribution.addProfileGenerator(generator_crystal)
    contribution.addProfileGenerator(generator_cluster)
    
    
    # If you have a multi-core computer (you probably do), run your refinement in parallel!
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
        generator_crystal.parallel(ncpu=ncpu, mapfunc=pool.map)
        generator_cluster.parallel(ncpu=ncpu, mapfunc=pool.map)
        

    # 13: Set the Fit Contribution profile to the Profile object.
    contribution.setProfile(profile, xname="r")

    # 14: Set an equation, based on your PDF generators. This is
    # a more complicated case, since we have two phases. The equation
    # here will be the sum of the contributions from each phase,
    # 'G_Si' and 'G_Ni' weighted by a refined scale term for each phase,
    # 's1_Si' and '(1 - s1_Si)'. We also include a general 's2'
    # to account for data scale.
    contribution.registerFunction(sphericalCF, name="fsphere")
    contribution.setEquation("s1*(s2*G1*fsphere + (1.0-s2)*G2)")

    # 15: Create the Fit Recipe object that holds all the details of the fit.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # 16: Add, initialize, and tag the two scale variables.
    
    recipe.addVar(contribution.s1, DATA_SCALE_I, tag="scale")
    recipe.addVar(contribution.s2, FCC_SCALE, tag="scale")
    recipe.restrain("s1",
                    lb=0.0,
                    scaled=True,
                    sig=0.00001)

    recipe.restrain("s2",
                    lb=0.65,
                    ub=1.0,
                    scaled=True,
                    sig=0.00001)
    recipe.addVar(contribution.psize,30,tag="diameter")
    recipe.restrain("psize",lb=15, ub=100,scaled=True,sig=0.00001)
    # 17:This is new, we want to ensure that the data scale parameter 's2'
    # is always positive, and the phase scale parameter 's1_Si' is always
    # bounded between zero and one, to avoid any negative PDF signals.
    # We do this by adding 'restraints'. Effectively a restrain will modify
    # our objective function such that if the parameter approaches a user
    # defined upper or lower bound, the objective function will increase,
    # driving the fit away from the boundary.
    

    # 18a: This is a bit new. We will again use the srfit function
    # constrainAsSpaceGroup to constrain the lattice and ADP parameters
    # according to the space group of each of the two phases.
    # We loop through generators composed of PDF Generators
    # and space groups specific to EACH of the TWO candidate phases.
    # We use 'enumerate' to create an iterating index 'i' such that each
    # parameter can get it's own unique name, without colliding parameters.
    phase_cluster1 = generator_cluster.phase

    lattice1 = phase_cluster1.getLattice()

    recipe.newVar("zoomscale", ZOOM_SCALE_I, tag="lat")

    recipe.constrain(lattice1.a, 'zoomscale')
    recipe.constrain(lattice1.b, 'zoomscale')
    recipe.constrain(lattice1.c, 'zoomscale')
    
    for name, generator, space_group in zip(["Au"],
                                            [generator_crystal],
                                            [sg1]):

        # 18b: Initialize the instrument parameters, Q_damp and Q_broad, and
        # assign Q_max and Q_min for each phase.
        generator.qdamp.value = QDAMP_I
        generator.qbroad.value = QBROAD_I
        generator.setQmax(QMAX)
        generator.setQmin(QMIN)

        # 18c: Get the symmetry equivalent parameters for each phase.
        spacegroupparams = constrainAsSpaceGroup(generator.phase,
                                                 space_group)
        # 18d: Loop over and constrain these parameters for each phase.
        # Each parameter name gets the loop index 'i' appeneded so there are not
        # parameter name collisions.
        for par in spacegroupparams.latpars:
            recipe.addVar(par,
                          name=f"{par.name}_{name}",
                          fixed=False,
                          tag='lat')
        for par in spacegroupparams.adppars:
            recipe.addVar(par,
                          name=f"{par.name}_{name}",
                          fixed=False,
                          tag='adp')
            recipe.restrain(f"{par.name}_{name}",lb=0,ub=0.5,scaled=True,sig=0.00001)
       
     
    atoms1 = phase_cluster1.getScatterers()
    if anis_adp_Flag==True:
        recipe.newVar("Au_U11",Au_U11,tag="adp");recipe.newVar("Au_U12",Au_U12,tag="adp");recipe.newVar("Au_U13",Au_U13,tag="adp")
        recipe.newVar("Au_U22",Au_U22,tag="adp");recipe.newVar("Au_U23",Au_U23,tag="adp")
        recipe.newVar("Au_U33",Au_U33,tag="adp")
        recipe.restrain("Au_U11",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Au_U12",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Au_U13",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Au_U22",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Au_U23",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Au_U33",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        for atom in atoms1:
            if atom.element.title()=="Au":
                recipe.constrain(atom.U11,"Au_U11");recipe.constrain(atom.U12,"Au_U12");recipe.constrain(atom.U13,"Au_U13")
                recipe.constrain(atom.U22,"Au_U22");recipe.constrain(atom.U23,"Au_U23")
                recipe.constrain(atom.U33,"Au_U33")
        
       
    if anis_adp_Flag==False:
        recipe.newVar("Au_Uiso", UISO_Au_I, tag="adp")
        for atom in atoms1:
            if atom.element.title() == "Au":
                recipe.constrain(atom.Uiso, "Au_Uiso")
        # 19: Add delta, but not instrumental parameters to Fit Recipe.
        # One for each phase.
    recipe.addVar(generator_crystal.delta2, name="Delta2_crystal",
                      value=DELTA2_I, tag="d2")
    recipe.addVar(generator_cluster.delta2, name="Delta2_cluster",
                      value=DELTA2_I, tag="d2")
    recipe.restrain("Delta2_crystal",lb=0.0,ub=20.0,scaled=True,sig=0.00001)
    recipe.restrain("Delta2_cluster",lb=0.0,ub=20.0,scaled=True,sig=0.00001)
        #recipe.restrain(f"Delta2_{name}",
        #                lb=0.0,
        #                ub=5.0,
        #                scaled=True,
        #                sig=0.00001)
    #19 bis define shape function for G(r) calculation using sphercialCF from diffpy
    #recipe.
    recipe.crystal.registerFunction(sphericalCF,name='recipe')
    # 20: Return the Fit Recipe object to be optimized.
    return recipe

    # End of function
