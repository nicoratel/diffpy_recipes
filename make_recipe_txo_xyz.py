def make_recipe_two_xyz(stru_path1,stru_path2, dat_path,anis_adp_Flag,fit_Qdamp_flag=False):
    """
    Creates and returns a Fit Recipe object

    Parameters
    ----------
    stru_path : string, The full path to the structure XYZ file to load.
    dat_path :  string, The full path to the PDF data to be fit.

    Returns
    ----------
    fitrecipe : The initialized Fit Recipe object using the datname and structure
                provided.
    """

    stru1 = Structure(filename=str(stru_path1))
    stru2 = Structure(filename=str(stru_path2))
    if anis_adp_Flag==True:
        stru1.anisotropy = True
        stru2.anisotropy=True
    # 9: Create a Profile object for the experimental dataset and
    # tell this profile the range and mesh of points in r-space.
    profile = Profile()
    parser = PDFParser()
    parser.parseFile(dat_path)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)

    # 10: Create a Debye PDF Generator object for the discrete structure model.
    generator_cluster1 = DebyePDFGenerator("G1")
    #generator_cluster1 = PDFGenerator("G1")
    generator_cluster1.setStructure(stru1, periodic=False)
    generator_cluster2 = DebyePDFGenerator("G2")
    #generator_cluster1 = PDFGenerator("G1")
    generator_cluster2.setStructure(stru2, periodic=False)
    # 11: Create a Fit Contribution object.
    contribution = FitContribution("cluster")
    contribution.addProfileGenerator(generator_cluster1)
    contribution.addProfileGenerator(generator_cluster2)
    # If you have a multi-core computer (you probably do),
    # run your refinement in parallel!
    # Here we just make sure not to overload your CPUs.
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
        generator_cluster1.parallel(ncpu=ncpu, mapfunc=pool.map)
    # 12: Set the Fit Contribution profile to the Profile object.
    contribution.setProfile(profile, xname="r")

    # 13: Set an equation, based on your PDF generators. 
    contribution.setEquation("s2*(s1*G1 + (1.0-s1)*G2)")

    # 14: Create the Fit Recipe object that holds all the details of the fit.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # 15: Initialize the instrument parameters, Q_damp and Q_broad, and
    # assign Q_max and Q_min.
    generator_cluster1.qdamp.value = QDAMP_I
    generator_cluster1.qbroad.value = QBROAD_I
    generator_cluster1.setQmax(QMAX)
    generator_cluster1.setQmin(QMIN)
    generator_cluster2.qdamp.value = QDAMP_I
    generator_cluster2.qbroad.value = QBROAD_I
    generator_cluster2.setQmax(QMAX)
    generator_cluster2.setQmin(QMIN)

    # 16: Add, initialize, and tag variables in the Fit Recipe object.
    # In this case we also add psize, which is the NP size.
    recipe.addVar(contribution.s1, SCALE_5shell, tag="scale")
    recipe.addVar(contribution.s2, DATA_SCALE, tag="scale")
    # 16b:This is new, we want to ensure that the data scale parameter 's2'
    # is always positive, and the phase scale parameter 's1' is always
    # bounded between zero and one, to avoid any negative PDF signals.
    # We do this by adding 'restraints'. Effectively a restrain will modify
    # our objective function such that if the parameter approaches a user
    # defined upper or lower bound, the objective function will increase,
    # driving the fit away from the boundary.
    recipe.restrain("s2",
                    lb=0.0,
                    scaled=True,
                    sig=0.00001)

    recipe.restrain("s1",
                    lb=0.0,
                    ub=1.0,
                    scaled=True,
                    sig=0.00001)

    # 17: Define a phase and lattice from the Debye PDF Generator
    # object and assign an isotropic lattice expansion factor tagged
    # "zoomscale" to the structure. 

    phase_cluster1 = generator_cluster1.phase

    lattice1 = phase_cluster1.getLattice()

    recipe.newVar("zoomscale", ZOOMSCALE_I, tag="lat")

    recipe.constrain(lattice1.a, 'zoomscale')
    recipe.constrain(lattice1.b, 'zoomscale')
    recipe.constrain(lattice1.c, 'zoomscale')
    phase_cluster2 = generator_cluster2.phase

    lattice2 = phase_cluster2.getLattice()

    #recipe.newVar("zoomscale1", ZOOMSCALE_I, tag="lat")

    recipe.constrain(lattice2.a, 'zoomscale')
    recipe.constrain(lattice2.b, 'zoomscale')
    recipe.constrain(lattice2.c, 'zoomscale')
    # 18: Initialize an atoms object and constrain the isotropic
    # Atomic Displacement Paramaters (ADPs) per element. 

    atoms1 = phase_cluster1.getScatterers()
    if anis_adp_Flag==True:
        recipe.newVar("Ag_U11",Ag_U11,tag="adp");recipe.newVar("Ag_U12",Ag_U12,tag="adp");recipe.newVar("Ag_U13",Ag_U13,tag="adp")
        recipe.newVar("Ag_U22",Ag_U22,tag="adp");recipe.newVar("Ag_U23",Ag_U23,tag="adp")
        recipe.newVar("Ag_U33",Ag_U33,tag="adp")
        recipe.restrain("Ag_U11",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Ag_U12",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Ag_U13",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Ag_U22",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Ag_U23",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        recipe.restrain("Ag_U33",lb=0.0,ub=1.0,scaled=True,sig=0.00001)
        for atom in atoms1:
            if atom.element.title()=="Ag":
                recipe.constrain(atom.U11,"Ag_U11");recipe.constrain(atom.U12,"Ag_U12");recipe.constrain(atom.U13,"Ag_U13")
                recipe.constrain(atom.U22,"Ag_U22");recipe.constrain(atom.U23,"Ag_U23")
                recipe.constrain(atom.U33,"Ag_U33")
    
   
    if anis_adp_Flag==False:
        recipe.newVar("Ag_Uiso", UISO_Ag_I, tag="adp")
        for atom in atoms1:
            if atom.element.title() == "Ag":
                recipe.constrain(atom.Uiso, "Ag_Uiso")

    # 19: Add and tag a variable for correlated motion effects, and Q damp
    recipe.addVar(generator_cluster1.delta2, name="Au_Delta2", value=DELTA2_I, tag="d2")
    if fit_Qdamp_flag:
        recipe.addVar(generator_cluster1.qdamp,
                      fixed=False,
                      name="Qdamp",
                      value=QDAMP_I,
                      tag="inst")
    return recipe
