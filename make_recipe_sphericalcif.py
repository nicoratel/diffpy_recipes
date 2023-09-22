def make_recipe(self,cif_path1, dat_path):
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
    
    PDF_RMIN=float(self.ui.rmin_lineEdit.text())
    PDF_RMAX=float(self.ui.rmax_lineEdit.text())
    PDF_RSTEP=float(self.ui.rstep_lineEdit.text())
    QBROAD_I=float(self.ui.qbroad_lineEdit.text())
    QDAMP_I=float(self.ui.qdamp_lineEdit.text())
    QMIN=float(self.ui.qmin_lineEdit.text())
    QMAX=float(self.ui.qmax_lineEdit.text())
    RUN_PARALLEL=True
    # 10: Create a Profile object for the experimental dataset and
    # tell this profile the range and mesh of points in r-space.
    profile = Profile()
    parser = PDFParser()
    parser.parseFile(dat_path)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)

    # 11a: Create a PDF Generator object for a periodic structure model
    # of phase 1.
    generator_crystal1 = DebyePDFGenerator("G1")
    generator_crystal1.setStructure(stru1, periodic=True)

    
    

    # 12: Create a Fit Contribution object. This is new, as we
    # need to tell the Fit Contribution about BOTH the phase
    # represented by 'generator_crystal1' AND the phase represented
    # by 'generator_crystal2'.
    contribution = FitContribution("crystal")
    contribution.addProfileGenerator(generator_crystal1)
    
    
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
        generator_crystal1.parallel(ncpu=ncpu, mapfunc=pool.map)
        

    # 13: Set the Fit Contribution profile to the Profile object.
    contribution.setProfile(profile, xname="r")

    # 14: Set an equation, based on your PDF generators. This is
    # a more complicated case, since we have two phases. The equation
    # here will be the sum of the contributions from each phase,
    # 'G_Si' and 'G_Ni' weighted by a refined scale term for each phase,
    # 's1_Si' and '(1 - s1_Si)'. We also include a general 's2'
    # to account for data scale.
    if self.ui.sphericalCF_checkBox.isChecked():
        contribution.registerFunction(sphericalCF, name="fsphere")
        contribution.setEquation("s1*G1*fsphere")
        
    else:
        contribution.setEquation("s1*G1")
    # 15: Create the Fit Recipe object that holds all the details of the fit.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # 16: Add, initialize, and tag the two scale variables.
    
    recipe.addVar(contribution.s1, float(self.ui.scaleFactor_lineEdit.text()), tag="scale")
    if self.ui.sphericalCF_checkBox.isChecked():
        recipe.addVar(contribution.psize,float(self.ui.diameter_lineEdit.text()),tag="radius")
        print("psize variable added to recipe")
   
    # 18a: This is a bit new. We will again use the srfit function
    # constrainAsSpaceGroup to constrain the lattice and ADP parameters
    # according to the space group of each of the two phases.
    # We loop through generators composed of PDF Generators
    # and space groups specific to EACH of the TWO candidate phases.
    # We use 'enumerate' to create an iterating index 'i' such that each
    # parameter can get it's own unique name, without colliding parameters.
    for name, generator, space_group in zip([""],
                                            [generator_crystal1],
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
        # 19: Add delta, but not instrumental parameters to Fit Recipe.
        # One for each phase.
        recipe.addVar(generator.delta2, name=f"Delta2_{name}",
                     value=float(self.ui.delta2_lineEdit.text()), tag="d2")
        recipe.restrain(f"Delta2_{name}",lb=0,ub=10,scaled=True,sig=0.00001)
        #recipe.restrain(f"Delta2_{name}",
        #                lb=0.0,
        #                ub=5.0,
        #                scaled=True,
        #                sig=0.00001)
    #19 bis define shape function for G(r) calculation using sphercialCF from diffpy
   
    recipe.crystal.registerFunction(sphericalCF,name='recipe')
    # 20: Return the Fit Recipe object to be optimized.
    return recipe
