def plot_results(recipe, figname):
    """
    Creates plots of the fitted PDF and residual, and writes them to disk
    as *.pdf files.

    Parameters
    ----------
    recipe :    The optimized Fit Recipe object containing the PDF data
                we wish to plot
    figname :   string, the location and name of the figure file to create

    Returns
    ----------
    None
    """
    r = recipe.cluster.profile.x

    g = recipe.cluster.profile.y
    gcalc = recipe.cluster.profile.ycalc
    diffzero = -0.65 * max(g) * np.ones_like(g)
    diff = g - gcalc + diffzero

    mpl.rcParams.update(mpl.rcParamsDefault)
    #plt.style.use(str(PWD.parent.parent.parent / "utils" / "billinge.mplstyle"))

    fig, ax1 = plt.subplots(1, 1)

    ax1.plot(r,
             g,
             ls="None",
             marker="o",
             ms=5,
             mew=0.2,
             mfc="None",
             label="G(r) Data")

    ax1.plot(r, gcalc, lw=1.3, label="G(r) Fit")
    ax1.plot(r, diff, lw=1.2, label="G(r) diff")
    ax1.plot(r, diffzero, lw=1.0, ls="--", c="black")

    ax1.set_xlabel(r"r($\mathrm{\AA}$)")
    ax1.set_ylabel(r"G($\mathrm{\AA}$$^{-2}$)")
    ax1.tick_params(axis="both",
                    which="major",
                    top=True,
                    right=True)

    ax1.set_xlim(r[0], r[-1])
    ax1.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(figname.parent / f"{figname.name}.pdf", format="pdf")

    # End of function