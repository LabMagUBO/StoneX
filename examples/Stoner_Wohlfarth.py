#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""
    Main program file example.

    StoneX V1.4.1 compatible.
    Date : 2 june 2016
"""


def init():
    """
        Import all required modules.
    """
    global sys, np, StoneX
    # Path to the module
    module_path = '/Users/zorg/These/Documents/Programmes/Python_Modules/StoneX_project/'

    # Define module path
    import sys
    sys.path.append(module_path)

    # Import numpy
    import numpy as np

    # Importing StoneX module
    import StoneX

    # Starting timer
    StoneX.timer(True)


def set_vsm():
    """
        Create the VSM from the StoneX module and return it.
    """
    # First, lets create the VSM
    vsm = StoneX.VSM()

    # Set the vsm parameters
    # Either B or H using vsm.B or vsm.H
    vsm.B = (100, 2, 65, 0.2, 'milli')     # B field = mu_0 * H
    vsm.phi = (0, 91, 5, 'deg')            # VSM direction
    vsm.T = (10, 50, 10, 'K')               # temperature sweep
    vsm._T = np.array([0.01])                # single or personalized array

    # Plotting control (True or False)
    vsm.plot_cycles = True
    vsm.plot_azimuthal = True
    vsm.plot_energyPath = False
    vsm.plot_energyLandscape = False    # Takes a lot of time and memory
    vsm.plot_T = False

    # Exporting data
    vsm.export_data = True

    # Displaying parameters
    StoneX.logger.info(vsm)

    return vsm


def set_domain():
    """
        Create a specific domain, with its parameters, and return it.
        Available models : Stoner_Wohlfarth, Meiklejohn_Bean, Garcia_Otero,
        Franco_Conde, Rotatable_AF, Double_MacroSpin
    """
    # Create the domain, with argument (model, 'folder_name')
    domain = StoneX.create_domain(
        StoneX.Stoner_Wohlfarth,
        'Stoner_Wohlfarth_example',
        mkdir=True
    )

    # Setting the domain's temperature
    # domain.T = 10

    # Theta discretization, with (step, 'unit')
    domain.theta = (0.1, 'deg')

    # Alpha step
    # domain.alpha = (2, 'deg')

    # Setting the physical parameters
    domain.Ms = 1740e3           # Magnetization in A/m, 1400kA/m for Co
    d_f = 8e-9                   # size of the particule, in meter
    domain.V_f = 4/3 * np.pi * (d_f / 2)**3
    domain.K_f = 50e3          # Uniaxial anisotropy, in J/m**3
    domain.gamma_f = 0         # Uniaxial anisotropy direction, in radians
    domain.K_bq = 0
    domain.gamma_bq = 0
    domain.K_iso = 0

    # AF parameters
    # domain.t_af = 0.5e-9             # Thickness of the AF shell, in m
    # d_af = d_f + domain.t_af * 2     # Global size of the F/AF core-shell, in m
    # domain.V_af = 4/3 * np.pi * ((d_af / 2)**3 - (d_f / 2)**3)   # AF volume
    # domain.J_ex = 10e-6  # 3e-7                      # Exchange energy
    # domain.S = 4 * np.pi * (d_f / 2)**2              # Exchange surface area
    # domain.K_af = 100e3                              # AF anisotropy

    StoneX.logger.info(domain)

    return domain


def set_sample(domain):
    """
        Creating a sample by combining multiple independent domains.
    """
    # Then we create a sample based on the domain
    # Distribution parameters
    N = 102                           # number of particules (average)
    # xMin = np.log(4)/np.log(10)
    # xMax = np.log(30) / np.log(10)
    xMin = 4
    xMax = 25
    nx = 100                               # number of domain
    # Arithmetic parameters
    mu = 8.29                              # mean particule diameter
    sigma = 2.99                           # 3nm standard deviation

    # LogNormal distribution
    # D = np.logspace(xMin, xMax, nx)
    D = np.linspace(xMin, xMax, nx)
    # D = np.linspace(-np.pi/2, np.pi/2, nx)              # flat distribution
    X, mu_log, sigma_log = StoneX.lognormale(D, mu, sigma)
    # X = np.zeros(D.size) + 1/N

    # Information about the distribution
    text = """Probability distribution : X
    Normale : mu = {}, sigma = {}
    Log-normale : m = {}, s = {}
    Effective Mean diameter: mean(X) = {}
    Number of particules : N = {}
    """.format(
        mu, sigma,
        mu_log, sigma_log,
        np.sum(D * X)/np.sum(X),
        np.sum(X * N)
        )
    StoneX.logger.info(text)

    # Creating the sample
    Density = N * X
    sample = StoneX.create_sample(domain, Density)

    # Change the domains' parameters according to the distribution
    # Domains data
    # diameter (nm), Density, V_f(m**3), V_af(m**3), Surface(m**2), Mag. Moment (A/m)
    domains_data = np.zeros((D.size, 6))
    domains_data[:, 0] = D          # Diameter (nm)
    R = D / 2 * 1e-9                # Radius (m)
    domains_data[:, 1] = Density
    domains_data[:, 2] = 4/3 * np.pi * R**3                     # V_f
    domains_data[:, 3] = 4/3 * np.pi * ((domain.t_af + R)**3 - R**3)   # V_af
    domains_data[:, 4] = 4 * np.pi * R**2                       # S
    domains_data[:, 5] = sample.domains[0].Ms                   # mu_s
    for i, d in enumerate(D):
        StoneX.logger.debug('Diameter {} : {}nm'.format(i, np.round(d, 2)))
        sample.domains[i].V_f = domains_data[i, 2]
        sample.domains[i].V_af = domains_data[i, 3]
        sample.domains[i].S = domains_data[i, 4]
        # logger.info(sample.domains[i])

    # Backing up the domain's values
    np.savetxt(
        '{}/domains_data.dat'.format(sample.name),
        domains_data,
        header='diameter (nm), Density, V_f(m**3), V_af(m**3), Surface(m**2), Mag. (A/m)'
    )

    return sample


def end(domain):
    """
        End processes
    """
    # Stopping timer
    StoneX.logger.info(StoneX.timer(False))

    # Saving the log file
    StoneX.shutil.copy(StoneX.main_file + '.log', domain.name + '/')


def main():
    """
        Main function.
    """
    # Importing modules
    init()

    # Setting the vsm
    vsm = set_vsm()

    # Set the domain
    domain = set_domain()

    # sample = set_sample(domain)

    # We can load into the VSM the sample or one domain only
    vsm.load(domain)

    vsm.measure()

    end(domain)


# If not imported, execute main()
if __name__ == '__main__':
    main()
