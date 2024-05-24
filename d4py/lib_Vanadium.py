#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:17:33 2022

@author: cuello
"""

import numpy as np  # NumPy library
import matplotlib.pyplot as plt  # For plotting
import lmfit as lm  # For fitting
import d4treat as d4
from scipy.interpolate import interp1d
from Element import elemento as elem  # specific Class created for load atoms parameters
from readParamD4 import generate_correction_factor, atomDisarm
import mc  # MonteCarlo Module


# --------1---------2---------3---------4---------5---------6---------7---------
class vanaObj:
    """
    Type: Class

    Object:
        fit vanadium data, for them first perform background subtraction

    Input:
        params: parameters obtained from the parameter txt file

    Output:
        An instance created with the following attributes and methods.

        self.beam (dict) Basic information about the beam
        self.can (dict) Basic information about the container
        self.binning (dict) Definition of the binning
        self.qmax
        self.vana (dict) Basic information about the vanadium
        self.environmentData (array) environment data as 3-column table
        self.containerData (array) container data as 3-column table
        self.vanadiumData (array) vanadium data as 3-column table
        self.environmentData_qbin (array) rebinned environment data as 3-column table
        self.containerData_qbin (array) rebinned container data as 3-column table
        self.vanadiumData_qbin (array) rebinned vanadium data as 3-column table

    Methods:

        __init__(self,params):
            load all the necessary information

        plotRaw(self):
            Only plot Raw and Rebinned data

        setting_bs_mc_d4vfit(self,params):
            calculate the background subtraction by the montecarlo method

        plot_bs_mc_results(self,energias,XS_model_sample,XS_model_can,xs_sample,xs_can,E0,qvals,single,s_aten,can,total):
            plot montecarlo results (is used by setting_bs_mc_d4vfit function)

        setting_model_d4vfit(self,a0Val=True,a1Val=True,a2Val=True,a0 = 1,a1 = 0.001,a2 = -0.0001,
                            Q0=25,dQ=10,x_fit_min=1,x_fit_max=24.5,y_fit_min=9780,y_fit_max=14000):
            allows fitting the vanadium data with a sigmoidal function multiplied by a polynomial of up to second order.

        plot_fit_results(self,vanaInc,fitdata,x,y,vana_i,vana_n,vana_inel):
            only plot fit results

    Author: Matías Ilarri
    Created: 01/01/2023
    Modified:18/01/2023
    """

    def __init__(self, params):

        # Basic information about the beam
        self.beam = d4.setBeamInfo(
            zeroAngle=params.zeroAngle,
            wavelength=params.wavelength,
            LohenSlit=params.LohenSlit,
            GamsSlit=params.GamsSlit,
            topFlag=params.topFlag,
            bottomFlag=params.bottomFlag,
        )
        # --------1---------2---------3---------4---------5---------6---------7-----
        # help(d4.setBeamInfo)
        # --------1---------2---------3---------4---------5---------6---------7-----
        # Basic information about the container
        self.can = d4.setCanInfo(
            material=params.material,
            shape=params.shape,
            outerDiam=params.outerDiam,
            innerDiam=params.innerDiam,
            height=params.height,
        )
        # --------1---------2---------3---------4---------5---------6---------7-----
        # Definition of the binning
        self.binning = d4.setBinInfo(
            AngularResolution=params.AngularResolution,
            AMin=params.AMin,
            AMax=params.AMax,
            AStep=params.AStep,
            QMin=params.QMin,
            QMax=params.QMax,
            QStep=params.QStep,
            RMin=params.RMin,
            RMax=params.RMax,
            RStep=params.RStep,
        )

        self.qmax = self.binning["Qbin"][1]

        auxElement = elem("V", table=0)
        auxVana = auxElement.__dict__

        self.vana = d4.setVanaInfo(
            IncXS=auxVana["sig_inc"],
            CohXS=auxVana["sig_coh"],
            ScaXS=auxVana["sig_sca"],
            AbsXS=auxVana["sig_abs"],
            CohSL=auxVana["re_bcoh"],
            molarM=auxVana["weight"],
            NAtoms=1,
            Diam=params.Diam,
            Height=self.beam["height"],
            density=params.density,
        )

        self.environmentData = d4.read_3col(params.environment)
        self.containerData = d4.read_3col(params.container)
        self.vanadiumData = d4.read_3col(params.vanadium)

        self.environmentData_qbin = d4.rebin(
            self.binning["Ares"],
            self.beam["wlength"],
            self.environmentData,
            *self.binning["Qbin"]
        )
        self.containerData_qbin = d4.rebin(
            self.binning["Ares"],
            self.beam["wlength"],
            self.containerData,
            *self.binning["Qbin"]
        )
        self.vanadiumData_qbin = d4.rebin(
            self.binning["Ares"],
            self.beam["wlength"],
            self.vanadiumData,
            *self.binning["Qbin"]
        )

        self.plotRaw()

        return

    def plotRaw(self):

        plt.figure(figsize=(9, 6))

        plt.subplot(2, 1, 1)
        plt.plot(
            self.environmentData[:, 0],
            self.environmentData[:, 1],
            "r-",
            label="Environment",
        )
        plt.plot(
            self.containerData[:, 0], self.containerData[:, 1], "b-", label="Container"
        )
        plt.plot(
            self.environmentData_qbin[:, 0],
            self.environmentData_qbin[:, 1],
            "g-",
            label="Environment(Rebinned)",
        )
        plt.plot(
            self.containerData_qbin[:, 0],
            self.containerData_qbin[:, 1],
            "y-",
            label="Container(Rebinned)",
        )

        plt.ylabel("$G(r)$")
        plt.legend(loc="best")
        plt.axis([None, None, None, None])
        plt.grid(True)
        plt.legend(loc="best")
        plt.title("Input data")
        plt.ylabel("Intensity (arb. units)")

        plt.subplot(2, 1, 2)
        plt.plot(
            self.vanadiumData[:, 0], self.vanadiumData[:, 1], "m-", label="Vanadium"
        )
        plt.plot(
            self.vanadiumData_qbin[:, 0],
            self.vanadiumData_qbin[:, 1],
            "y-",
            label="Vanadium(Rebinned)",
        )
        plt.ylabel("$g(r)$")
        plt.legend(loc="best")
        plt.axis([None, None, None, None])
        plt.grid(True)
        plt.legend(loc="best")
        plt.xlabel(r"$Q (1/\AA)$")
        plt.ylabel("Intensity (arb. units)")

        plt.show()
        return

    def setting_bs_mc_d4vana(self, params):

        # ----------Variable Area------------
        # -----------------------------------
        E0 = 327  # meV
        Eeff = 25.3  # meV
        # -----------------------------------
        # -----------------------------------

        folder_generated = "generated_files"
        # Creating and saving the vanadium background-subtracted diffractogram.
        self.van_qbin_sub = d4.wsum2(
            1.0, self.vanadiumData_qbin, -1.0, self.environmentData_qbin
        )

        # ia=0
        # while self.van_qbin_sub[ia,1] == 0:
        #     ia=ia+1
        # j=ia
        # while self.van_qbin_sub[j,1] < np.mean(self.van_qbin_sub[ia:,1]):
        #     j=j+1

        n_points = 400
        # sample
        x = self.van_qbin_sub[:, 0]
        y = self.van_qbin_sub[:, 1]
        f = interp1d(x, y)

        new_x = np.linspace(x.min(), x.max(), n_points)
        new_y = f(new_x)

        # can
        # x_bj, y_bj = self.environmentData_qbin[:, 0], self.environmentData_qbin[:, 1]
        x_can, y_can = self.containerData_qbin[:, 0], np.ones(
            len(self.environmentData_qbin[:, 1])
        )
        #         x_canB, y_canB = self.containerData_qbin[:,0],self.containerData_qbin[:,1]

        #         x_can = x_canB
        #        y_can = y_canB - y_bj

        f2 = interp1d(x_can, y_can)
        new_y_c = f2(new_x)

        header = """Q    DQMUESTRA    DQCAN   (distrib. muestra-can en Q)   , iteracion 9
2.6371248E-02  6.7652303E-04 -5.0012457E-05
3.6499664E-02  5.6225456E-04 -1.8670474E-05"""
        np.savetxt(
            folder_generated + "/" + "Mcd8dq.dat",
            np.array((new_x, new_y, new_y_c)).T,
            fmt="%.7e",
            header=header,
            delimiter="  ",
        )

        cvals = []
        sca_XS = []
        avals = []
        abs_XS = []

        U0 = np.log10(E0 * 1e-3)

        can_elements = {}
        comp_sample = {}
        comp_sample[params.material] = 1
        can_atoms_sample = atomDisarm(comp_sample)

        can_natoms_sample = d4.getNofAtoms(can_atoms_sample)
        can_c_sample = {}
        for key in can_atoms_sample.keys():
            can_element = elem("{}".format(key), table=0)
            can_elements[can_element.symbol] = can_element.__dict__

        for key, value in can_atoms_sample.items():
            can_elements[key]["MassNbr"] = d4.getMassNumber(
                float(can_elements[key]["weight"])
            )
            can_c_sample[key] = value / can_natoms_sample

        for key in can_elements:
            cvals.append(can_c_sample[key])
            sca_XS.append(can_elements[key]["sig_sca"])
            avals.append(can_elements[key]["weight"])
            abs_XS.append(can_elements[key]["sig_abs"])

        xs_sample = d4.XS_model(
            Eval=E0, Eeff=Eeff, composition_vec=cvals, bound_xs_vec=sca_XS, A_vec=avals
        )

        can_cvals = cvals
        can_sca_XS = sca_XS
        can_avals = avals
        can_abs_XS = abs_XS

        xs_can = xs_sample

        mac_xs_sample = (
            self.vana["den_aac"] * xs_can
        )  # va effden_acc? atomos/A-3 * barns / atomos
        mac_xs_van = mac_xs_sample

        p_scat_sample = d4.scattering_probability(E0, xs_sample, cvals, abs_XS)
        p_scat_can = d4.scattering_probability(E0, xs_can, can_cvals, can_abs_XS)

        # print(U0, xs_sample, mac_xs_sample, p_scat_sample, xs_can, mac_xs_van, p_scat_can)

        energias = np.linspace(1, 400, 200)
        # XS_model_sample = [
        #     d4.XS_model(
        #         Eval=energia,
        #         Eeff=Eeff,
        #         composition_vec=cvals,
        #         bound_xs_vec=sca_XS,
        #         A_vec=avals,
        #     )
        #     for energia in energias
        # ]

        XS_model_can = [
            d4.XS_model(
                Eval=energia,
                Eeff=Eeff,
                composition_vec=can_cvals,
                bound_xs_vec=can_sca_XS,
                A_vec=can_avals,
            )
            for energia in energias
        ]

        # ----------Variable Area------------
        # -----------------------------------
        max_count = 6  # number of iterations of Monte carlo routine
        ncic = 5
        nt = 1000
        r2 = float(params.Diam) / 2
        r1 = r2 - 0.00001  # interno del container
        h = 5
        dh = 0.01
        ah = 5
        # -----------------------------------
        # -----------------------------------
        # print(mc.mc2022.__doc__)

        new_y2 = new_y.copy()
        print(13 * "=+" + " " + "Waiting" + " " + 13 * "+=")
        for i in range(0, max_count):
            # np.savetxt('generated_files/Mcd8dq.dat', np.array((new_x, new_y2, new_y_c)).T, fmt='%.7e',
            #            header= header, delimiter='  ')
            angulo, single, s_aten, can, total = mc.mc2022(
                "generated_files/salida",
                ncic=ncic,
                nt=nt,
                e0=E0,
                utab=U0,
                setot=xs_sample,
                sem=mac_xs_sample,
                pscat=p_scat_sample,
                setotc=xs_can,
                semc=mac_xs_van,
                pscatc=p_scat_can,
                r1=r1 / 10,
                r2=r2 / 10,
                h=h / 10,
                dh=dh,
                ah=ah,
            )
            qvals = d4.ang2q(angulo, wlength=float(params.wavelength))
            global_factor = generate_correction_factor(
                new_x, qvals, single[:, 0], s_aten, total
            )
            new_y2 = new_y * global_factor
            print(15 * "=+" + " " + str(i) + " " + 15 * "+=")

        print(60 * "-")
        print("MonteCarlo Simulation is finished")
        print(60 * "-")

        # Creating and saving the vanadium background-subtracted diffractogram.
        x = new_x
        y = new_y2
        f = interp1d(x, y)

        self.van_qbin_sub_mc = self.van_qbin_sub
        self.van_qbin_sub_mc[:, 1] = f(self.van_qbin_sub[:, 0])

        heading = [
            "Background subtracted vanadium data",
            "result =  vanadium - MTbellar",
            "  Q(1/A)       Intensity            Error",
        ]
        file = "vana_qbin_sub.dat"
        d4.saveFile_3col(file, self.van_qbin_sub_mc, heading)
        print(70 * "-")
        self.plot_bs_mc_results(
            energias, XS_model_can, xs_can, E0, qvals, single, s_aten, can, total
        )
        return

    def setting_model_d4vana(
        self,
        a0Val=True,
        a1Val=False,
        a2Val=True,
        a3Val=False,
        a4Val=True,
        a0=1,
        a1=0.0,
        a2=0.001,
        a3=0.0,
        a4=-0.0001,
        x_fit_min=1,
        x_fit_max=24.5,
        y_fit_min=9780,
        y_fit_max=14000,
    ):
        # ----------Variable Area------------
        # -----------------------------------

        #         a0Val=True               #if xxVal es "True" then the code adjunst xx parameter
        #         a1Val=True               #if xxVal is "False" then the code does not adjust xx parameter
        #         a2Val=True
        #         a0 = 1                   # the 3 parameters of a second order polynomial
        #         a1 = 0.001
        #         a2 = -0.0001
        #         x_fit_min=1
        #         x_fit_max=24.5
        #         y_fit_min=9780
        #         y_fit_max=14000
        # -----------------------------------
        # -----------------------------------

        print()
        print(10 * "-", " Vanadium Q-dependence ", 10 * "-")

        # A=self.vana['MassNbr']
        # i=0
        # while self.van_qbin_sub_mc[i,1] == 0:
        #     i=i+1
        # lowQ=np.mean(self.van_qbin_sub_mc[i:,1])
        # --------1---------2---------3---------4---------5---------6---------7-----
        # Definition of the model to fit vanadium data in Q-scale
        # A: the mass number for vanadium
        # the 3 parameters of a second order polynomial
        # then lowQ limit, the inflection point of the sigmoidal and the width of this function
        vanaInc = d4.polyQ4(self.van_qbin_sub_mc[:, 0], a0, a1, a2, a3, a4)
        # --------1---------2---------3---------4---------5---------6---------7-----
        # Here we are defining the limits of the fit.
        # We provide x1,x2,y1,y2, defining a rectangle.
        # Only data inside this rectangle are used for the fitting
        #
        x = [x_fit_min, x_fit_max]
        y = [y_fit_min, y_fit_max]
        fitdata = d4.fit_range(x[0], x[1], y[0], y[1], self.van_qbin_sub_mc)

        print()
        print(10 * "-", " Fitting vanadium data ", 10 * "-")

        fitdata = d4.fit_range(
            x_fit_min, x_fit_max, y_fit_min, y_fit_max, self.van_qbin_sub_mc
        )
        # fitdata = d4.fit_range(0.36,23,15000,18000,van_qbin_sub)
        gmodel = lm.Model(d4.polyQ4)
        print("parameter names: {}".format(gmodel.param_names))
        print("independent variables: {}".format(gmodel.independent_vars))
        gmodel.make_params(a0=a0, a1=a1, a2=a2, a3=a3, a4=a4)
        # gmodel.set_param_hint('A',vary=False)
        # gmodel.set_param_hint('lowQ',vary=True)
        # gmodel.set_param_hint('Q0',vary=True)
        # gmodel.set_param_hint('dQ',vary=True)
        gmodel.set_param_hint("a0", vary=a0Val)
        gmodel.set_param_hint("a1", vary=a1Val)
        gmodel.set_param_hint("a2", vary=a2Val)
        gmodel.set_param_hint("a3", vary=a3Val)
        gmodel.set_param_hint("a4", vary=a4Val)
        result = gmodel.fit(
            fitdata[:, 1], x=fitdata[:, 0], a0=a0, a1=a1, a2=a2, a3=a3, a4=a4
        )
        print(result.fit_report())
        # print(result.best_values)

        # vana_inel is a tuple with the results of the fitted parameters of the sigmoidal function
        # vana_inel = (result.best_values['A'],result.best_values['lowQ'],result.best_values['Q0'],
        #              result.best_values['dQ'])
        # vana_res is a tuple with the results of the fitted parameters of the polynomial
        vana_fit = (
            result.best_values["a0"],
            result.best_values["a1"],
            result.best_values["a2"],
            result.best_values["a3"],
            result.best_values["a4"],
        )

        # Total incoherent (or self) contribution.
        # This includes the sigmoidal and the polynomial
        vanaInc = d4.polyQ4(self.van_qbin_sub_mc[:, 0], *vana_fit)

        # This is the inelastic effect described by the sigmoidal function
        # vana_i = d4.inelastic(self.van_qbin_sub_mc[:,0],*vana_inel)

        # The first element of the tuple is the lowQ parameter
        print("----> lowQ inelasticity = {:10.3f}".format(vana_fit[0]))
        print(
            "----> Model at Q={:4.2f} 1/A = {:10.3f}".format(
                self.van_qbin_sub_mc[0, 0], vanaInc[0]
            )
        )
        print(
            "----> Factor for resolution = {:10.3f}".format(
                vanaInc[0] / result.best_values["a0"]
            )
        )

        # This is the contribution of the instrument resolution function and other effects
        vana_r = d4.polyQ4(self.van_qbin_sub_mc[:, 0], *vana_fit)
        # Here this contribution is renormalised to the low-Q limit of the fit.
        vana_n = vanaInc[0] / result.best_values["a0"] * vana_r
        # print(vanaInc[0]/result.best_values['a0'])

        # --------1---------2---------3---------4---------5---------6---------7-----
        # Renormalisation of the vanadium data, thus the scale is now absolute, i.e., barns/sterad/atom
        print()
        print(10 * "-", " Renormalisation of the vanadium data ", 10 * "-")

        # y_van_qbin_sub_nor = np.array(y_van_qbin_sub) * vana['SelfQ0'] /vana_n
        # e_van_qbin_sub_nor = np.array(e_van_qbin_sub) * vana['SelfQ0'] /vana_n

        self.van_qbin_sub_nor = self.van_qbin_sub_mc.copy()
        self.van_qbin_sub_nor[:, 1] = (
            self.van_qbin_sub_mc[:, 1] * self.vana["SelfQ0"] / vana_n
        )
        self.van_qbin_sub_nor[:, 2] = (
            self.van_qbin_sub_mc[:, 2] * self.vana["SelfQ0"] / vana_n
        )

        print(
            "Low-Q limit for vanadium= {:.6f} barns/sterad/atom".format(
                self.vana["SelfQ0"]
            )
        )
        print(
            "High-Q limit for vanadium= {:.6f} barns/sterad/atom".format(
                self.vana["FreeIncXS"] / 4 / np.pi
            )
        )

        self.plot_fit_results(vanaInc, fitdata, x, y, vana_r)

        heading = [
            "Renormalisation of the vanadium data",
            "result =  vanadium",
            "  Q(1/A)       Intensity            Error",
        ]
        file = "van_qbin_sub_nor.dat"
        d4.saveFile_3col(file, self.van_qbin_sub_nor, heading)

        # Saving vana_n.
        fp = open("savedVana_nToFit.txt", "w")
        for value in vana_n:
            fp.write(str(value) + " ")
        fp.close()

        return

    def plot_fit_results(self, vanaInc, fitdata, x, y, vana_inel):
        i = 0
        while self.van_qbin_sub_mc[i, 1] == 0:
            i = i + 1

        # Plot of the experimental data and the proposed model
        plt.figure(figsize=(8, 5))
        plt.plot(
            self.van_qbin_sub_mc[i:, 0],
            self.van_qbin_sub_mc[i:, 1],
            "g-",
            label="(V+E)*Fmc",
        )
        plt.plot(self.van_qbin_sub_mc[i:, 0], vanaInc[i:], "r-", label="Initial model")
        plt.legend(loc="best")
        plt.title("Vanadium: Rebinned data in $Q$ scale")
        plt.xlabel(r"$Q (1/\AA)$")
        plt.ylabel("Intensity (arb. units)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(8, 5))
        plt.plot(
            self.van_qbin_sub_mc[i:, 0],
            self.van_qbin_sub_mc[i:, 1],
            "bo",
            label="experimental data",
        )
        plt.plot(fitdata[:, 0], fitdata[:, 1], "y+", label="fitted data")
        plt.axvline(x=x[0], color="r", linestyle="--")
        plt.axvline(x=x[1], color="r", linestyle="--")
        plt.axhline(y=y[0], color="r", linestyle="--")
        plt.axhline(y=y[1], color="r", linestyle="--")
        plt.legend(loc="best")
        plt.title("Vanadium: Fitting the self scattering")
        plt.xlabel(r"$Q (1/\AA)$")
        plt.ylabel("Intensity (arb. units)")
        plt.grid(True)
        plt.show()

        plt.figure(figsize=(8, 5))
        plt.plot(
            self.van_qbin_sub_mc[i:, 0],
            self.van_qbin_sub_mc[i:, 1],
            "bo",
            label="experimental data",
        )
        plt.plot(fitdata[:, 0], fitdata[:, 1], "y+", label="fitted data")
        plt.plot(self.van_qbin_sub_mc[:, 0], vanaInc, "r-", label="best fit")
        # plt.plot(self.van_qbin_sub_mc[:,0], vana_i, 'g-', label='inelastic component')
        # plt.plot(self.van_qbin_sub_mc[:,0], vana_n, 'm-', label='normalisation data')
        plt.legend(loc="best")
        plt.title("Vanadium: Fitting the self scattering")
        plt.xlabel(r"$Q (1/\AA)$")
        plt.ylabel("Intensity (arb. units)")
        plt.grid(True)
        plt.show()

        plt.figure(figsize=(8, 5))
        plt.errorbar(
            self.van_qbin_sub_nor[i:, 0],
            self.van_qbin_sub_nor[i:, 1],
            yerr=None,
            label="Vanadium",
        )
        plt.plot(
            self.van_qbin_sub_nor[:, 0],
            vanaInc * self.vana["SelfQ0"] / vana_inel[1],
            "r-",
            label="Model",
        )
        plt.legend(loc="best")
        plt.title("Vanadium: Normalised data")
        plt.xlabel(r"$Q$ (${\AA}^{-1}$)")
        plt.ylabel(r"$d\sigma/d\Omega$ (barns/sterad/atom)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        print(self.vana["SelfQ0"], vana_inel[1])

        return

    def plot_bs_mc_results(
        self, energias, XS_model_can, xs_can, E0, qvals, single, s_aten, can, total
    ):

        plt.figure(figsize=(6, 4))
        plt.plot(energias, XS_model_can, label="Vanadium model")
        plt.vlines(
            E0, 3, 7, linestyles="--", colors="red", alpha=0.5, label="D4 $E_0$ energy"
        )
        plt.scatter(E0, xs_can, s=50)
        plt.scatter(E0, xs_can)
        plt.tick_params(axis="both", direction="in")
        plt.xlabel("Energy (meV)")
        plt.ylabel("Total Cross section (barns)")
        plt.ylim(3, 6)
        plt.grid()
        plt.legend(loc="best")
        plt.title("Vanadium Model vs D4 $E_0$ energy")
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(6, 4))
        plt.plot(qvals, single[:, 0], lw=2, label="Single scattering", alpha=0.9)
        plt.plot(qvals, single[:, 10], lw=2, label="Multiple scattering", alpha=0.9)
        plt.plot(qvals, s_aten, lw=2, label="unattenuated single scatt.", alpha=0.9)
        plt.plot(qvals, can * 25, lw=2, label="Scattering on can", alpha=0.9)
        plt.plot(qvals, total, lw=2, label="Total scattering", alpha=0.9)
        plt.xlim(0, 20)
        plt.legend(loc="upper right")
        plt.grid()
        plt.tick_params(axis="both", direction="in")
        plt.title("Monte Carlo results")
        plt.xlabel(r"Q ($1/\AA$)")
        plt.ylabel("Intensity (arb. units)")
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(9, 5))
        plt.plot(
            self.vanadiumData_qbin[:, 0],
            self.vanadiumData_qbin[:, 1] * 0.7,
            label="Vanadium",
        )
        plt.plot(
            self.environmentData_qbin[:, 0],
            self.environmentData_qbin[:, 1],
            label="Environment",
        )
        plt.plot(
            self.environmentData_qbin[:, 0],
            (self.vanadiumData_qbin[:, 1] - self.environmentData_qbin[:, 1]) * 0.7,
            label="V-E",
        )
        plt.plot(
            self.van_qbin_sub_mc[:, 0],
            self.van_qbin_sub_mc[:, 1],
            label="Vanadium corrected with ATT and MS",
        )
        plt.xlabel(r"Q ($1/\AA$)")
        plt.ylabel("Intensity (arb. units)")
        plt.title("Compare and results Background Subtraction")
        plt.grid()
        plt.tick_params(axis="both", direction="in")
        plt.legend(loc="best")
        plt.show()

        return
