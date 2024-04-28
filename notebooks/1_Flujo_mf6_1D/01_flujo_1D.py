
import os   # Intefaces con el sistema operativo.
import numpy as np # Manejo de arreglos numéricos multidimensionales
import matplotlib.pyplot as plt # Graficación

# Biblioteca y módulos de flopy
import flopy
from flopy.plot.styles import styles

import vflow

mesh = vflow.MeshDis(
    nrow = 1,    # Number of rows
    ncol = 120,  # Number of columns
    nlay = 1,    # Number of layers
    row_length = 12.0,    # Length of system ($cm$)
    column_length = 0.1,  # Length of system ($cm$)
    top = 1.0,   # Top of the model ($cm$)
    bottom = 0,  # Layer bottom elevation ($cm$)
)

tparams = dict(nper = 1,  # Number of periods
               total_time = 120.0,  # Simulation time ($s$)
               nstp = 1.0,   # Number of time steps
               tsmult = 1.0  # Multiplier for the length of successive time steps.
              )

model_units = {"time": "seconds",
               "length": "centimeters"}

pparams = dict(hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
               specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
               source_concentration = 1.0  # Source concentration (unitless)
              )

os_params = dict(ws = os.getcwd(), # Ruta de donde estamos actualmente
                 mf6_exe = '/home/jovyan/GMMC/WMA/mf6/bin/mf6', # Ejecutable
                 name = 'flow', # Nombre de la simulación
                 head_file = "flow.hds", 
                 budget_file = "flow.bud"
                )

print(' Ruta hacia el ejecutable de MODFLOW 6 : {} '.format(os_params["mf6_exe"]))
print(' Ruta actual : {}'.format(os_params["ws"]))
print(' Nombre de esta simulación: {}'.format(os_params["name"]))

def build_gwf(mesh, tparams, pparams, model_units, os_params):
    sim = flopy.mf6.MFSimulation(sim_name=os_params["name"], 
                                 sim_ws=os_params["ws"], 
                                 exe_name=os_params["mf6_exe"])
    flopy.mf6.ModflowTdis(sim, 
                          nper=tparams["nper"], 
                          perioddata=((tparams["total_time"], 
                                       tparams["nstp"], 
                                       tparams["tsmult"]),), 
                          time_units=model_units["time"])
    
    flopy.mf6.ModflowIms(sim)
    
    gwf = flopy.mf6.ModflowGwf(sim, 
                               modelname=os_params["name"],
                               save_flows=True)
    
    
    flopy.mf6.ModflowGwfdis(gwf,
                            length_units=model_units["length"],
                            nlay=mesh.nlay,
                            nrow=mesh.nrow,
                            ncol=mesh.ncol,
                            delr=mesh.delr,
                            delc=mesh.delc,
                            top=mesh.top,
                            botm=mesh.bottom,
    )
    
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    
    flopy.mf6.ModflowGwfnpf(gwf,
                            save_specific_discharge=True,
                            save_saturation=True,
                            icelltype=0,
                            k=pparams["hydraulic_conductivity"],
    )
    
    flopy.mf6.ModflowGwfchd(gwf, 
                            stress_period_data=[[(0, 0, mesh.ncol - 1), 1.0]]) 
    
    q = pparams["specific_discharge"] * mesh.delc * mesh.delr * mesh.top
    aux = pparams["source_concentration"]
    flopy.mf6.ModflowGwfwel(gwf,
                            stress_period_data=[[(0, 0, 0), q, aux,]],
                            pname="WEL-1",
                            auxiliary=["CONCENTRATION"],
    )
    
    flopy.mf6.ModflowGwfoc(gwf,
                           head_filerecord=os_params["head_file"],
                           budget_filerecord=os_params["budget_file"],
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, gwf

def plot_results(gwf, mesh, os_params):
    # Obtenemos los resultados de la carga hidráulica
    head = flopy.utils.HeadFile(os.path.join(os_params["ws"], os_params["head_file"])).get_data()
    
    # Obtenemos los resultados del BUDGET
    bud  = flopy.utils.CellBudgetFile(os.path.join(os_params["ws"], os_params["budget_file"]),
                                      precision='double')
    # Obtenemos las velocidades
    spdis = bud.get_data(text='DATA-SPDIS')[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
    
    with styles.USGSPlot():
        plt.rcParams['font.family'] = 'DeJavu Sans'
        x, _, _ = mesh.get_coords()
        plt.figure(figsize=(10,3))
        plt.plot(x, head[0, 0], marker=".", ls ="-", mec="blue", mfc="none", markersize="1", label = 'Head')
        plt.xlim(0, 12)
        plt.xticks(ticks=np.linspace(0, mesh.row_length,13))
        plt.xlabel("Distance (cm)")
        plt.ylabel("Normalized Concentration (unitless)")
        plt.legend()
        plt.grid()
        plt.savefig('head.pdf')


sim, gwf = build_gwf(mesh, tparams, pparams, model_units, os_params)

sim.write_simulation()

sim.run_simulation()

plot_results(gwf, mesh, os_params)
    