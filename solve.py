# from gurobipy import Model
# from gurobipy import GRB
# from gurobipy import LinExpr
# from gurobipy import tuplelist
import math
import numpy as np

desiredSpeed = 60; # m/s
acceleration_distance = 700; # m
# pod_mass = 200; # kg, w/ LIM, estimated
# stator_tooth_height = 0.05; # m

n_poles_choices = [1,2,3,4,5]

# stator_tooth_thickness = 0.008; # m, tooth thickness
# stator_yoke_height = 0.008; # m, yoke thickness
stator_width = 0.127; # m
# airgap = 0.005; # m
# voltage = 105; # VAC, line-to-line
packEff = 0.80;
dutyCycle = 0.1;
slip = 0.05;
powerFactor = 0.57;
forceFactor = 0.25;
coilSurfAreaFactor = 1;
aspect_ratio = 4;

operating_temp = 100 + 273.15; # Kelvin
envTemp = 30 + 273.15; # Kelvin

magPermVac = 0.000001256637061; # H/m
secondary_thickness = 0.009525; # m
conductorRes = 0.0000000168; # Ohm-m
surfResAl = 0.0000000282; # Ohm-m
alpha = 0.00404; # thermal coefficient of res. for conductor
statorYoungs = 200000000000; # Pa
conductorDens = 8960; # kg/m3 for Copper
stator_density = 7800; # kg/m3 for Steel/Iron
insulator_density = 1420; # kg/m3 for PEI
epsilon = 0.9; # surface emmistivity of conductor
lengthCutoff = 1.5; # m

conductor_diameters = [0.00064262, 0.000724, 0.000813, 0.000912, 0.001024, 0.001151, 0.00129, 0.00145, 0.001628, 0.001829, 0.002052, 0.002304, 0.002588, 0.002906, 0.003264]; # 22 to 8awg, every 1awg
wire_diameters = [0.00067564, 0.000757, 0.000851, 0.000947, 0.001062, 0.001191, 0.001331, 0.001491, 0.001674, 0.001872, 0.002096, 0.00235, 0.002634, 0.002951, 0.003307];
stefanBoltzmann = 0.000000056703;


#
stator_tooth_thickness_choices = np.linspace(0.005,0.02,10)
stator_tooth_height_choices = np.linspace(0.06,0.08,5)
# stator_yoke_height_choices = np.linspace(0.005,0.02,10)
# airgap_choices = np.linspace(0.001,0.005,5)
voltage_choices = np.linspace(100,110,10)

# n_poles_choices = range(2,8,2)
# n_long_turns_choices = range(10,21)

# stator_tooth_height_choices = [0.06]
# stator_tooth_thickness_choices = [0.00833]
stator_yoke_height_choices = [0.001]
airgap_choices = [0.001]
n_poles_choices =[4,6,8]
n_long_turns_choices = [10, 12, 14,16,18,20]
voltage_choices = [105.55]

mass_weight = 0.2
top_speed_weight =  0.8

# --------------------------- SOLVER ---------------------------

# reqForce = pod_mass*((desiredSpeed**2)/(2*acceleration_distance));

design_score_old = -1000

counter = 0

for stator_tooth_height in stator_tooth_height_choices:
    for stator_tooth_thickness in stator_tooth_thickness_choices:
        for stator_yoke_height in stator_yoke_height_choices:
            for airgap in airgap_choices:
                for voltage in voltage_choices:
                    for wire_index, wire_diameter in enumerate(wire_diameters):
                        # stator_slot_pitch = (wire_diameter/packEff) + stator_tooth_thickness;
                        # init_stator_slot_pitch = stator_slot_pitch;
                        for n_long_turns in n_long_turns_choices:
                            stator_slot_pitch = n_long_turns*(wire_diameter/packEff) + stator_tooth_thickness
                            for n_poles in n_poles_choices:
                                conductor_cross_area = math.pi * (conductor_diameters[wire_index]/2)**2;
                                n_vertical_turns = math.floor((stator_tooth_height/(wire_diameter/packEff))/3);

                                if stator_tooth_height <((wire_diameter/packEff)*12):
                                    stator_tooth_height = ((wire_diameter/packEff)*12)
                                    print("stator tooth height adjusted")

                                n_turns_per_coil = n_long_turns * n_vertical_turns;
                                turnLength = (2*stator_width) + (4*stator_slot_pitch);
                                coil_wire_length = turnLength * n_turns_per_coil;

                                rtResistance = ((conductorRes * coil_wire_length)/conductor_cross_area) * n_poles;
                                envResistance = (rtResistance * (1+(alpha*(envTemp-293.15))));
                                opResistance = (rtResistance * (1+(alpha*(operating_temp-293.15))));
                                coilSurfArea = ((turnLength*n_vertical_turns*math.pi*wire_diameter) + (turnLength*n_long_turns*math.pi*wire_diameter))*coilSurfAreaFactor;
                                calculated_operating_temp = ((((((voltage**2/opResistance)*dutyCycle)/(epsilon*coilSurfArea*n_poles))+(stefanBoltzmann*envTemp**4))/stefanBoltzmann)**0.25);

                                stator_length = (5+3*(n_poles-1)) * stator_slot_pitch + stator_tooth_thickness;

                                if (stator_length>=lengthCutoff) or (calculated_operating_temp >= operating_temp) or n_long_turns >= (aspect_ratio*n_vertical_turns):
                                    # print("rejected")
                                    continue
                                else:
                                    toothGap = stator_slot_pitch - stator_tooth_thickness;
                                    windingFactor = 1;
                                    magnetic_air_gap = airgap + secondary_thickness;
                                    gamma = (4/math.pi) * (((toothGap/(2*magnetic_air_gap)) * math.atan(toothGap/(2*magnetic_air_gap))) - math.log((1+ (toothGap/(2*magnetic_air_gap))**2)**0.5));
                                    k_carter = stator_slot_pitch/(stator_slot_pitch - (gamma*magnetic_air_gap));
                                    effective_air_gap = k_carter * magnetic_air_gap;
                                    eq_stator_width = stator_width + magnetic_air_gap;
                                    tau = stator_length/n_poles;
                                    fMax = math.floor(desiredSpeed/(2*tau*(1-slip)));

                                    Xm = (24*magPermVac*math.pi*fMax*eq_stator_width*windingFactor*((n_poles*n_turns_per_coil)**2)*tau)/((math.pi**2)*n_poles*effective_air_gap); #PF = per frequency
                                    goodness_factor = (2*magPermVac*fMax*tau**2)/(math.pi*(surfResAl/secondary_thickness)*effective_air_gap);
                                    rSubTwo = Xm/goodness_factor;
                                    rotorPhaseCurrent = (voltage/opResistance)/((1+ (1/((slip*goodness_factor)**2)))**0.5);
                                    syncVelo = 2*fMax*tau;

                                    force = ((3*rotorPhaseCurrent**2*rSubTwo)/(syncVelo*slip)) * forceFactor;

                                    # calculation of objectives
                                    mass = n_poles * 3 * coil_wire_length * conductor_cross_area * conductorDens; # add mass of conductor
                                    mass += (n_poles * 3 * coil_wire_length * ((math.pi * (wire_diameter / 2) ** 2) - conductor_cross_area) * insulator_density);  # add mass of insulation
                                    mass += (stator_width * (stator_tooth_height + stator_yoke_height) * (stator_tooth_thickness + stator_slot_pitch * (5 + 3 * (n_poles - 1))) * stator_density); # add stator as giant block
                                    mass +- ((5 + 3 * (n_poles - 1)) * stator_density * stator_width * (stator_slot_pitch - stator_tooth_thickness) * stator_tooth_height); # subtract mass of slots

                                    pod_mass = 60 + mass
                                    top_speed = math.sqrt(2*force*acceleration_distance/pod_mass)

                                    design_score_new = top_speed_weight*top_speed - mass_weight * mass

                                    if design_score_new > design_score_old:
                                        best_design_specs = {
                                        "params":
                                            {
                                            "stator_tooth_height":stator_tooth_height,
                                            "stator_tooth_thickness": stator_tooth_thickness,
                                            "stator_yoke_height":stator_yoke_height,
                                            "air gap": airgap,
                                            "voltage": voltage,
                                            "n_poles": n_poles,
                                            "n_long_turns": n_long_turns,
                                            "n_vertical_turns": n_vertical_turns,
                                            "operating_temp": operating_temp,
                                            "turns_per_coil": n_turns_per_coil,
                                            "wire diameter": wire_diameter,
                                            "wire gauge": 22 - (wire_index),
                                            # "stator_slot_pitch": stator_slot_pitch,

                                            "tooth_pitch":stator_slot_pitch,
                                            "force": force
                                            },
                                        "objectives":
                                            {
                                            "lim_mass": mass,
                                            "speed":top_speed,
                                            "mass weight": mass_weight,
                                            "speed weight": top_speed_weight
                                            }
                                        }

                                        design_score_old = design_score_new
                                        # print("top speed: {}, mass: {}".format(top_speed, mass))
                                        print(best_design_specs)
                                    #
                                    # if :
                                    #     n_poles = n_poles +1;
                                    #     n_long_turns = 1;
                                    #     stator_slot_pitch = init_stator_slot_pitch;
                                    #
                                    # else:
                                    #     n_long_turns = n_long_turns + 1;
                                    #     stator_slot_pitch = stator_slot_pitch + (wire_diameter/packEff);
                                    counter += 1
                                    # print(counter)

print("Number of combinations:", counter)
print(best_design_specs)
