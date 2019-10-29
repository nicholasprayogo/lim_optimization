%Single-Sided Linear Induction Motor Design
%Ryerson International Hyperloop
%Created: October 26, 2019 - Artin Sarkezians
%Last Modified: October 27, 2019 - Artin Sarkezians

%--------------------------- SOLVER INPUTS ---------------------------

desiredSpeed = 60; %m/s
accelDist = 700; %m
podMass = 100; %kg, w/ LIM, estimated
statorH = 0.05; %m
statorT = 0.008; %m, tooth thickness
statorB = 0.008; %m, yoke thickness
statorWdt = 0.127; %m
airGap = 0.005; %m
voltage = 105; %VAC, line-to-line
packEff = 0.80;
dutyCycle = 0.1;
slip = 0.05;
powerFactor = 0.57;
forceFactor = 0.25;
coilSurfAreaFactor = 1;
aspectRatio = 4;
opTemp = 100 + 273.15; %Kelvin
envTemp = 30 + 273.15; %Kelvin

magPermVac = 0.000001256637061; %H/m
secondaryThk = 0.009525; %m
conductorRes = 0.0000000168; %Ohm-m
surfResAl = 0.0000000282; %Ohm-m
alpha = 0.00404; %thermal coefficient of res. for conductor
statorYoungs = 200000000000; %Pa
conductorDens = 8960; %kg/m3 for Copper
statorDens = 7800; %kg/m3 for Steel/Iron
insulatorDens = 1420; %kg/m3 for PEI
epsilon = 0.9; %surface emmistivity of conductor
lengthCutoff = 1.5; %m

conductorDia = [0.00064262, 0.000724, 0.000813, 0.000912, 0.001024, 0.001151, 0.00129, 0.00145, 0.001628, 0.001829, 0.002052, 0.002304, 0.002588, 0.002906, 0.003264]; %22 to 8awg, every 1awg
wireDia = [0.00067564, 0.000757, 0.000851, 0.000947, 0.001062, 0.001191, 0.001331, 0.001491, 0.001674, 0.001872, 0.002096, 0.00235, 0.002634, 0.002951, 0.003307];
stefanBoltzmann = 0.000000056703;

%--------------------------- SOLVER ---------------------------

reqForce = podMass*((desiredSpeed^2)/(2*accelDist));

for wireGA=1:15
    longTurns = 1;
    numPoles = 1;
    conductorCA = pi * (conductorDia(wireGA)/2)^2;
    verticalTurns = floor((statorH/(wireDia(wireGA)/packEff))/3);
    statorX = (wireDia(wireGA)/packEff) + statorT;
    initStatorX = statorX;

    %adjust stator tooth height for 4 minimum turns
    if statorH < (wireDia(wireGA)/packEff)*12
       statorH = (wireDia(wireGA)/packEff)*12;
       msgbox('stator tooth height was adjusted');
    end

    while true
        numTurns = longTurns * verticalTurns;
        turnLength = (2*statorWdt) + (4*statorX);
        coilWireLength = turnLength * numTurns;

        rtResistance = ((conductorRes * coilWireLength)/conductorCA) * numPoles;
        envResistance = (rtResistance * (1+(alpha*(envTemp-293.15))));
        opResistance = (rtResistance * (1+(alpha*(opTemp-293.15))));
        coilSurfArea = ((turnLength*verticalTurns*pi*wireDia(wireGA)) + (turnLength*longTurns*pi*wireDia(wireGA)))*coilSurfAreaFactor;
        calcOpTemp = ((((((voltage^2/opResistance)*dutyCycle)/(epsilon*coilSurfArea*numPoles))+(stefanBoltzmann*envTemp^4))/stefanBoltzmann)^0.25);

        statorLength = (5+3*(numPoles-1)) * statorX + statorT;

        toothGap = statorX - statorT;
        windingFactor = 1; %will always be 1 for 3 phases, 1 slot/pole/phase
        magAirGap = airGap + secondaryThk;
        gamma = (4/pi) * (((toothGap/(2*magAirGap)) * atan(toothGap/(2*magAirGap))) - log((1+ (toothGap/(2*magAirGap))^2)^0.5));
        carter = statorX/(statorX - (gamma*magAirGap));
        effAirGap = carter * magAirGap;
        eqStatorWdt = statorWdt + magAirGap;
        tau = statorLength/numPoles;
        fMax = floor(desiredSpeed/(2*tau*(1-slip)));
        xSubM = (24*magPermVac*pi*fMax*eqStatorWdt*windingFactor*((numPoles*numTurns)^2)*tau)/((pi^2)*numPoles*effAirGap); %PF = per frequency
        goodness = (2*magPermVac*fMax*tau^2)/(pi*(surfResAl/secondaryThk)*effAirGap);
        rSubTwo = xSubM/goodness;
        rotorPhaseCurrent = (voltage/opResistance)/((1+ (1/((slip*goodness)^2)))^0.5);
        syncVelo = 2*fMax*tau;
        force = ((3*rotorPhaseCurrent^2*rSubTwo)/(syncVelo*slip)) * forceFactor;

        mass = numPoles * 3 * coilWireLength * conductorCA * conductorDens; %add mass of conductor
        mass = mass + (numPoles * 3 * coilWireLength * ((pi * (wireDia(wireGA) / 2) ^ 2) - conductorCA) * insulatorDens);  %add mass of insulation
        mass = mass + (statorWdt * (statorH + statorB) * (statorT + statorX * (5 + 3 * (numPoles - 1))) * statorDens); %add stator as giant block
        mass = mass - ((5 + 3 * (numPoles - 1)) * statorDens * statorWdt * (statorX - statorT) * statorH); %subtract mass of slots

        if statorLength >= lengthCutoff
            break
        elseif force >= reqForce && calcOpTemp <= opTemp
            break
        else
            if longTurns >= (aspectRatio*verticalTurns)
                numPoles = numPoles +1;
                longTurns = 1;
                statorX = initStatorX;
            else
                longTurns = longTurns + 1;
                statorX = statorX + (wireDia(wireGA)/packEff);
            end
        end
    end

    if force >= reqForce && calcOpTemp <= opTemp
        break
    end
end

fprintf('\nOutput Force: %f N\n', force)
fprintf('Trajectory Required Force: %f N\n', reqForce)
fprintf('0-100km/h Acceleration Time: %f s\n', (27.777*podMass)/force)
fprintf('LIM Power-On Time: %f s\n\n', (desiredSpeed*podMass)/force)

fprintf('Operating Temp: %f C\n', calcOpTemp-273.15)
fprintf('Wire Gauge: %d AWG \n', 22-(wireGA-1))
fprintf('Turns per Coil: %d\n', longTurns*verticalTurns)
fprintf('Height-wise Turns: %d\n', verticalTurns)
fprintf('Length-wise Turns: %d\n', longTurns)
fprintf('Coil Wire Length: %f m\n', coilWireLength)
fprintf('Resistance per Phase at R.T: %f Ohm\n', rtResistance)
fprintf('Pesistance per Phase at O.T: %f Ohm\n', opResistance)
fprintf('Total Number of Coils: %d \n\n', numPoles*3)

fprintf('Number of Poles: %d \n', numPoles)
fprintf('Stator Length: %f m\n', statorLength)
fprintf('Maximum Operating Freq: %d Hz\n', fMax)
fprintf('LIM mass: %f kg\n\n', mass)

%--------------------------- TRAJECTORY ---------------------------

%plot the motor speed at variable frequency/voltage in 0.1s increments
