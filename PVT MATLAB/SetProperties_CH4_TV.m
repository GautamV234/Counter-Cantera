% code for finding values through T V input
function [outputArg1] = SetProperties_CH4_TV(inputArg1,inputArg2)

inpT = inputArg1;
inpV = inputArg2;
% inpT = 105;
% inpV = 0.92;
% assuming input is gotten
subsetA = xlsread('CH4.xlsx', 'satCH4_Tsat');
Pressure_Vector = subsetA(:,2);
VG_Vector = subsetA(:,5);
VF_Vector = subsetA(:,3);
Temperature_Vector = subsetA(:,1);
[rowCount,colCount] = size(VG_Vector);

% error Handling
maxTemperature = Temperature_Vector(rowCount ,1);
% add logic later for checking min from other table also
minTemperature = Temperature_Vector(1 , 1);
%maxVolume = VG_Vector(1,1);
minVolume = VF_Vector(1,1);
% error hanndling regarding inpV > maxVolume to be done in the superHeated Function
if inpT < minTemperature || inpT> maxTemperature || inpV<minVolume
    % wrong input
    fprintf("The input values are out of range. Try different values");
    return
end



% defining output variables
p = -1;
v=-1;
T=-1;
u=-1;
h=-1;
s=-1; 
x=-1;

% interpolating to get the volume val
%error handling

volG = interp1(Temperature_Vector,VG_Vector,inpT);
P_sat = interp1(Temperature_Vector,Pressure_Vector,inpT);
% error handling
% checking if vol>vg
if inpV > volG 
%     fprintf("Must Change Table");
    subsetB = xlsread('CH4.xlsx', 'supHeatCH4');
    Pressure_SuperHeated_Vector = subsetB(:,1);
    Temperature_SuperHeated_Vector =subsetB(:,2);
    V_SuperHeated_Vector = subsetB(:,3);
    U_SuperHeated_Vector = subsetB(:,4);
    H_SuperHeated_Vector = subsetB(:,5);
    S_SuperHeated_Vector = subsetB(:,6);
    % error Handling for checking if inpV> MaxV_SuperHeated
    [rowCount_SH , colCount_SH] = size(V_SuperHeated_Vector);
    rowCount_SH;
%     min_Difference1 = 4480000;
%     min_Difference2 = 4480000;
    lower_Temperature = 95;
%     checking if pressure value exists in the table
    temp_Temperature = lower_Temperature; % temporary Temperature lol
    

% getting all pressure and volume values below Psat of given Temperature
pressure_final_vector = [];
volume_final_vector = [];
pressure_final_vector(end+1) = P_sat;
volume_final_vector(end+1) = volG;
volume_temporary_vector = zeros(10,1);
temperature_temporary_vector = zeros(10,1);
% fprintf("instantiated\n");
% adding pressure vals
i=1;
while i<=340 
    p1 = Pressure_SuperHeated_Vector(i,1);
    % write custom code for checking if the pressure already exists in pressure_final_vector
    len = length(pressure_final_vector);
    l1 = 1;
    flagTemp = 0;
    while l1<=len
        if pressure_final_vector(1,l1)==p1
            flagTemp = 1;
            break;
        end
        l1 = l1+1;
    end
    if flagTemp==0
        if p1<P_sat
             pressure_final_vector(end+1) = p1;
        end
    end
    i=i+1;
end

% now finding other parameters correspondingly

lenOfPressureFinalVector = length(pressure_final_vector);
numOfPressures = lenOfPressureFinalVector;
tempNumOfPressures = numOfPressures;
j=1;
while tempNumOfPressures >0
    currentPressure = pressure_final_vector(1,j);
    % now find all values of T and V corresponding to this Pressure
    k=1;
    while k<=340
        if Pressure_SuperHeated_Vector(k,1)== currentPressure
%             fprintf("Value of pressure \n")
            Pressure_SuperHeated_Vector(k,1);
            tempoindex = k;
            for z = 1:10
                
                volume_temporary_vector(z,1) = V_SuperHeated_Vector(tempoindex,1);
                temperature_temporary_vector(z,1) = Temperature_SuperHeated_Vector(tempoindex,1);
                tempoindex = tempoindex + 1;
            end
%             volume_temporary_vector(end+1) =  V_SuperHeated_Vector(k,1);
%             temperature_temporary_vector(end+1) =  Temperature_SuperHeated_Vector(k,1);
        end
        k = k+10;
    end
%     volumeTempVecTranspose = sort(volume_temporary_vector);
%     temperatureTempVecTranspose = sort(temperature_temporary_vector);
%     volumeTempVecTranspose
%     temperatureTempVecTranspose
    dupliFlag = 0;
    for n=1:length(volume_temporary_vector)-1
        if volume_temporary_vector(n,1)== volume_temporary_vector(n+1,1) || temperature_temporary_vector(n,1)== temperature_temporary_vector(n+1,1)
%             fprintf("Same elem exits\n");
            dupliFlag = 1;
        end
    end
    if dupliFlag==0
        sz = size(volume_temporary_vector);
        if length(volume_temporary_vector)~=0
            volRecieved = interp1(temperature_temporary_vector,volume_temporary_vector,inpT);
            volume_final_vector(end+1) = interp1(temperature_temporary_vector,volume_temporary_vector,inpT);
        end
    end
   

    tempNumOfPressures = tempNumOfPressures - 1;
    j = j+1;
end
p = interp1(sort(volume_final_vector' , 'descend'),sort(pressure_final_vector'),inpV);

% code for other values u,h,s to be calculated like P,V

lower_Pressure = 20000;
temp_Pressure = lower_Pressure;
u=1;
flag = 0;
while u<340
    if temp_Pressure == p
        flag = 1;
        break;
    end
    u=u+1;
    temp_Pressure = Pressure_SuperHeated_Vector(u,1);%% incrementing temp_Pressure
end

if flag==1
    %   add new code here
    %   add code for when pressure is present in the table
%     fprintf("\nPressure value lies in the table!");
    temp_SH_Vector = [];
    volume_SH_Vector = [];
    u_SH_Vector = [];
    hech_SH_Vector =[];
    s_SH_Vector = [];
    for o=u:u+10
% tempTemp1(end+1)=Temperature_SuperHeated_Vector(k,1);
        %temp_SH_Vector(end+1) = Temperature_SuperHeated_Vector(o,1);
        volume_SH_Vector(end+1) = V_SuperHeated_Vector(o,1);
        u_SH_Vector(end+1) = U_SuperHeated_Vector(o,1);
        hech_SH_Vector(end+1) = H_SuperHeated_Vector(o,1);
        s_SH_Vector(end+1) = S_SuperHeated_Vector(o,1);
    end
  
    v= inpV;
    T= inpT;
    u= interp1(volume_SH_Vector', u_SH_Vector', inpV);
    h=interp1(volume_SH_Vector', hech_SH_Vector', inpV);
    s= interp1(volume_SH_Vector', s_SH_Vector', inpV);
    x= 1;
   
% code done for this part        
end

if flag==0
%      fprintf("\nPressure is intermediate");
    % add code for when pressure is intermediate between 2 pressures of the table
    lowerBoundPressure = 20000;
    upperBoundPressure = 20000;
    a=0;
    while(lowerBoundPressure<p && a<340)
        a = a+1;
        lowerBoundPressure = Pressure_SuperHeated_Vector(a,1);
    end
    lowerBoundPressure = Pressure_SuperHeated_Vector(a-1,1);
    upperBoundPressure = Pressure_SuperHeated_Vector(a,1);
    
    % found both lower and upper bound pressurers 
    %tempTemp1 =[];
    tempVol1 = [];
    tempU1 = [];
    tempH1 = [];
    tempS1 = [];
    %tempTemp2 =[];
    tempVol2 = [];
    tempU2 = [];
    tempH2 = [];
    
    tempS2 = [];
    for k=a-10:a-1
        %tempTemp1(end+1)=Temperature_SuperHeated_Vector(k,1);
        tempVol1(end+1)=V_SuperHeated_Vector(k,1);
        tempU1(end+1) = U_SuperHeated_Vector(k,1);
        tempH1(end+1) = H_SuperHeated_Vector(k,1);
        tempS1(end+1) = S_SuperHeated_Vector(k,1);
      
    end
    
    %lowerTemp = tempTemp1';
    lowerVol = tempVol1';
    lowerEnergy = tempU1';
    lowerEnthalpy = tempH1';
    lowerEntropy = tempS1';
    
    for k=a:a+9
        %tempTemp2(end+1)=Temperature_SuperHeated_Vector(k,1);
        tempVol2(end+1)=V_SuperHeated_Vector(k,1);
        tempU2(end+1) = U_SuperHeated_Vector(k,1);
        tempH2(end+1) = H_SuperHeated_Vector(k,1);
        tempS2(end+1) = S_SuperHeated_Vector(k,1);
    end
    %upperTemp = tempTemp2';
    upperVol = tempVol2';
    upperEnergy = tempU2';
    upperEnthalpy = tempH2';
    upperEntropy = tempS2';
   

    


[~, index] = sort(lowerVol);
% Internal Energy
F = griddedInterpolant(lowerVol(index), lowerEnergy(index), 'makima', 'makima');
intEnergy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEnergy(index), 'makima', 'makima');
intEnergy_ou = F(inpV);  % change Notice

IntEnergy = intEnergy_ou + (intEnergy_ou -intEnergy_ol)*(p-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);

% Enthalpy
F = griddedInterpolant(lowerVol(index), lowerEnthalpy(index), 'makima', 'makima');
enthalpy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEnthalpy(index), 'makima', 'makima');
enthalpy_ou = F(inpV); %change Notice


Enthalpy = enthalpy_ou + (enthalpy_ou -enthalpy_ol)*(p-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);

% Entropy 
F = griddedInterpolant(lowerVol(index), lowerEntropy(index), 'makima', 'makima');
entropy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEntropy(index), 'makima', 'makima');
entropy_ou = F(inpV); %change Notice

Entropy = entropy_ou + (entropy_ou -entropy_ol)*(p-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);

%     t1 = interp1(lowerVol , lowerTemp,'linear','extrap' , inpV)
%     t2 = interp1(upperVol , upperTemp , inpV,'linear','extrap' , inpV)
%     testTemp = interp1(Pressure_SuperHeated_Vector,Temperature_SuperHeated_Vector,inpP)
% custom method for linear interpolation
    %vMax = V_SuperHeated_Vector(10,1);
%     Temperature_Temporary_Vector = rand(10,1);  
%assigning all values to the returning variables
T = inpT;
u = IntEnergy;
h = Enthalpy;
s = Entropy;

v = inpV;
x = 1;

end

end
if inpV <= volG
%     fprintf("Right Table for saturated region");
    VF_Vector = subsetA(:,3);
    UF_Vector = subsetA(:,6);
    UG_Vector = subsetA(:,8);
    HF_Vector = subsetA(:,9);
    HG_Vector = subsetA(:,11);
    SF_Vector = subsetA(:,12);
    SG_Vector = subsetA(:,14);
    % adding main logic
    % interpolating to get all values
    %finding Quality Factor (x)
    volF = interp1(Temperature_Vector, VF_Vector,inpT);
    x = (inpV - volF)/(volG - volF);
   
    uF = interp1(Temperature_Vector,UF_Vector,inpT); 
    uG = interp1(Temperature_Vector ,UG_Vector,inpT);
    hF = interp1(Temperature_Vector ,HF_Vector,inpT);
    hG = interp1(Temperature_Vector ,HG_Vector,inpT);
    sF = interp1(Temperature_Vector ,SF_Vector,inpT);
    sG = interp1(Temperature_Vector ,SG_Vector,inpT);
    
    % finding all required vals
    T = inpT;
    u = x*uG + (1-x)*uF;
    h = x*hG + (1-x)*hF;
    s = x*sG + (1-x)*sF;
    p = P_sat;
    v = inpV;
end

% if T>228.666
%     fprintf("Input out of range");
% end
% if T<=228.666
    outputArg1 = [p;v;T;u;h;s;x];
    disp(['Pressure = ', num2str(p),' Pa']);
    disp(['Volume= ', num2str(v),' m^3/kg']);
    disp(['Temperature = ', num2str(T),' K']);
    disp(['Internal Energy = ', num2str(u),' J/kg']);
    disp(['Enthalpy = ', num2str(h),' J/kg']);
    disp(['Entropy = ', num2str(s),' J/kg-K']);
    disp(['Vapour Fraction = ', num2str(x)]);
end
