% function for P,V
function [outputArg1] = SetProperties_CH4_PV(inputArg1,inputArg2)

inpP = inputArg1;
inpV = inputArg2;


% inpP = 732817;
% inpV = 0.125;

% assuming input is gotten

subsetA = xlsread('CH4.xlsx', 'satCH4_Psat');
Pressure_Vector = subsetA(:,1);
VG_Vector = subsetA(:,5);
VF_Vector = subsetA(:,3);
Temperature_Vector = subsetA(:,2);
[rowCount,colCount] = size(VG_Vector);

% error Handling
maxPressure = Pressure_Vector(rowCount ,1);
% add logic later for checking min from other table also
minPressure = Pressure_Vector(1 , 1);
%maxVolume = VG_Vector(1,1);
minVolume = VF_Vector(1,1);
% error hanndling regarding inpV > maxVolume to be done in the superHeated Function
if inpP < minPressure || inpP> maxPressure || inpV<minVolume
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
volG = interp1(Pressure_Vector,VG_Vector,inpP);
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
%     min_Difference1 = 4480000;
%     min_Difference2 = 4480000;
    lower_Pressure = 20000;
%     checking if pressure value exists in the table
    temp_Pressure = lower_Pressure; % temporay pressure lol
    

% checking if pressure is there
u=1;
flag = 0;
while u<340
    if temp_Pressure == inpP
        flag = 1;
        break;
    end
    u=u+1;
    temp_Pressure = Pressure_SuperHeated_Vector(u,1);%% incrementing temp_Pressure
end

if flag==1
    %   add code for when pressure is present in the table
%     fprintf("\nPressure value lies in the table!");
    temp_SH_Vector = [];
    volume_SH_Vector = [];
    u_SH_Vector = [];
    hech_SH_Vector =[];
    s_SH_Vector = [];
    for o=u:u+10
% tempTemp1(end+1)=Temperature_SuperHeated_Vector(k,1);
        temp_SH_Vector(end+1) = Temperature_SuperHeated_Vector(o,1);
        volume_SH_Vector(end+1) = V_SuperHeated_Vector(o,1);
        u_SH_Vector(end+1) = U_SuperHeated_Vector(o,1);
        hech_SH_Vector(end+1) = H_SuperHeated_Vector(o,1);
        s_SH_Vector(end+1) = S_SuperHeated_Vector(o,1);
    end
    p = inpP;
    v= inpV;
    T= interp1(volume_SH_Vector', temp_SH_Vector', inpV);
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
    while(lowerBoundPressure<inpP && a<340)
        a = a+1;
        lowerBoundPressure = Pressure_SuperHeated_Vector(a,1);
    end
    lowerBoundPressure = Pressure_SuperHeated_Vector(a-1,1);
    upperBoundPressure = Pressure_SuperHeated_Vector(a,1);
    
    % found both lower and upper bound pressurers 
    tempTemp1 =[];
    tempVol1 = [];
    tempU1 = [];
    tempH1 = [];
    tempS1 = [];
    tempTemp2 =[];
    tempVol2 = [];
    tempU2 = [];
    tempH2 = [];
    
    tempS2 = [];
    for k=a-10:a-1
        tempTemp1(end+1)=Temperature_SuperHeated_Vector(k,1);
        tempVol1(end+1)=V_SuperHeated_Vector(k,1);
        tempU1(end+1) = U_SuperHeated_Vector(k,1);
        tempH1(end+1) = H_SuperHeated_Vector(k,1);
        tempS1(end+1) = S_SuperHeated_Vector(k,1);
      
    end
    
    lowerTemp = tempTemp1';
    lowerVol = tempVol1';
    lowerEnergy = tempU1';
    lowerEnthalpy = tempH1';
    lowerEntropy = tempS1';
    
    for k=a:a+9
        tempTemp2(end+1)=Temperature_SuperHeated_Vector(k,1);
        tempVol2(end+1)=V_SuperHeated_Vector(k,1);
        tempU2(end+1) = U_SuperHeated_Vector(k,1);
        tempH2(end+1) = H_SuperHeated_Vector(k,1);
        tempS2(end+1) = S_SuperHeated_Vector(k,1);
    end
    upperTemp = tempTemp2';
    upperVol = tempVol2';
    upperEnergy = tempU2';
    upperEnthalpy = tempH2';
    upperEntropy = tempS2';
   

    

% Temperature
[~, index] = sort(lowerTemp);
F = griddedInterpolant(lowerVol(index), lowerTemp(index), 'makima', 'makima');
Temperature_ol = F(inpV);
F = griddedInterpolant(upperVol(index), upperTemp(index), 'makima', 'makima');

Temperature_ou = F(inpV);
T_out = Temperature_ol + (Temperature_ou-Temperature_ol)*(inpP-lowerBoundPressure)/(upperBoundPressure-lowerBoundPressure);

% finding all other relevant paramters in a similar way

% Internal Energy
F = griddedInterpolant(lowerVol(index), lowerEnergy(index), 'makima', 'makima');
intEnergy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEnergy(index), 'makima', 'makima');
intEnergy_ou = F(inpV);  % change Notice

IntEnergy = intEnergy_ou + (intEnergy_ou -intEnergy_ol)*(inpP-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);

% Enthalpy
F = griddedInterpolant(lowerVol(index), lowerEnthalpy(index), 'makima', 'makima');
enthalpy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEnthalpy(index), 'makima', 'makima');
enthalpy_ou = F(inpV); %change Notice


Enthalpy = enthalpy_ou + (enthalpy_ou -enthalpy_ol)*(inpP-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);

% Entropy 
F = griddedInterpolant(lowerVol(index), lowerEntropy(index), 'makima', 'makima');
entropy_ol = F(inpV);

F = griddedInterpolant(upperVol(index), upperEntropy(index), 'makima', 'makima');
entropy_ou = F(inpV); %change Notice

Entropy = entropy_ou + (entropy_ou -entropy_ol)*(inpP-lowerBoundPressure)/(upperBoundPressure - lowerBoundPressure);
  
%assigning all values to the returning variables
T = T_out;
u = IntEnergy;
h = Enthalpy;
s = Entropy;
p = inpP;
v = inpV;
x = 1;

end

end
if inpV <= volG
%     fprintf("Right Table Hooray");
    %VF_Vector = subsetA(:,3);
    UF_Vector = subsetA(:,6);
    UG_Vector = subsetA(:,8);
    HF_Vector = subsetA(:,9);
    HG_Vector = subsetA(:,11);
    SF_Vector = subsetA(:,12);
    SG_Vector = subsetA(:,14);
    % adding main logic
    % interpolating to get all values
    %finding Quality Factor (x)
    volF = interp1(Pressure_Vector, VF_Vector,inpP);
    x = (inpV - volF)/(volG - volF);
   
    uF = interp1(Pressure_Vector,UF_Vector,inpP); 
    uG = interp1(Pressure_Vector ,UG_Vector,inpP);
    hF = interp1(Pressure_Vector ,HF_Vector,inpP);
    hG = interp1(Pressure_Vector ,HG_Vector,inpP);
    sF = interp1(Pressure_Vector ,SF_Vector,inpP);
    sG = interp1(Pressure_Vector ,SG_Vector,inpP);
    
    % finding all required vals
    T = interp1(Pressure_Vector,Temperature_Vector,inpP);
    u = x*uG + (1-x)*uF;
    h = x*hG + (1-x)*hF;
    s = x*sG + (1-x)*sF;
    p = inpP;
    v = inpV;
end

if T<=228.666
    outputArg1 = [p;v;T;u;h;s;x];
    disp(['Pressure = ', num2str(p),' Pa']);
    disp(['Volume= ', num2str(v),' m^3/kg']);
    disp(['Temperature = ', num2str(T),' K']);
    disp(['Internal Energy = ', num2str(u),' J/kg']);
    disp(['Enthalpy = ', num2str(h),' J/kg']);
    disp(['Entropy = ', num2str(s),' J/kg-K']);
    disp(['Vapour Fraction = ', num2str(x)]);
else
    fprintf("The input values are out of range. Try different values.");
end


end

     
     
    



