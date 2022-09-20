% nessma gaber ibrahim kandil sec 4 bn 40 
% resonant frequencies in RWG
clc
clear 
close all

% User input 
cond = 'Enter the conductivity of the medium:';
cond = input(cond); 
eps_r = 'Enter the relative permitivity of the medium:';
eps_r = input(eps_r); 
loss_tang = 'Enter the loss tangent of the medium:';
loss_tang = input(loss_tang); 
% Dimensions Assume b < d < a.
a = 'Enter the 1st dimension a in meters:';
a = input(a); 
b = 'Enter the 2nd dimension b in meters:';
b = input(b); 
d = 'Enter the 3rd dimension d in meters:';
d = input(d);



% Resonant Frequencies.
% create an two arrays to save frequencies values and the corresponding mode. 
if (b < d) && (d < a)    
    c = 3 * 10^8;
    modes_arr = [];
    f_arr = []; 
    for m = 0 : 2
        for n = 0 : 2
            for l = 0 : 2
            % condition to extract the modes that doesn’t occur. 
                if ~(m == 0 && n == 0) && ~(m == 0 && l == 0) && ~(l == 0 && n == 0)
                    f_r = (c/(2* sqrt(eps_r))) * sqrt((m/a)^2 +(n/b)^2 +(l/d)^2);
                    f_arr = [f_arr ;f_r];
                    modes_arr = [modes_arr ;[m n l]]; 
                end
            end
        end
    end 
    % get the min frequency from the array and its position.
    [M1,I1] = min(f_arr);
    first_res_freq = M1;
    first_res_mode = modes_arr(I1,:);
    % delete the min frequency from the array and its position. 
    f_arr = f_arr([1:(I1-1),(I1+1):end]);
    modes_arr = modes_arr([1:(I1-1),(I1+1):end],:); 
    % get the second min frequency from the array and its position.
    [M2,I2] = min(f_arr);
    second_res_freq = M2;
    second_res_mode = modes_arr(I2,:);


    % Quality Factor.
    Qd = 1/loss_tang;

    mu0 = 4*pi*10^-7;
    eps0 = 8.854*10^-12;

    Rs1 = sqrt((pi*mu0*first_res_freq)/cond);
    K1 = (2*pi*first_res_freq*sqrt(eps_r))/c;

    Rs2 = sqrt((pi*mu0*second_res_freq)/cond);
    K2 = (2*pi*second_res_freq*sqrt(eps_r))/c ;

    eta = (120*pi)/sqrt(eps_r);

    
    Qc101 = (((K1*a*d)^3)*b*eta)/((2*(pi^2)*Rs1)*(2*(a^3)*b+ 2*(d^3)*b+ (a^3)*d+ (d^3)*a));

    Q101 = (Qc101*Qd)/(Qc101+Qd);
    % Fisrt Resonance BandWidth.
    first_res_BW =(2*pi*first_res_freq)/Q101;
        
    if second_res_mode(1)== 2 && second_res_mode(2)== 0 &&second_res_mode(3)== 1
        Qc201 = (((K2*a*d)^3)*b*eta)/((2*(pi^2)*Rs2)*(2*(a^3)*b+ 8*(d^3)*b+ (a^3)*d+ 4*(d^3)*a));
        Q201 = (Qc201*Qd)/(Qc201+Qd);
        % Second Resonance BandWidth.
        second_res_BW =(2*pi*second_res_freq)/Q201;
    elseif second_res_mode(1)== 1 && second_res_mode(2)== 1 &&second_res_mode(3)== 0
        Qc110=(eta*pi^2*a^3* b^3 *d*((1/a)^2+(1/b)^2 )^2)/(2*Rs2*K2*(2*b^3* d+ 2*a^3 *d+a^3 * b+a*b^3));    
        Q110 = (Qc110*Qd)/(Qc110+Qd);
        % Second Resonance BandWidth.
        second_res_BW =(2*pi*second_res_freq)/Q110;
        
    end
    X = sprintf('\nFirst Resonance mode = %d %d %d',first_res_mode);
    disp(X);
    X = sprintf('First Resonance Frequency in HZ= %d',first_res_freq);
    disp(X);
    X = sprintf('First Resonance Frequency dielectric quality factor = %d',Qd);
    disp(X);
    X = sprintf('First Resonance Frequency conductor quality factor = %d',Qc101);
    disp(X); 
    X = sprintf('First Resonance Frequency total quality factor = %d',Q101);
    disp(X);
    X = sprintf('First Resonance Frequency bandwidth in radians= %d',first_res_BW);
    disp(X);
    
    X = sprintf('\nSecond Resonance mode = %d %d %d',second_res_mode);
    disp(X);
    X = sprintf('Second Resonance Frequency in HZ= %d',second_res_freq);
    disp(X);
    X = sprintf('Second Resonance Frequency dielectric quality factor = %d',Qd);
    disp(X);
    if second_res_mode(1)== 2 && second_res_mode(2)== 0 &&second_res_mode(3)== 1
        X = sprintf('Second Resonance Frequency conductor quality factor = %d',Qc201);
        disp(X);
        X = sprintf('Second Resonance Frequency total quality factor = %d',Q201);
        disp(X);
        X = sprintf('Second Resonance Frequency bandwidth in radians= %d',second_res_BW);
        disp(X);
    elseif second_res_mode(1)== 1 && second_res_mode(2)== 1 &&second_res_mode(3)== 0
        X = sprintf('Second Resonance Frequency conductor quality factor = %d',Qc110);
        disp(X);
        X = sprintf('Second Resonance Frequency total quality factor = %d',Q110);
        disp(X);
        X = sprintf('Second Resonance Frequency bandwidth in radians= %d',second_res_BW);
        disp(X);
    end

else
    disp('Invalid Dimensions');
    
end













